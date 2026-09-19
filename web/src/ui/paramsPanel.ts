// Parameters left panel (v1: Mill / Media / Lifters / Slurry / Simulation groups, schema-driven
// from src/params/schema.ts, plus a Derived group). See docs/PLAN.md ss4.3. Replaces the earlier
// <dialog>-based modal (ui/paramsModal.ts) with a persistent <aside> panel that sits beside the
// canvas, symmetric with the metrics panel on the right. A Display tab and full validation UX
// beyond field-level range checks land in M5; for now Apply always resets the simulation (see
// schema.ts's v1 simplification note).

import type { ParamsJson } from "../protocol";
import {
  criticalSpeedRpm,
  effectiveMedia,
  fluidParticleCountEstimate,
  percentCriticalOf,
  rpmOf,
  substepDisplacementOverDiameter,
  trueBallCount,
} from "../params/derived";
import { fromDisplayValue, GROUPS, getPath, SCHEMA, toDisplayValue, withPath } from "../params/schema";

// pbf.rs's default fluid lattice resolution (crates/mill-core/src/params.rs
// SimulationParams::default). Not exposed in SCHEMA (out of scope for this panel), so it's
// hardcoded here purely as a static estimate input for the "Est. fluid particles" readout.
const DEFAULT_FLUID_RESOLUTION = 40;

const GROUPS_STORAGE_KEY = "milldynamics.paramsGroups";

const DEFAULT_GROUP_OPEN: Record<string, boolean> = {
  Mill: true,
  Media: true,
  Lifters: false,
  Slurry: false,
  Simulation: false,
  Derived: true,
};

function readGroupOpenPrefs(): Record<string, boolean> {
  try {
    const raw = localStorage.getItem(GROUPS_STORAGE_KEY);
    if (raw === null) return { ...DEFAULT_GROUP_OPEN };
    const parsed = JSON.parse(raw) as Record<string, unknown>;
    const result: Record<string, boolean> = { ...DEFAULT_GROUP_OPEN };
    for (const key of Object.keys(result)) {
      if (typeof parsed[key] === "boolean") result[key] = parsed[key] as boolean;
    }
    return result;
  } catch {
    return { ...DEFAULT_GROUP_OPEN };
  }
}

function writeGroupOpenPrefs(prefs: Record<string, boolean>): void {
  try {
    localStorage.setItem(GROUPS_STORAGE_KEY, JSON.stringify(prefs));
  } catch {
    // ignore (private browsing / disabled storage)
  }
}

export interface ParamsPanel {
  /** Refreshes the form to reflect `params` (called on every worker "ready"). */
  setParams(params: ParamsJson): void;
  /** Surfaces a worker-reported error in the panel's status/error UI. */
  showError(message: string): void;
}

interface FieldFailure {
  path: string;
  message: string;
}

export function createParamsPanel(container: HTMLElement, onApply: (params: ParamsJson) => void): ParamsPanel {
  const groupOpenPrefs = readGroupOpenPrefs();

  const form = document.createElement("form");
  form.className = "params-form";
  form.noValidate = true;
  container.appendChild(form);

  // Sticky header: Apply / Revert / status indicator.
  const header = document.createElement("div");
  header.className = "params-panel-header";

  const applyBtn = document.createElement("button");
  applyBtn.type = "submit";
  applyBtn.textContent = "Apply (resets simulation)";

  const revertBtn = document.createElement("button");
  revertBtn.type = "button";
  revertBtn.textContent = "Revert";
  revertBtn.addEventListener("click", () => setParams(currentParams));

  const status = document.createElement("span");
  status.className = "params-status";
  status.textContent = "Applied";

  header.append(applyBtn, revertBtn, status);
  form.appendChild(header);

  function setStatus(state: "applied" | "edited" | "error"): void {
    status.classList.remove("is-edited", "is-error");
    if (state === "applied") {
      status.textContent = "Applied";
    } else if (state === "edited") {
      status.textContent = "Edited — not applied";
      status.classList.add("is-edited");
    } else {
      status.textContent = "Error";
      status.classList.add("is-error");
    }
  }

  const errorBanner = document.createElement("div");
  errorBanner.className = "params-error";
  errorBanner.hidden = true;
  form.appendChild(errorBanner);

  const inputs = new Map<string, HTMLInputElement | HTMLSelectElement>();
  const rowsByPath = new Map<string, HTMLElement>();
  const detailsByGroup = new Map<string, HTMLDetailsElement>();
  let currentParams: ParamsJson = {};
  let mediaWarn: HTMLParagraphElement | null = null;

  for (const group of GROUPS) {
    const details = document.createElement("details");
    details.className = "params-group";
    details.open = groupOpenPrefs[group] ?? false;
    details.addEventListener("toggle", () => {
      groupOpenPrefs[group] = details.open;
      writeGroupOpenPrefs(groupOpenPrefs);
    });
    detailsByGroup.set(group, details);

    const summary = document.createElement("summary");
    summary.textContent = group;
    details.appendChild(summary);

    for (const field of SCHEMA.filter((f) => f.group === group)) {
      const row = document.createElement("label");
      row.className = "params-row";

      const labelText = document.createElement("span");
      labelText.textContent = field.unit ? `${field.label} (${field.unit})` : field.label;
      row.appendChild(labelText);

      let input: HTMLInputElement | HTMLSelectElement;
      if (field.type === "select") {
        const select = document.createElement("select");
        for (const opt of field.options ?? []) {
          const optionEl = document.createElement("option");
          optionEl.value = opt.value;
          optionEl.textContent = opt.label;
          select.appendChild(optionEl);
        }
        input = select;
      } else if (field.type === "boolean") {
        const checkbox = document.createElement("input");
        checkbox.type = "checkbox";
        input = checkbox;
      } else {
        const numberInput = document.createElement("input");
        numberInput.type = "number";
        if (field.min !== undefined) numberInput.min = String(field.min);
        if (field.max !== undefined) numberInput.max = String(field.max);
        if (field.step !== undefined) numberInput.step = String(field.step);
        input = numberInput;
      }
      input.addEventListener("input", () => {
        setStatus("edited");
        refreshDerived();
      });
      inputs.set(field.path, input);
      rowsByPath.set(field.path, row);
      row.appendChild(input);
      details.appendChild(row);
    }

    if (group === "Media") {
      mediaWarn = document.createElement("p");
      mediaWarn.className = "params-warn";
      mediaWarn.hidden = true;
      details.appendChild(mediaWarn);
    }

    form.appendChild(details);
  }

  const derivedDetails = document.createElement("details");
  derivedDetails.className = "params-group";
  derivedDetails.open = groupOpenPrefs.Derived ?? true;
  derivedDetails.addEventListener("toggle", () => {
    groupOpenPrefs.Derived = derivedDetails.open;
    writeGroupOpenPrefs(groupOpenPrefs);
  });
  const derivedSummary = document.createElement("summary");
  derivedSummary.textContent = "Derived";
  derivedDetails.appendChild(derivedSummary);

  const derived = document.createElement("dl");
  derived.className = "params-derived";
  derivedDetails.appendChild(derived);
  form.appendChild(derivedDetails);

  function refreshDerived(): void {
    const diameterInput = inputs.get("mill.diameter_m");
    const modeInput = inputs.get("mill.speed_mode");
    const valueInput = inputs.get("mill.speed_value");
    const ballDiameterInput = inputs.get("media.ball_diameter_m");
    const mediaFillInput = inputs.get("media.fill_fraction");
    const packingFractionInput = inputs.get("media.packing_fraction_2d");
    const mediaDensityInput = inputs.get("media.density_kg_m3");
    const maxBallsInput = inputs.get("simulation.max_balls");
    const substepsInput = inputs.get("simulation.substeps");
    const slurryFillInput = inputs.get("slurry.fill_fraction");
    if (
      !diameterInput ||
      !modeInput ||
      !valueInput ||
      !ballDiameterInput ||
      !mediaFillInput ||
      !packingFractionInput ||
      !mediaDensityInput ||
      !maxBallsInput ||
      !substepsInput ||
      !slurryFillInput
    ) {
      derived.replaceChildren();
      return;
    }
    const diameterM = Number(diameterInput.value);
    const mill = { diameter_m: diameterM, speed_mode: modeInput.value, speed_value: Number(valueInput.value) };
    if (!(diameterM > 0)) {
      derived.replaceChildren();
      if (mediaWarn) mediaWarn.hidden = true;
      return;
    }

    const ballDiameterField = SCHEMA.find((f) => f.path === "media.ball_diameter_m");
    const ballDiameterM = ballDiameterField
      ? fromDisplayValue(ballDiameterField, Number(ballDiameterInput.value))
      : Number(ballDiameterInput.value);
    const media = {
      ball_diameter_m: ballDiameterM,
      fill_fraction: Number(mediaFillInput.value),
      packing_fraction_2d: Number(packingFractionInput.value),
      density_kg_m3: Number(mediaDensityInput.value),
    };
    const maxBalls = Number(maxBallsInput.value);
    const substeps = Number(substepsInput.value);
    const slurryFill = Number(slurryFillInput.value);

    const nTrue = trueBallCount(diameterM, media);
    const eff = effectiveMedia(diameterM, media, maxBalls);
    const fluidCount = fluidParticleCountEstimate(diameterM / 2, DEFAULT_FLUID_RESOLUTION, slurryFill);
    const substepDisp = substepDisplacementOverDiameter(mill, substeps, eff.diameterM);
    // Mirrors pbf.rs's `dx = drum_radius_m / resolution`, same dx used internally by
    // fluidParticleCountEstimate above. A ratio > 1 means a single fluid particle is wider than a
    // ball -- the "unresolved ball<->fluid coupling" regime (docs/PHYSICS.md ss6) where the
    // coupling's per-particle taper weighting is a poor approximation (many fluid particles
    // overlapping one ball, or one fluid particle spanning several balls); this was part of the
    // fluidised-charge bug's root cause.
    const fluidDxM = diameterM / 2 / DEFAULT_FLUID_RESOLUTION;
    const fluidSpacingVsBall = eff.diameterM > 0 ? fluidDxM / eff.diameterM : 0;

    const rows: [string, string, boolean?][] = [
      ["Critical speed (Nc)", `${criticalSpeedRpm(diameterM).toFixed(1)} rpm`],
      ["Current speed", `${rpmOf(mill).toFixed(1)} rpm (${percentCriticalOf(mill).toFixed(0)}% Nc)`],
      ["True ball count (N_true)", nTrue.toFixed(0)],
      ["Simulated balls (N_sim)", String(eff.ballCount)],
      ["Coarse-graining (k)", eff.scaleFactor.toFixed(3)],
      ["Effective ball diameter", `${(eff.diameterM * 1000).toFixed(2)} mm`],
      ["Fluid spacing vs ball diameter", fluidSpacingVsBall.toFixed(2), fluidSpacingVsBall > 1.0],
      ["Est. fluid particles", String(fluidCount)],
      ["Sub-step displacement / diameter", substepDisp.toFixed(3), substepDisp > 0.3],
    ];

    derived.replaceChildren();
    for (const [label, value, warn] of rows) {
      const dt = document.createElement("dt");
      dt.textContent = label;
      const dd = document.createElement("dd");
      dd.textContent = value;
      dd.classList.toggle("is-warn", Boolean(warn));
      derived.append(dt, dd);
    }

    if (mediaWarn) {
      if (eff.scaleFactor > 1) {
        mediaWarn.hidden = false;
        mediaWarn.textContent =
          `Coarse-graining is active (k = ${eff.scaleFactor.toFixed(3)}). Ball diameter has no ` +
          `effect on the simulation at this Max balls -- the solver runs ${eff.ballCount} balls of ` +
          `${(eff.diameterM * 1000).toFixed(2)} mm regardless. Raise Max balls above ${nTrue.toFixed(0)} ` +
          `to simulate the true diameter.`;
      } else {
        mediaWarn.hidden = true;
      }
    }
  }

  function clearInvalid(): void {
    for (const row of rowsByPath.values()) {
      row.classList.remove("is-invalid");
    }
    errorBanner.hidden = true;
    errorBanner.textContent = "";
  }

  function validate(): FieldFailure[] {
    const failures: FieldFailure[] = [];
    for (const field of SCHEMA) {
      if (field.type !== "number") continue;
      const input = inputs.get(field.path);
      if (!input) continue;
      const raw = Number((input as HTMLInputElement).value);
      if ((input as HTMLInputElement).value.trim() === "" || !Number.isFinite(raw)) {
        failures.push({ path: field.path, message: `${field.label}: required` });
        continue;
      }
      if (field.min !== undefined && raw < field.min) {
        failures.push({ path: field.path, message: `${field.label}: must be >= ${field.min}` });
        continue;
      }
      if (field.max !== undefined && raw > field.max) {
        failures.push({ path: field.path, message: `${field.label}: must be <= ${field.max}` });
      }
    }
    return failures;
  }

  form.addEventListener("submit", (event) => {
    event.preventDefault();
    const failures = validate();
    if (failures.length > 0) {
      for (const row of rowsByPath.values()) {
        row.classList.remove("is-invalid");
      }
      let firstRow: HTMLElement | null = null;
      for (const failure of failures) {
        const row = rowsByPath.get(failure.path);
        if (!row) continue;
        row.classList.add("is-invalid");
        if (!firstRow) firstRow = row;
        const details = row.closest("details");
        if (details) details.open = true;
      }
      errorBanner.hidden = false;
      errorBanner.textContent = failures.map((f) => f.message).join("; ");
      setStatus("error");
      firstRow?.scrollIntoView({ block: "center" });
      return;
    }

    clearInvalid();
    let next = currentParams;
    for (const [path, input] of inputs) {
      const field = SCHEMA.find((f) => f.path === path);
      let value: unknown;
      if (field?.type === "number") value = fromDisplayValue(field, Number(input.value));
      else if (field?.type === "boolean") value = (input as HTMLInputElement).checked;
      else value = input.value;
      next = withPath(next, path, value);
    }
    onApply(next);
  });

  function setParams(params: ParamsJson): void {
    currentParams = params;
    for (const [path, input] of inputs) {
      const field = SCHEMA.find((f) => f.path === path);
      const value = getPath(params, path);
      if (value === undefined) continue;
      if (field?.type === "boolean") {
        (input as HTMLInputElement).checked = Boolean(value);
      } else {
        input.value = String(field ? toDisplayValue(field, value) : value);
      }
    }
    clearInvalid();
    refreshDerived();
    setStatus("applied");
  }

  function showError(message: string): void {
    errorBanner.hidden = false;
    errorBanner.textContent = message;
    setStatus("error");
  }

  refreshDerived();

  return { setParams, showError };
}
