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
  interstitialFilling,
  percentCriticalOf,
  rpmOf,
  substepDisplacementOverDiameter,
  trueBallCount,
} from "../params/derived";
import { QUALITY_PRESETS } from "../params/presets";
import { fromDisplayValue, GROUPS, getPath, SCHEMA, toDisplayValue, withPath } from "../params/schema";

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
  // Not "resets simulation" any more: main.ts's onApply callback hot-applies via "setParams"
  // (Simulation::set_params) unless the changed field actually needs one
  // (params/schema.ts's `paramsChangeRequiresReset`), in which case it resets -- this label can't
  // know which in advance without duplicating that decision here, so it stays deliberately
  // non-committal; see the browser console for which happened on a given Apply.
  applyBtn.textContent = "Apply";

  const revertBtn = document.createElement("button");
  revertBtn.type = "button";
  revertBtn.textContent = "Revert";
  revertBtn.addEventListener("click", () => setParams(currentParams));

  const status = document.createElement("span");
  status.className = "params-status is-applied";
  status.textContent = "Applied";

  header.append(applyBtn, revertBtn, status);
  form.appendChild(header);

  // Tracks whether setStatus() has ever run before -- guards the "applied" flash (below) so it
  // doesn't fire on the very first setParams() call during panel initialization, only on later
  // transitions that represent an actual user-triggered apply/revert settling.
  let statusInitialized = false;
  // Timeout id for the header's brief "applied" flash, so a rapid second "applied" transition
  // (e.g. two quick hot-applies) restarts the fade instead of leaving two overlapping timeouts.
  let appliedFlashTimeoutId: ReturnType<typeof setTimeout> | undefined;

  function setStatus(state: "applied" | "edited" | "error" | "applying"): void {
    status.classList.remove("is-applied", "is-edited", "is-error", "is-applying");
    // Header itself (not just the status span) picks up the highlight: the header is
    // `position: sticky` so it stays on-screen while the panel is scrolled, but a small
    // grey-on-dark status span is still easy to miss in peripheral vision while typing further
    // down -- a background tint on the whole sticky bar is not (user feedback: the Apply
    // button/status is easy to lose track of while editing fields lower in the panel).
    header.classList.remove("is-edited");
    if (state === "applied") {
      status.textContent = "Applied";
      status.classList.add("is-applied");
      // Brief flash on the header to confirm the apply landed -- same peripheral-vision reasoning
      // as the "is-edited" tint above, but skipped on the very first call (initial panel
      // population before the user has done anything).
      if (statusInitialized) {
        clearTimeout(appliedFlashTimeoutId);
        header.classList.add("is-applied-flash");
        appliedFlashTimeoutId = setTimeout(() => header.classList.remove("is-applied-flash"), 600);
      }
    } else if (state === "edited") {
      status.textContent = "Edited — not applied";
      status.classList.add("is-edited");
      header.classList.add("is-edited");
    } else if (state === "applying") {
      status.textContent = "Applying…";
      status.classList.add("is-applying");
    } else {
      status.textContent = "Error";
      status.classList.add("is-error");
    }
    // Disabled while an apply is in flight (hot-apply resolves synchronously so this never
    // renders a visible frame on that path; the async "init" reset path is the actual target --
    // it can take noticeably long for large ball counts / high slurry resolution).
    applyBtn.disabled = state === "applying";
    statusInitialized = true;
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
  let presetSelect: HTMLSelectElement | null = null;
  let presetNote: HTMLParagraphElement | null = null;

  /**
   * Resets the Quality preset dropdown to "(custom)" once `simulation.max_balls`/`resolution` no
   * longer match the selected preset's values -- e.g. the user picked "Accuracy" and then hand-
   * edited "Max balls" from 600 to 500. Without this the dropdown kept showing "Accuracy" after a
   * manual edit made it wrong (reported by an external review: the displayed selection and the
   * actual param values could silently disagree).
   */
  function syncPresetSelectWithFields(): void {
    if (!presetSelect || presetSelect.value === "") return;
    const preset = QUALITY_PRESETS.find((p) => p.id === presetSelect!.value);
    if (!preset) return;
    const maxBallsInput = inputs.get("simulation.max_balls") as HTMLInputElement | undefined;
    const resolutionInput = inputs.get("simulation.resolution") as HTMLInputElement | undefined;
    const matches =
      maxBallsInput !== undefined &&
      resolutionInput !== undefined &&
      Number(maxBallsInput.value) === preset.maxBalls &&
      Number(resolutionInput.value) === preset.resolution;
    if (!matches) {
      presetSelect.value = "";
      if (presetNote) presetNote.hidden = true;
    }
  }

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

    if (group === "Simulation") {
      const presetRow = document.createElement("label");
      presetRow.className = "params-row";
      const presetLabel = document.createElement("span");
      presetLabel.textContent = "Quality preset";
      presetRow.appendChild(presetLabel);

      presetSelect = document.createElement("select");
      const placeholderOption = document.createElement("option");
      placeholderOption.value = "";
      placeholderOption.textContent = "(custom)";
      presetSelect.appendChild(placeholderOption);
      for (const preset of QUALITY_PRESETS) {
        const opt = document.createElement("option");
        opt.value = preset.id;
        opt.textContent = preset.label;
        presetSelect.appendChild(opt);
      }
      presetRow.appendChild(presetSelect);
      details.appendChild(presetRow);

      presetNote = document.createElement("p");
      presetNote.className = "params-note";
      presetNote.hidden = true;
      details.appendChild(presetNote);

      // Sets simulation.max_balls/resolution's *displayed* (not-yet-applied) values, exactly as
      // if the user had typed them in directly -- Apply still commits (and still resets, both
      // fields are `resetRequired`, see schema.ts). Deliberately not auto-applied: the user should
      // see and confirm the derived-values readout (fluid particle count, fluid spacing vs. ball
      // diameter) below before committing to a reset.
      presetSelect.addEventListener("change", () => {
        const preset = QUALITY_PRESETS.find((p) => p.id === presetSelect!.value);
        if (!preset) {
          if (presetNote) presetNote.hidden = true;
          return;
        }
        const maxBallsInput = inputs.get("simulation.max_balls");
        const resolutionInput = inputs.get("simulation.resolution");
        // Both values are set *before* either "input" event is dispatched: the generic field
        // listener's `syncPresetSelectWithFields` call (added for the dropdown-desync fix above)
        // compares both fields against the preset every time either one changes, so dispatching
        // them one at a time would see max_balls already updated but resolution still stale,
        // read that as a mismatch, and immediately reset the dropdown back to "(custom)" right
        // after this handler had just set it.
        if (maxBallsInput) maxBallsInput.value = String(preset.maxBalls);
        if (resolutionInput) resolutionInput.value = String(preset.resolution);
        if (maxBallsInput) maxBallsInput.dispatchEvent(new Event("input", { bubbles: true }));
        if (resolutionInput) resolutionInput.dispatchEvent(new Event("input", { bubbles: true }));
        if (presetNote) {
          presetNote.textContent = preset.note;
          presetNote.hidden = false;
        }
      });
    }

    for (const field of SCHEMA.filter((f) => f.group === group)) {
      const row = document.createElement("label");
      row.className = "params-row";

      // min/max (when both are defined) are already in display units per schema.ts's own comment
      // on `displayScale` -- no conversion needed here. `noValidate` on the form means the
      // browser never surfaces the input's min/max attributes (set below) to the user on its own,
      // so the valid range is otherwise invisible; show it in the label and as a hover tooltip.
      const rangeSuffix = field.min !== undefined && field.max !== undefined ? ` [${field.min}–${field.max}]` : "";
      const labelText = document.createElement("span");
      labelText.textContent = (field.unit ? `${field.label} (${field.unit})` : field.label) + rangeSuffix;
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
        if (field.min !== undefined && field.max !== undefined) {
          numberInput.title = `Valid range: ${field.min}–${field.max}`;
        }
        input = numberInput;
      }
      input.addEventListener("input", () => {
        setStatus("edited");
        refreshDerived();
        if (field.path === "simulation.max_balls" || field.path === "simulation.resolution") {
          syncPresetSelectWithFields();
        }
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
    const resolutionInput = inputs.get("simulation.resolution");
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
      !resolutionInput ||
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
    const resolution = Number(resolutionInput.value);
    const slurryFill = Number(slurryFillInput.value);

    const nTrue = trueBallCount(diameterM, media);
    const eff = effectiveMedia(diameterM, media, maxBalls);
    const fluidCount = fluidParticleCountEstimate(diameterM / 2, resolution, slurryFill);
    const substepDisp = substepDisplacementOverDiameter(mill, substeps, eff.diameterM);
    // Mirrors pbf.rs's `dx = drum_radius_m / resolution`, same dx used internally by
    // fluidParticleCountEstimate above. A ratio > 1 means a single fluid particle is wider than a
    // ball -- the "unresolved ball<->fluid coupling" regime (docs/PHYSICS.md ss6) where the
    // coupling's per-particle taper weighting is a poor approximation (many fluid particles
    // overlapping one ball, or one fluid particle spanning several balls); this was part of the
    // fluidised-charge bug's root cause, and (h = 2*dx, so this ratio is h/(2*ball_diameter), a
    // fixed multiple of the h/r wetting-resolution ratio) directly explains why a thin wetting
    // film looks like a few oversized blobs rather than a smooth coating once this exceeds ~0.5.
    const fluidDxM = diameterM / 2 / resolution;
    const fluidSpacingVsBall = eff.diameterM > 0 ? fluidDxM / eff.diameterM : 0;
    const u = interstitialFilling(slurryFill, media);

    const rows: [string, string, boolean?][] = [
      ["Critical speed (Nc)", `${criticalSpeedRpm(diameterM).toFixed(1)} rpm`],
      ["Current speed", `${rpmOf(mill).toFixed(1)} rpm (${percentCriticalOf(mill).toFixed(0)}% Nc)`],
      ["True ball count (N_true)", nTrue.toFixed(0)],
      ["Simulated balls (N_sim)", String(eff.ballCount)],
      ["Coarse-graining (k)", eff.scaleFactor.toFixed(3)],
      ["Effective ball diameter", `${(eff.diameterM * 1000).toFixed(2)} mm`],
      ["Interstitial filling (U)", u.toFixed(3)],
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

    const ballDiameterRow = rowsByPath.get("media.ball_diameter_m");
    if (mediaWarn) {
      const coarseGrained = eff.scaleFactor > 1;
      if (coarseGrained) {
        mediaWarn.hidden = false;
        mediaWarn.textContent =
          `Coarse-graining is active (k = ${eff.scaleFactor.toFixed(3)}). Ball diameter has no ` +
          `effect on the simulation at this Max balls -- the solver runs ${eff.ballCount} balls of ` +
          `${(eff.diameterM * 1000).toFixed(2)} mm regardless. Raise Max balls above ${nTrue.toFixed(0)} ` +
          `to simulate the true diameter.`;
      } else {
        mediaWarn.hidden = true;
      }
      // Dim (not disable) the ball-diameter input itself so the warning is hard to miss even
      // without reading the banner text -- the value stays editable for when the user later
      // raises Max balls.
      ballDiameterRow?.classList.toggle("is-ineffective", coarseGrained);
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
    setStatus("applying");
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
