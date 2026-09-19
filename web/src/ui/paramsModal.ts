// Native <dialog>-based parameters modal (v1: Mill / Media / Lifters / Slurry / Simulation
// groups, schema-driven from src/params/schema.ts). See docs/PLAN.md ss4.3. A Display tab and
// full validation UX/presets land in M5; for now Apply always resets the simulation (see
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
// SimulationParams::default). Not exposed in SCHEMA (out of scope for this modal), so it's
// hardcoded here purely as a static estimate input for the "Est. fluid particles" readout.
const DEFAULT_FLUID_RESOLUTION = 40;

export interface ParamsModal {
  /** Opens the modal, pre-filled from `params`. */
  open(params: ParamsJson): void;
}

export function createParamsModal(onApply: (params: ParamsJson) => void): ParamsModal {
  const dialog = document.createElement("dialog");
  dialog.className = "params-modal";

  const form = document.createElement("form");
  form.method = "dialog";
  dialog.appendChild(form);

  const title = document.createElement("h2");
  title.textContent = "Parameters";
  form.appendChild(title);

  const inputs = new Map<string, HTMLInputElement | HTMLSelectElement>();
  let currentParams: ParamsJson = {};

  for (const group of GROUPS) {
    const fieldset = document.createElement("fieldset");
    const legend = document.createElement("legend");
    legend.textContent = group;
    fieldset.appendChild(legend);

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
      input.addEventListener("input", () => refreshDerived());
      inputs.set(field.path, input);
      row.appendChild(input);
      fieldset.appendChild(row);
    }
    form.appendChild(fieldset);
  }

  const derived = document.createElement("dl");
  derived.className = "params-derived";
  form.appendChild(derived);

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
    const timeScaleInput = inputs.get("simulation.time_scale");
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
      !timeScaleInput ||
      !slurryFillInput
    ) {
      derived.replaceChildren();
      return;
    }
    const diameterM = Number(diameterInput.value);
    const mill = { diameter_m: diameterM, speed_mode: modeInput.value, speed_value: Number(valueInput.value) };
    if (!(diameterM > 0)) {
      derived.replaceChildren();
      return;
    }

    const ballDiameterField = SCHEMA.find((f) => f.path === "media.ball_diameter_m");
    const ballDiameterM = ballDiameterField ? fromDisplayValue(ballDiameterField, Number(ballDiameterInput.value)) : Number(ballDiameterInput.value);
    const media = {
      ball_diameter_m: ballDiameterM,
      fill_fraction: Number(mediaFillInput.value),
      packing_fraction_2d: Number(packingFractionInput.value),
      density_kg_m3: Number(mediaDensityInput.value),
    };
    const maxBalls = Number(maxBallsInput.value);
    const substeps = Number(substepsInput.value);
    const timeScale = Number(timeScaleInput.value);
    const slurryFill = Number(slurryFillInput.value);

    const nTrue = trueBallCount(diameterM, media);
    const eff = effectiveMedia(diameterM, media, maxBalls);
    const fluidCount = fluidParticleCountEstimate(diameterM / 2, DEFAULT_FLUID_RESOLUTION, slurryFill);
    const substepDisp = substepDisplacementOverDiameter(mill, timeScale, substeps, eff.diameterM);

    const rows: [string, string, boolean?][] = [
      ["Critical speed (Nc)", `${criticalSpeedRpm(diameterM).toFixed(1)} rpm`],
      ["Current speed", `${rpmOf(mill).toFixed(1)} rpm (${percentCriticalOf(mill).toFixed(0)}% Nc)`],
      ["True ball count (N_true)", nTrue.toFixed(0)],
      ["Simulated balls (N_sim)", String(eff.ballCount)],
      ["Coarse-graining (k)", eff.scaleFactor.toFixed(3)],
      ["Effective ball diameter", `${(eff.diameterM * 1000).toFixed(2)} mm`],
      ["Est. fluid particles", String(fluidCount)],
      ["Sub-step displacement / diameter", substepDisp.toFixed(3), substepDisp > 0.5],
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
  }

  const footer = document.createElement("div");
  footer.className = "params-footer";
  const cancelBtn = document.createElement("button");
  cancelBtn.type = "button";
  cancelBtn.textContent = "Cancel";
  cancelBtn.addEventListener("click", () => dialog.close());
  const applyBtn = document.createElement("button");
  applyBtn.type = "submit";
  applyBtn.textContent = "Apply (resets simulation)";
  footer.append(cancelBtn, applyBtn);
  form.appendChild(footer);

  form.addEventListener("submit", (event) => {
    event.preventDefault();
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
    dialog.close();
  });

  document.body.appendChild(dialog);

  return {
    open(params: ParamsJson) {
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
      refreshDerived();
      dialog.showModal();
    },
  };
}
