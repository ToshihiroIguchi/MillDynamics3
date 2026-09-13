// Native <dialog>-based parameters modal (v1: Mill / Media / Simulation groups, schema-driven from
// src/params/schema.ts). See docs/PLAN.md ss4.3. Slurry/Lifters/Display tabs and full validation
// UX/presets land in M2/M3/M5; for now Apply always resets the simulation (see schema.ts's v1
// simplification note).

import type { ParamsJson } from "../protocol";
import { criticalSpeedRpm, percentCriticalOf, rpmOf } from "../params/derived";
import { GROUPS, getPath, SCHEMA, withPath } from "../params/schema";

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

  const derived = document.createElement("p");
  derived.className = "params-derived";
  form.appendChild(derived);

  function refreshDerived(): void {
    const diameterInput = inputs.get("mill.diameter_m");
    const modeInput = inputs.get("mill.speed_mode");
    const valueInput = inputs.get("mill.speed_value");
    if (!diameterInput || !modeInput || !valueInput) return;
    const diameterM = Number(diameterInput.value);
    const mill = { diameter_m: diameterM, speed_mode: modeInput.value, speed_value: Number(valueInput.value) };
    if (!(diameterM > 0)) {
      derived.textContent = "";
      return;
    }
    derived.textContent =
      `Critical speed: ${criticalSpeedRpm(diameterM).toFixed(1)} rpm — ` +
      `current: ${rpmOf(mill).toFixed(1)} rpm (${percentCriticalOf(mill).toFixed(0)}% Nc)`;
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
      const value = field?.type === "number" ? Number(input.value) : input.value;
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
        const value = getPath(params, path);
        if (value !== undefined) {
          input.value = String(value);
        }
      }
      refreshDerived();
      dialog.showModal();
    },
  };
}
