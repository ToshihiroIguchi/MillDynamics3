// Parameters left panel (v1: Mill / Media / Lifters / Slurry / Numerical accuracy groups, schema-driven
// from src/params/schema.ts, plus a Derived group). See docs/PLAN.md ss4.3. Replaces the earlier
// <dialog>-based modal (ui/paramsModal.ts) with a persistent <aside> panel that sits beside the
// canvas, symmetric with the metrics panel on the right. A Display tab and full validation UX
// beyond field-level range checks land in M5; for now Apply always resets the simulation (see
// schema.ts's v1 simplification note).

import type { ParamsJson } from "../protocol";
import {
  type CoarseGrainingMode,
  criticalSpeedRpm,
  effectiveFluidResolution,
  effectiveMedia,
  effectiveSubsteps,
  fluidParticleCountEstimate,
  interstitialFilling,
  percentCriticalOf,
  rpmOf,
  substepDisplacementOverDiameter,
  TRUE_COUNT_BLOCK_THRESHOLD,
  TRUE_COUNT_WARN_THRESHOLD,
  trueBallCount,
  trueCountGuard,
} from "../params/derived";
import { createHelpIcon } from "./helpIcon";
import { MEDIA_MATERIALS } from "../params/materials";
import { QUALITY_PRESETS, PRESET_FIELD_PATHS, presetFieldValues, presetMatchesFields } from "../params/presets";
import {
  fromDisplayValue,
  GROUPS,
  getPath,
  paramsChangeRequiresReset,
  SCHEMA,
  toDisplayValue,
  withPath,
} from "../params/schema";

const GROUPS_STORAGE_KEY = "milldynamics.paramsGroups";

const DEFAULT_GROUP_OPEN: Record<string, boolean> = {
  Mill: true,
  Media: true,
  Lifters: false,
  Slurry: false,
  "Numerical accuracy": false,
  Derived: true,
  // Collapsed "Advanced" sub-sections, keyed `${group}/Advanced` (readGroupOpenPrefs only keeps
  // keys listed here).
  "Media/Advanced": false,
  "Lifters/Advanced": false,
  "Slurry/Advanced": false,
  "Numerical accuracy/Advanced": false,
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
    // "Numerical accuracy" was called "Simulation" in earlier versions: honour the old stored key
    // when the new one is absent, so a user's open/closed choice survives the rename.
    if (typeof parsed["Numerical accuracy"] !== "boolean" && typeof parsed.Simulation === "boolean") {
      result["Numerical accuracy"] = parsed.Simulation;
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
  // main.ts's onApply callback hot-applies via "setParams" (Simulation::set_params) unless a
  // changed field needs a full reset (params/schema.ts's `paramsChangeRequiresReset`). While the
  // form is being edited, `setStatus("edited", ...)` below labels this button with which of the
  // two the pending edit will do ("Apply (live)" / "Apply & reset"); otherwise it reads "Apply".
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

  /**
   * `resetLabel` is only used for the "edited" state: the label of the first changed field that
   * needs a full simulation reset, or `null` when the pending edit can be applied live.
   */
  function setStatus(state: "applied" | "edited" | "error" | "applying", resetLabel: string | null = null): void {
    status.classList.remove("is-applied", "is-edited", "is-error", "is-applying");
    // Header itself (not just the status span) picks up the highlight: the header is
    // `position: sticky` so it stays on-screen while the panel is scrolled, but a small
    // grey-on-dark status span is still easy to miss in peripheral vision while typing further
    // down -- a background tint on the whole sticky bar is not (user feedback: the Apply
    // button/status is easy to lose track of while editing fields lower in the panel).
    header.classList.remove("is-edited");
    applyBtn.textContent = "Apply";
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
      if (resetLabel !== null) {
        status.textContent = `Edited — reset needed: ${resetLabel}`;
        applyBtn.textContent = "Apply & reset";
      } else {
        status.textContent = "Edited — not applied";
        applyBtn.textContent = "Apply (live)";
      }
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
  // Advanced sub-sections: the <details> plus the field paths it holds (to hide it when empty).
  const advancedByGroup = new Map<string, { details: HTMLDetailsElement; paths: string[] }>();
  let currentParams: ParamsJson = {};
  let mediaWarn: HTMLParagraphElement | null = null;
  let presetSelect: HTMLSelectElement | null = null;
  let presetNote: HTMLParagraphElement | null = null;
  let materialSelect: HTMLSelectElement | null = null;
  let materialNote: HTMLParagraphElement | null = null;
  let coarseGrainingWarn: HTMLParagraphElement | null = null;
  let coarseModeNotice: HTMLParagraphElement | null = null;

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
    const matches = presetMatchesFields(preset, (path) => inputs.get(path)?.value);
    if (!matches) {
      presetSelect.value = "";
      if (presetNote) presetNote.hidden = true;
    }
  }

  /**
   * Keeps the Media material dropdown in sync with `media.density_kg_m3`: selects the material
   * whose density matches the current field value (and shows its note), or falls back to
   * "(custom)" when no material matches -- e.g. the user picked "Stainless steel" and then
   * hand-edited "Media density" away from 7700, or a freshly-loaded params object's density
   * happens to equal a known material's value. Unlike `syncPresetSelectWithFields` (which only
   * resets an already-selected preset), this searches all of `MEDIA_MATERIALS` every time since
   * there is a single backing field to match against, so it also handles selecting the right
   * material on initial load.
   */
  function syncMaterialSelectWithFields(): void {
    if (!materialSelect) return;
    const densityInput = inputs.get("media.density_kg_m3") as HTMLInputElement | undefined;
    const material =
      densityInput !== undefined ? MEDIA_MATERIALS.find((m) => m.densityKgM3 === Number(densityInput.value)) : undefined;
    if (material) {
      materialSelect.value = material.id;
      if (materialNote) {
        materialNote.textContent = material.note;
        materialNote.hidden = false;
      }
    } else {
      materialSelect.value = "";
      if (materialNote) materialNote.hidden = true;
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

    if (group === "Numerical accuracy") {
      const presetRow = document.createElement("label");
      presetRow.className = "params-row";
      const presetLabel = document.createElement("span");
      presetLabel.className = "params-label";
      presetLabel.append("Quality preset");
      const presetHelp = createHelpIcon(
        "One-click (Max balls, Slurry resolution) pair trading speed against accuracy. Apply is still needed, and both fields reset the simulation.",
      );
      if (presetHelp) presetLabel.appendChild(presetHelp);
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
        // All fields are set *before* any "input" event is dispatched, for the reason above.
        const targets = presetFieldValues(preset)
          .map(({ path, value }) => ({ input: inputs.get(path), value }))
          .filter((t): t is { input: HTMLInputElement | HTMLSelectElement; value: string } => t.input !== undefined);
        for (const { input, value } of targets) input.value = value;
        for (const { input } of targets) {
          input.dispatchEvent(new Event("input", { bubbles: true }));
        }
        if (presetNote) {
          presetNote.textContent = preset.note;
          presetNote.hidden = false;
        }
      });
    }

    if (group === "Media") {
      const materialRow = document.createElement("label");
      materialRow.className = "params-row";
      const materialLabel = document.createElement("span");
      materialLabel.className = "params-label";
      materialLabel.append("Media material");
      const materialHelp = createHelpIcon(
        "Fills in Media density from a common ball material. Choose (custom) to type your own density.",
      );
      if (materialHelp) materialLabel.appendChild(materialHelp);
      materialRow.appendChild(materialLabel);

      materialSelect = document.createElement("select");
      const materialPlaceholderOption = document.createElement("option");
      materialPlaceholderOption.value = "";
      materialPlaceholderOption.textContent = "(custom)";
      materialSelect.appendChild(materialPlaceholderOption);
      for (const material of MEDIA_MATERIALS) {
        const opt = document.createElement("option");
        opt.value = material.id;
        opt.textContent = material.label;
        materialSelect.appendChild(opt);
      }
      materialRow.appendChild(materialSelect);
      details.appendChild(materialRow);

      materialNote = document.createElement("p");
      materialNote.className = "params-note";
      materialNote.hidden = true;
      details.appendChild(materialNote);

      // Sets media.density_kg_m3's *displayed* (not-yet-applied) value, exactly as if the user had
      // typed it in directly -- Apply still commits (and still resets; the field is
      // `resetRequired`, see schema.ts).
      materialSelect.addEventListener("change", () => {
        const material = MEDIA_MATERIALS.find((m) => m.id === materialSelect!.value);
        if (!material) {
          if (materialNote) materialNote.hidden = true;
          return;
        }
        const densityInput = inputs.get("media.density_kg_m3");
        if (densityInput) densityInput.value = String(material.densityKgM3);
        if (densityInput) densityInput.dispatchEvent(new Event("input", { bubbles: true }));
        if (materialNote) {
          materialNote.textContent = material.note;
          materialNote.hidden = false;
        }
      });
    }

    const groupFields = SCHEMA.filter((f) => f.group === group);
    let advancedDetails: HTMLDetailsElement | null = null;
    if (groupFields.some((f) => f.advanced)) {
      const advancedKey = `${group}/Advanced`;
      const sub = document.createElement("details");
      sub.className = "params-subgroup";
      sub.open = groupOpenPrefs[advancedKey] ?? false;
      sub.addEventListener("toggle", () => {
        groupOpenPrefs[advancedKey] = sub.open;
        writeGroupOpenPrefs(groupOpenPrefs);
      });
      const subSummary = document.createElement("summary");
      subSummary.textContent = "Advanced";
      sub.appendChild(subSummary);
      advancedDetails = sub;
      advancedByGroup.set(group, { details: sub, paths: [] });
    }

    for (const field of groupFields) {
      const row = document.createElement("label");
      row.className = "params-row";

      // min/max (when both are defined) are already in display units per schema.ts's own comment
      // on `displayScale` -- no conversion needed here. `noValidate` on the form means the
      // browser never surfaces the input's min/max attributes (set below) to the user on its own,
      // so the valid range is otherwise invisible; show it in the label and as a hover tooltip.
      const rangeSuffix = field.min !== undefined && field.max !== undefined ? ` [${field.min}–${field.max}]` : "";
      const labelText = document.createElement("span");
      labelText.className = "params-label";
      labelText.append((field.unit ? `${field.label} (${field.unit})` : field.label) + rangeSuffix);
      if (field.resetRequired) {
        // Marks fields whose change only takes effect through a full simulation reset (schema.ts's
        // `resetRequired`); the tooltip note is appended here so it never needs repeating per field.
        const glyph = document.createElement("span");
        glyph.className = "reset-glyph";
        glyph.textContent = "↻";
        glyph.title = "Changing this restarts the simulation.";
        labelText.appendChild(glyph);
      }
      const helpIcon = createHelpIcon(
        field.help && field.resetRequired ? `${field.help} Changing this restarts the simulation.` : field.help,
      );
      if (helpIcon) labelText.appendChild(helpIcon);
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
        setStatus("edited", paramsChangeRequiresReset(currentParams, readFormParams()));
        refreshDerived();
        if (PRESET_FIELD_PATHS.has(field.path)) {
          syncPresetSelectWithFields();
        }
        if (field.path === "media.density_kg_m3") {
          syncMaterialSelectWithFields();
        }
      });
      inputs.set(field.path, input);
      rowsByPath.set(field.path, row);
      row.appendChild(input);
      if (field.advanced && advancedDetails) {
        advancedDetails.appendChild(row);
        advancedByGroup.get(group)?.paths.push(field.path);
      } else {
        details.appendChild(row);
      }
      if (field.path === "simulation.coarse_graining_mode") {
        // Persistent amber notice directly under the mode row; shown whenever mode != off.
        coarseModeNotice = document.createElement("p");
        coarseModeNotice.className = "params-warn";
        coarseModeNotice.hidden = true;
        details.appendChild(coarseModeNotice);
      }
    }

    // After the basic rows, before the group-wide warn/note paragraphs below.
    if (advancedDetails) details.appendChild(advancedDetails);

    if (group === "Media") {
      mediaWarn = document.createElement("p");
      mediaWarn.className = "params-warn";
      mediaWarn.hidden = true;
      details.appendChild(mediaWarn);
    }

    if (group === "Numerical accuracy") {
      coarseGrainingWarn = document.createElement("p");
      coarseGrainingWarn.className = "params-warn";
      coarseGrainingWarn.hidden = true;
      details.appendChild(coarseGrainingWarn);
    }

    if (group === "Slurry") {
      // Static, always-visible one-line reminder (not a threshold-triggered alert like `mediaWarn`
      // above, so `params-note` rather than `params-warn`). The longer explanation lives in the
      // "Slurry resolution" field's help text (schema.ts) and docs/PHYSICS.md section 9.
      const slurryInfiltrationNote = document.createElement("p");
      slurryInfiltrationNote.className = "params-note";
      slurryInfiltrationNote.textContent =
        "Slurry barely seeps into a settled bed in this 2D model -- expected, not a bug (docs/PHYSICS.md §9).";
      const slurryAbout = document.createElement("details");
      slurryAbout.className = "note-details";
      const slurryAboutSummary = document.createElement("summary");
      slurryAboutSummary.textContent = "About slurry infiltration";
      slurryAbout.append(slurryAboutSummary, slurryInfiltrationNote);
      details.appendChild(slurryAbout);
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

  function setRowEnabled(path: string, enabled: boolean): void {
    const row = rowsByPath.get(path);
    const input = inputs.get(path);
    if (input) (input as HTMLInputElement).disabled = !enabled;
    row?.classList.toggle("is-ineffective", !enabled);
  }

  /**
   * Shows/hides a field's row by dependency. A row marked invalid stays visible so a failed
   * validation never points at an invisible field (values are always sent to the core regardless).
   */
  function setRowVisible(path: string, visible: boolean): void {
    const row = rowsByPath.get(path);
    if (row) row.hidden = !visible && !row.classList.contains("is-invalid");
  }

  /** Hides an Advanced sub-section when every row inside it is hidden. */
  function refreshAdvancedVisibility(): void {
    for (const { details: sub, paths } of advancedByGroup.values()) {
      sub.hidden = paths.every((p) => rowsByPath.get(p)?.hidden === true);
    }
  }

  function trueCountBlockMessage(nTrue: number): string {
    return (
      `True ball count ${nTrue.toFixed(0)} exceeds ${TRUE_COUNT_BLOCK_THRESHOLD} with Coarse-graining Off: ` +
      `too slow to run. Increase Ball diameter or lower Fill fraction, or, if you really want a ` +
      `coarse-grained run, choose Coarse-graining Auto or Manual explicitly.`
    );
  }

  /** Error text when the form's true ball count is too large to run with coarse-graining Off, else null. */
  function trueCountError(): string | null {
    const dIn = inputs.get("mill.diameter_m");
    const bIn = inputs.get("media.ball_diameter_m");
    const fIn = inputs.get("media.fill_fraction");
    const pIn = inputs.get("media.packing_fraction_2d");
    const mIn = inputs.get("simulation.coarse_graining_mode");
    if (!dIn || !bIn || !fIn || !pIn || !mIn) return null;
    const dField = SCHEMA.find((f) => f.path === "mill.diameter_m");
    const bField = SCHEMA.find((f) => f.path === "media.ball_diameter_m");
    const drum = dField ? fromDisplayValue(dField, Number(dIn.value)) : Number(dIn.value);
    const ball = bField ? fromDisplayValue(bField, Number(bIn.value)) : Number(bIn.value);
    const n = trueBallCount(drum, {
      ball_diameter_m: ball,
      fill_fraction: Number(fIn.value),
      packing_fraction_2d: Number(pIn.value),
      density_kg_m3: 0,
    });
    return trueCountGuard(mIn.value, n) === "block" ? trueCountBlockMessage(n) : null;
  }

  /**
   * Applies the dependency rules that hide fields with no effect under the current settings.
   * Runs before refreshDerived's early returns so invalid rows are revealed even when the derived
   * readout cannot be computed. Missing inputs count as "visible".
   */
  function refreshVisibility(): void {
    const liftersCountInput = inputs.get("lifters.count");
    const slurryEnabledInput = inputs.get("slurry.enabled") as HTMLInputElement | undefined;
    const coarseMode = inputs.get("simulation.coarse_graining_mode")?.value || "off";
    // setRowVisible keeps invalid rows visible. Max balls stays visible in Manual (greyed out)
    // because the Manual warning in refreshDerived refers to it.
    const noLifters = liftersCountInput !== undefined && Number(liftersCountInput.value) === 0;
    for (const path of ["lifters.height_m", "lifters.base_width_m", "lifters.top_width_m", "lifters.phase_deg"]) {
      setRowVisible(path, !noLifters);
    }
    setRowVisible("simulation.max_balls", coarseMode !== "off");
    setRowVisible("simulation.coarse_graining_k", coarseMode === "manual");
    const slurryOff = slurryEnabledInput !== undefined && !slurryEnabledInput.checked;
    for (const field of SCHEMA) {
      if (field.group === "Slurry" && field.path !== "slurry.enabled") setRowVisible(field.path, !slurryOff);
    }
    refreshAdvancedVisibility();
  }

  function refreshDerived(): void {
    refreshVisibility();
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
    const coarseModeInput = inputs.get("simulation.coarse_graining_mode");
    const coarseKInput = inputs.get("simulation.coarse_graining_k");
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
      !slurryFillInput ||
      !coarseModeInput ||
      !coarseKInput
    ) {
      derived.replaceChildren();
      return;
    }
    const drumDiameterField = SCHEMA.find((f) => f.path === "mill.diameter_m");
    const diameterM = drumDiameterField
      ? fromDisplayValue(drumDiameterField, Number(diameterInput.value))
      : Number(diameterInput.value);
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
    const coarseGrainingMode = (coarseModeInput.value || "off") as CoarseGrainingMode;
    const coarseGrainingK = Number(coarseKInput.value);

    // Dependent rows: Max balls only matters in Auto, Manual k only in Manual. They are greyed out
    // and disabled (still readable; readFormParams reads their values regardless) otherwise.
    setRowEnabled("simulation.max_balls", coarseGrainingMode === "auto");
    setRowEnabled("simulation.coarse_graining_k", coarseGrainingMode === "manual");
    if (coarseModeNotice) {
      coarseModeNotice.hidden = coarseGrainingMode === "off";
      coarseModeNotice.textContent =
        "Coarse-grained run: results are not physically reliable. Kinetic energy, per-impact stress " +
        "statistics and the shoulder angle are not preserved and power draw is off by 17-39 % at " +
        "large D/d. Use only for qualitative charge-shape checks; set Coarse-graining to Off for any measurement.";
    }

    const nTrue = trueBallCount(diameterM, media);
    const eff = effectiveMedia(diameterM, media, maxBalls, coarseGrainingMode, coarseGrainingK);
    const effectiveResolution = effectiveFluidResolution(
      diameterM,
      media,
      maxBalls,
      resolution,
      coarseGrainingMode,
      coarseGrainingK,
    );
    const fluidCount = fluidParticleCountEstimate(diameterM / 2, effectiveResolution, slurryFill);
    const substepDisp = substepDisplacementOverDiameter(mill, substeps, eff.diameterM);
    // Mirrors pbf.rs's `dx = drum_radius_m / resolution`, same dx used internally by
    // fluidParticleCountEstimate above. A ratio > 1 means a single fluid particle is wider than a
    // ball -- the "unresolved ball<->fluid coupling" regime (docs/PHYSICS.md ss6) where the
    // coupling's per-particle taper weighting is a poor approximation (many fluid particles
    // overlapping one ball, or one fluid particle spanning several balls); this was part of the
    // fluidised-charge bug's root cause, and (h = 2*dx, so this ratio is h/(2*ball_diameter), a
    // fixed multiple of the h/r wetting-resolution ratio) directly explains why a thin wetting
    // film looks like a few oversized blobs rather than a smooth coating once this exceeds ~0.5.
    // Uses `effectiveResolution` (Params::effective_fluid_resolution), the resolution actually used
    // to seed the fluid lattice, not the raw `resolution` input -- see the new "Effective slurry
    // resolution" row below for when the two differ.
    const fluidDxM = diameterM / 2 / effectiveResolution;
    const fluidSpacingVsBall = eff.diameterM > 0 ? fluidDxM / eff.diameterM : 0;
    // Uses `effectiveSubsteps` (Params::effective_substeps), the sub-step count actually used by
    // the solver, not the raw `substeps` input -- see the "Effective sub-steps" row below for when
    // the two differ (coarse-graining shrinking `d_eff` enough to threaten the XPBD `< 1`
    // displacement-ratio criterion, docs/PHYSICS.md §9).
    const effectiveSubstepCount = effectiveSubsteps(mill, media, maxBalls, substeps, coarseGrainingMode, coarseGrainingK);
    const u = interstitialFilling(slurryFill, media);
    const coarseGrainingActive = eff.scaleFactor > 1;

    const rows: [string, string, boolean | undefined, string][] = [
      [
        "Critical speed (Nc)",
        `${criticalSpeedRpm(diameterM).toFixed(1)} rpm`,
        undefined,
        "Speed at which a ball at the wall just stops falling and centrifuges, from the drum diameter.",
      ],
      [
        "Current speed",
        `${rpmOf(mill).toFixed(1)} rpm (${percentCriticalOf(mill).toFixed(0)}% Nc)`,
        undefined,
        "The drum speed entered above, shown both in rpm and as a percentage of critical speed.",
      ],
      [
        "True ball count (N_true)",
        nTrue.toFixed(0),
        undefined,
        "Number of balls a real mill with these Media settings would hold.",
      ],
      [
        "Simulated balls (N_sim)",
        String(eff.ballCount),
        undefined,
        "Number of balls actually simulated; lower than N_true when coarse-graining is active.",
      ],
      [
        "Coarse-graining",
        coarseGrainingActive ? `ON (k = ${eff.scaleFactor.toFixed(3)})` : "OFF (k = 1, true size)",
        undefined,
        "Whether balls are enlarged by a factor k to keep the count manageable while preserving charge mass and footprint.",
      ],
      [
        "Effective ball diameter",
        `${(eff.diameterM * 1000).toFixed(2)} mm`,
        undefined,
        "Diameter of the balls actually simulated (true diameter times k).",
      ],
      [
        "Interstitial filling (U)",
        u.toFixed(3),
        undefined,
        "Slurry fill divided by the void space of the ball bed (J x void fraction). About 1 means the slurry just fills the gaps; above 1 leaves a pool over the bed.",
      ],
      [
        "Fluid spacing vs ball diameter",
        fluidSpacingVsBall.toFixed(2),
        fluidSpacingVsBall > 1.0,
        "Slurry particle spacing divided by the effective ball diameter. Above 1 a fluid particle is wider than a ball, so ball-fluid coupling is poorly resolved; raise Slurry resolution.",
      ],
      [
        "Effective slurry resolution",
        `${effectiveResolution}`,
        effectiveResolution > resolution,
        "Slurry resolution actually used. The solver raises it above your setting when the fluid lattice would otherwise be wider than a ball.",
      ],
      [
        "Est. fluid particles",
        String(fluidCount),
        undefined,
        "Estimated number of slurry particles, the main driver of computation cost.",
      ],
      [
        "Effective sub-steps",
        `${effectiveSubstepCount}`,
        effectiveSubstepCount > substeps,
        "Sub-steps per frame actually used. The solver raises it above your setting (up to 16) when the balls are small enough that they would otherwise move too far per sub-step; this includes the default 2 mm balls (8 to 16).",
      ],
      [
        "Sub-step displacement / diameter",
        substepDisp.toFixed(3),
        substepDisp > 0.3,
        "Distance the drum wall moves per sub-step relative to the ball diameter. Should stay well below 1 for stable contacts.",
      ],
    ];

    derived.replaceChildren();
    for (const [label, value, warn, help] of rows) {
      const dt = document.createElement("dt");
      dt.append(label);
      const helpIcon = createHelpIcon(help);
      if (helpIcon) dt.appendChild(helpIcon);
      const dd = document.createElement("dd");
      dd.textContent = value;
      dd.classList.toggle("is-warn", Boolean(warn));
      derived.append(dt, dd);
    }

    const ballDiameterRow = rowsByPath.get("media.ball_diameter_m");
    if (mediaWarn) {
      if (coarseGrainingActive && coarseGrainingMode === "auto") {
        mediaWarn.hidden = false;
        mediaWarn.textContent =
          `Coarse-graining is active (k = ${eff.scaleFactor.toFixed(3)}). Ball diameter has no ` +
          `effect on the simulation at this Max balls -- the solver runs ${eff.ballCount} balls of ` +
          `${(eff.diameterM * 1000).toFixed(2)} mm regardless. Raise Max balls above ${nTrue.toFixed(0)} ` +
          `to simulate the true diameter.`;
      } else if (coarseGrainingActive) {
        mediaWarn.hidden = false;
        mediaWarn.textContent =
          `Coarse-graining is active (k = ${eff.scaleFactor.toFixed(3)}, mode: ${coarseGrainingMode}). ` +
          `Ball diameter has no effect on the simulation -- the solver runs ${eff.ballCount} balls of ` +
          `${(eff.diameterM * 1000).toFixed(2)} mm regardless. Change "Coarse-graining" in the ` +
          `Numerical accuracy group to simulate the true diameter.`;
      } else {
        mediaWarn.hidden = true;
      }
      // Dim (not disable) the ball-diameter input itself so the warning is hard to miss even
      // without reading the banner text -- the value stays editable for when the user later
      // changes Max balls / the coarse-graining mode.
      ballDiameterRow?.classList.toggle("is-ineffective", coarseGrainingActive);
    }

    // Coarse-graining Off/Manual can push N_sim or the effective ball diameter into a regime Auto
    // would never produce (Auto always targets/stays under Max balls and never grows the effective
    // diameter past the true one relative to the drum) -- surface that explicitly rather than
    // letting a user silently run an unreasonably heavy or unphysical configuration.
    if (coarseGrainingWarn) {
      const messages: string[] = [];
      const nSimVsTarget = maxBalls > 0 ? eff.ballCount / maxBalls : 0;
      const guard = trueCountGuard(coarseGrainingMode, nTrue);
      if (guard === "warn") {
        messages.push(
          `Coarse-graining is Off: the solver will simulate the true ball count (${eff.ballCount}). ` +
            `Above ${TRUE_COUNT_WARN_THRESHOLD} balls the simulation runs much slower than real time.`,
        );
      } else if (guard === "block") {
        messages.push(trueCountBlockMessage(nTrue));
      }
      if (coarseGrainingMode === "manual") {
        if (nSimVsTarget > 2) {
          messages.push(
            `Manual k = ${coarseGrainingK.toFixed(3)} leaves N_sim = ${eff.ballCount}, ` +
              `${nSimVsTarget.toFixed(1)}x above your Max balls target (${maxBalls}). Raise k to shrink N_sim.`,
          );
        }
        if (eff.diameterM >= diameterM) {
          messages.push(
            `Manual k = ${coarseGrainingK.toFixed(3)} makes the effective ball diameter ` +
              `(${(eff.diameterM * 1000).toFixed(1)} mm) as large as or larger than the drum itself ` +
              `(${(diameterM * 1000).toFixed(1)} mm) -- unphysical. Lower k.`,
          );
        } else if (eff.ballCount < 3) {
          messages.push(
            `Manual k = ${coarseGrainingK.toFixed(3)} leaves only ${eff.ballCount} simulated ball(s), too ` +
              `coarse to represent a charge. Lower k.`,
          );
        }
      }
      coarseGrainingWarn.hidden = messages.length === 0;
      coarseGrainingWarn.textContent = messages.join(" ");
    }
  }

  /** Builds the params object the form currently describes (SI units) on top of `currentParams`. */
  function readFormParams(): ParamsJson {
    let next = currentParams;
    for (const [path, input] of inputs) {
      const field = SCHEMA.find((f) => f.path === path);
      let value: unknown;
      if (field?.type === "number") value = fromDisplayValue(field, Number(input.value));
      else if (field?.type === "boolean") value = (input as HTMLInputElement).checked;
      else value = input.value;
      next = withPath(next, path, value);
    }
    return next;
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
    const countError = failures.length === 0 ? trueCountError() : null;
    if (countError) {
      errorBanner.hidden = false;
      errorBanner.textContent = countError;
      setStatus("error");
      rowsByPath.get("simulation.coarse_graining_mode")?.scrollIntoView({ block: "center" });
      return;
    }
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
        // Open every ancestor <details> (group and its Advanced sub-section).
        for (let d = row.closest("details"); d; d = d.parentElement?.closest("details") ?? null) {
          d.open = true;
        }
      }
      // Reveal invalid rows that a dependency rule would otherwise hide.
      refreshDerived();
      errorBanner.hidden = false;
      errorBanner.textContent = failures.map((f) => f.message).join("; ");
      setStatus("error");
      firstRow?.scrollIntoView({ block: "center" });
      return;
    }

    clearInvalid();
    refreshDerived();
    setStatus("applying");
    onApply(readFormParams());
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
    syncMaterialSelectWithFields();
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
