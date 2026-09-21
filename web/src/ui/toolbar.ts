// Toolbar (Play/Pause, Step, Reset, Parameters, Panel). Icon buttons with a native `title`
// tooltip on hover and an `aria-label` carrying the same accessible name the old text label had
// (tests locate these buttons by that name, e.g. `getByRole("button", { name: "Panel" })`) --
// see docs/PLAN.md ss4.4.

export interface ToolbarCallbacks {
  onTogglePlay(): boolean; // returns the new "running" state, so the icon can reflect it
  onStep(): void;
  onReset(): void;
  onToggleParams(): boolean; // returns the new "params panel visible" state, reflected via aria-pressed
  initialParamsVisible: boolean;
  onTogglePanel(): boolean; // returns the new "panel visible" state, reflected via aria-pressed
  initialPanelVisible: boolean;
}

/** 18x18 stroke-based icons, one per toolbar action. Kept inline (no icon font/library
 * dependency) since the toolbar only ever needs these five. */
const ICONS = {
  play: '<svg viewBox="0 0 18 18" width="18" height="18" fill="currentColor" aria-hidden="true"><path d="M5 3.5v11l9-5.5z"/></svg>',
  pause:
    '<svg viewBox="0 0 18 18" width="18" height="18" fill="currentColor" aria-hidden="true"><rect x="4" y="3.5" width="3.5" height="11"/><rect x="10.5" y="3.5" width="3.5" height="11"/></svg>',
  step: '<svg viewBox="0 0 18 18" width="18" height="18" fill="currentColor" aria-hidden="true"><path d="M4 3.5v11l7-5.5z"/><rect x="12.5" y="3.5" width="2" height="11"/></svg>',
  reset:
    '<svg viewBox="0 0 18 18" width="18" height="18" fill="none" stroke="currentColor" stroke-width="1.6" stroke-linecap="round" aria-hidden="true"><path d="M14 9a5 5 0 1 1-1.6-3.68"/><path d="M14 3.5v3.2h-3.2" stroke-linejoin="round"/></svg>',
  parameters:
    '<svg viewBox="0 0 18 18" width="18" height="18" fill="none" stroke="currentColor" stroke-width="1.6" stroke-linecap="round" aria-hidden="true"><line x1="3" y1="5" x2="15" y2="5"/><line x1="3" y1="9" x2="15" y2="9"/><line x1="3" y1="13" x2="15" y2="13"/><circle cx="7" cy="5" r="1.6" fill="currentColor" stroke="none"/><circle cx="12" cy="9" r="1.6" fill="currentColor" stroke="none"/><circle cx="6" cy="13" r="1.6" fill="currentColor" stroke="none"/></svg>',
  panel:
    '<svg viewBox="0 0 18 18" width="18" height="18" fill="none" stroke="currentColor" stroke-width="1.6" aria-hidden="true"><rect x="2.5" y="3" width="13" height="12" rx="1.5"/><line x1="11.5" y1="3" x2="11.5" y2="15"/></svg>',
};

/** Creates a toolbar button showing only an icon, with `label` as both the tooltip (native
 * `title`, shown on hover) and the accessible name (`aria-label`, what assistive tech and
 * `getByRole("button", { name })` see). */
function createIconButton(icon: string, label: string): HTMLButtonElement {
  const btn = document.createElement("button");
  btn.type = "button";
  btn.className = "toolbar-icon-btn";
  btn.innerHTML = icon;
  btn.title = label;
  btn.setAttribute("aria-label", label);
  return btn;
}

export function createToolbar(callbacks: ToolbarCallbacks): HTMLElement {
  const el = document.createElement("div");
  el.className = "toolbar";

  const playPauseBtn = createIconButton(ICONS.pause, "Pause");
  playPauseBtn.title = "Pause the simulation (Space)";
  playPauseBtn.addEventListener("click", () => {
    const running = callbacks.onTogglePlay();
    playPauseBtn.innerHTML = running ? ICONS.pause : ICONS.play;
    const label = running ? "Pause" : "Play";
    playPauseBtn.setAttribute("aria-label", label);
    playPauseBtn.title = running ? "Pause the simulation (Space)" : "Resume the simulation (Space)";
  });

  const stepBtn = createIconButton(ICONS.step, "Step");
  stepBtn.title = "Advance one simulation step";
  stepBtn.addEventListener("click", () => callbacks.onStep());

  const resetBtn = createIconButton(ICONS.reset, "Reset");
  resetBtn.title = "Reset the simulation to its initial state";
  resetBtn.addEventListener("click", () => callbacks.onReset());

  const paramsBtn = createIconButton(ICONS.parameters, "Parameters");
  paramsBtn.title = "Show/hide the parameters panel";
  paramsBtn.setAttribute("aria-pressed", String(callbacks.initialParamsVisible));
  paramsBtn.addEventListener("click", () => {
    const visible = callbacks.onToggleParams();
    paramsBtn.setAttribute("aria-pressed", String(visible));
  });

  const panelBtn = createIconButton(ICONS.panel, "Panel");
  panelBtn.title = "Show/hide the metrics panel";
  panelBtn.setAttribute("aria-pressed", String(callbacks.initialPanelVisible));
  panelBtn.addEventListener("click", () => {
    const visible = callbacks.onTogglePanel();
    panelBtn.setAttribute("aria-pressed", String(visible));
  });

  el.append(playPauseBtn, stepBtn, resetBtn, paramsBtn, panelBtn);
  return el;
}
