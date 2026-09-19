// Toolbar (v1: Play/Pause, Step, Reset, Parameters). Presets/Screenshot land in M5 (see
// docs/PLAN.md ss4.4).

export interface ToolbarCallbacks {
  onTogglePlay(): boolean; // returns the new "running" state, so the button label can reflect it
  onStep(): void;
  onReset(): void;
  onToggleParams(): boolean; // returns the new "params panel visible" state, reflected via aria-pressed
  initialParamsVisible: boolean;
  onTogglePanel(): boolean; // returns the new "panel visible" state, reflected via aria-pressed
  initialPanelVisible: boolean;
}

export function createToolbar(callbacks: ToolbarCallbacks): HTMLElement {
  const el = document.createElement("div");
  el.className = "toolbar";

  const playPauseBtn = document.createElement("button");
  playPauseBtn.type = "button";
  playPauseBtn.textContent = "Pause";
  playPauseBtn.addEventListener("click", () => {
    const running = callbacks.onTogglePlay();
    playPauseBtn.textContent = running ? "Pause" : "Play";
  });

  const stepBtn = document.createElement("button");
  stepBtn.type = "button";
  stepBtn.textContent = "Step";
  stepBtn.addEventListener("click", () => callbacks.onStep());

  const resetBtn = document.createElement("button");
  resetBtn.type = "button";
  resetBtn.textContent = "Reset";
  resetBtn.addEventListener("click", () => callbacks.onReset());

  const paramsBtn = document.createElement("button");
  paramsBtn.type = "button";
  paramsBtn.textContent = "Parameters";
  paramsBtn.setAttribute("aria-pressed", String(callbacks.initialParamsVisible));
  paramsBtn.addEventListener("click", () => {
    const visible = callbacks.onToggleParams();
    paramsBtn.setAttribute("aria-pressed", String(visible));
  });

  const panelBtn = document.createElement("button");
  panelBtn.type = "button";
  panelBtn.textContent = "Panel";
  panelBtn.setAttribute("aria-pressed", String(callbacks.initialPanelVisible));
  panelBtn.addEventListener("click", () => {
    const visible = callbacks.onTogglePanel();
    panelBtn.setAttribute("aria-pressed", String(visible));
  });

  el.append(playPauseBtn, stepBtn, resetBtn, paramsBtn, panelBtn);
  return el;
}
