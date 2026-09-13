// Toolbar (v1: Play/Pause, Step, Reset, Parameters). Presets/Screenshot land in M5 (see
// docs/PLAN.md ss4.4).

export interface ToolbarCallbacks {
  onTogglePlay(): boolean; // returns the new "running" state, so the button label can reflect it
  onStep(): void;
  onReset(): void;
  onOpenParams(): void;
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
  paramsBtn.addEventListener("click", () => callbacks.onOpenParams());

  el.append(playPauseBtn, stepBtn, resetBtn, paramsBtn);
  document.body.appendChild(el);
  return el;
}
