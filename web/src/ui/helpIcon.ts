// Small "?" help affordance shared by the params panel and the metrics panel: a focusable span
// whose native `title` tooltip carries the help text (keyboard users get it via aria-label).

/** Builds a `<span class="help-icon">?</span>` showing `help` as its tooltip, or `null` if there
 * is no help text to show. */
export function createHelpIcon(help: string | undefined): HTMLSpanElement | null {
  if (!help) return null;
  const icon = document.createElement("span");
  icon.className = "help-icon";
  icon.textContent = "?";
  icon.title = help;
  icon.setAttribute("aria-label", help);
  icon.tabIndex = 0;
  // Inside a <label>, a click would otherwise forward to (and e.g. toggle) the row's input.
  icon.addEventListener("click", (event) => event.preventDefault());
  return icon;
}
