// Software version string for outputs (PDF report, CSV, batch exports). `__APP_COMMIT__` is a Vite
// `define` (vite.config.ts, short git commit hash); it is undefined under plain vitest, so guard it.

declare const __APP_COMMIT__: string;

export const APP_COMMIT: string = typeof __APP_COMMIT__ === "string" ? __APP_COMMIT__ : "unknown";

export function versionLine(): string {
  return `MillDynamics3 commit ${APP_COMMIT}`;
}
