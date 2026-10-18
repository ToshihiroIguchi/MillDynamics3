import { execSync } from "node:child_process";
import { defineConfig } from "vite";

/** Short git commit hash for the version stamp in exported outputs; "unknown" outside a checkout. */
function gitCommit(): string {
  try {
    return execSync("git rev-parse --short HEAD", { stdio: ["ignore", "pipe", "ignore"] })
      .toString()
      .trim();
  } catch {
    return "unknown";
  }
}

export default defineConfig(({ command }) => ({
  // GitHub Pages serves this project as https://toshihiroiguchi.github.io/MillDynamics3/, a
  // subpath, so a production build needs asset URLs rooted there; dev/preview keep "/" so local
  // serving is unaffected.
  base: command === "build" ? "/MillDynamics3/" : "/",
  define: {
    __APP_COMMIT__: JSON.stringify(gitCommit() || "unknown"),
  },
  worker: {
    // Dedicated worker (see src/worker.ts) is imported as an ES module.
    format: "es",
  },
  build: {
    target: "es2022",
  },
}));
