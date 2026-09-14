import { defineConfig } from "vite";

export default defineConfig(({ command }) => ({
  // GitHub Pages serves this project as https://toshihiroiguchi.github.io/MillDynamics3/, a
  // subpath, so a production build needs asset URLs rooted there; dev/preview keep "/" so local
  // serving is unaffected.
  base: command === "build" ? "/MillDynamics3/" : "/",
  worker: {
    // Dedicated worker (see src/worker.ts) is imported as an ES module.
    format: "es",
  },
  build: {
    target: "es2022",
  },
}));
