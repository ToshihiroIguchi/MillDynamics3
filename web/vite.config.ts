import { defineConfig } from "vite";

export default defineConfig({
  worker: {
    // Dedicated worker (see src/worker.ts) is imported as an ES module.
    format: "es",
  },
  build: {
    target: "es2022",
  },
});
