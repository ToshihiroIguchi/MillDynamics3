import { expect, test } from "@playwright/test";

test("speed badge dot and bar stay put while the badge text changes length", async ({ page }) => {
  await page.goto("/");
  await page.locator(".hud", { hasText: /rpm/ }).waitFor();

  // The badge is right-anchored; with the indicator box (dot/bar) last, a longer or shorter text may only
  // resize the text box leftwards (ui/hud.ts), never move the indicator.
  const spread = await page.evaluate(async () => {
    const badge = document.querySelector(".hud-speed-badge") as HTMLElement;
    const bar = badge.querySelector(".hud-speed-bar") as HTMLElement;
    const dot = badge.querySelector(".hud-speed-dot") as HTMLElement;
    const barLefts: number[] = [];
    const dotLefts: number[] = [];
    const t0 = performance.now();
    while (performance.now() - t0 < 3000) {
      barLefts.push(bar.getBoundingClientRect().left);
      dotLefts.push(dot.getBoundingClientRect().left);
      await new Promise((r) => requestAnimationFrame(r));
    }
    const range = (xs: number[]) => Math.max(...xs) - Math.min(...xs);
    return { bar: range(barLefts), dot: range(dotLefts) };
  });

  expect(spread.bar).toBeLessThan(0.5);
  expect(spread.dot).toBeLessThan(0.5);
});
