import { expect, test } from "@playwright/test";
import { openMetricsPanel } from "../helpers/panels";

test("Generate report now downloads a PDF immediately, with no end time required", async ({ page }) => {
  await page.addInitScript(() => localStorage.removeItem("milldynamics.panel"));
  await page.goto("/");
  await openMetricsPanel(page);
  // The Report section is a collapsed <details> by default.
  await page.locator("details.metrics-report > summary").click();

  const downloadPromise = page.waitForEvent("download");
  await page.getByRole("button", { name: "Generate report now" }).click();
  const download = await downloadPromise;

  expect(download.suggestedFilename()).toMatch(/^milldynamics-report-.*\.pdf$/);
});
