import { expect, test } from "@playwright/test";
import { openMetricsPanel } from "../helpers/panels";

test("Export CSV downloads a metrics CSV file once history has accumulated", async ({ page }) => {
  await page.addInitScript(() => localStorage.removeItem("milldynamics.panel"));
  await page.goto("/");
  await openMetricsPanel(page);

  // Let history.length grow past 0 so the Export CSV button becomes enabled.
  await page.waitForTimeout(3000);

  const downloadPromise = page.waitForEvent("download");
  await page.getByRole("button", { name: "Export CSV" }).click();
  const download = await downloadPromise;

  expect(download.suggestedFilename()).toMatch(/^milldynamics-metrics-.*\.csv$/);
});
