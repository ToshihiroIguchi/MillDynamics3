import { expect, test } from "@playwright/test";

test("Generate report now downloads a PDF immediately, with no end time required", async ({ page }) => {
  await page.addInitScript(() => localStorage.removeItem("milldynamics.panel"));
  await page.goto("/");

  const downloadPromise = page.waitForEvent("download");
  await page.getByRole("button", { name: "Generate report now" }).click();
  const download = await downloadPromise;

  expect(download.suggestedFilename()).toMatch(/^milldynamics-report-.*\.pdf$/);
});
