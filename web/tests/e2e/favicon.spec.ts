import { test, expect } from "@playwright/test";

test("favicon assets are configured in head and return 200 OK", async ({ page, request }) => {
  await page.goto("/");

  // Check link tags
  const iconLink = page.locator('link[rel="icon"][type="image/x-icon"]');
  await expect(iconLink).toHaveAttribute("href", "/favicon.ico");

  const svgLink = page.locator('link[rel="icon"][type="image/svg+xml"]');
  await expect(svgLink).toHaveAttribute("href", "/favicon.svg");

  const pngLink = page.locator('link[rel="icon"][sizes="32x32"]');
  await expect(pngLink).toHaveAttribute("href", "/favicon-32x32.png");

  // Verify that HTTP GET for the favicons returns 200 OK
  const icoRes = await request.get("/favicon.ico");
  expect(icoRes.status()).toBe(200);
  expect(icoRes.headers()["content-type"]).toContain("image");

  const svgRes = await request.get("/favicon.svg");
  expect(svgRes.status()).toBe(200);

  const pngRes = await request.get("/favicon-32x32.png");
  expect(pngRes.status()).toBe(200);
});
