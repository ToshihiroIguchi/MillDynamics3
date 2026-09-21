import os
import math
import base64
from PIL import Image

repo_dir = r"c:\Users\toshi\python\MillDynamics3"
public_dir = os.path.join(repo_dir, "web", "public")
os.makedirs(public_dir, exist_ok=True)

src_path = r"C:\Users\toshi\.gemini\antigravity-ide\brain\f27656f0-d00a-4645-9633-35481a14ff3d\.user_uploaded\media_1789984809840.png"
img = Image.open(src_path).convert("RGBA")

# Center: (79.0, 76.0), radius: 57.5 in original 155x155
cx, cy = 79.0, 76.0
r = 57.5

# High-resolution rendering canvas: 512x512
scale = 4
large = img.resize((img.width * scale, img.height * scale), Image.Resampling.LANCZOS)
lcx, lcy, lr = cx * scale, cy * scale, r * scale

out_size = 512
out_circle = Image.new("RGBA", (out_size, out_size), (0, 0, 0, 0))

target_r = 250.0  # 6px transparent padding around circumference
ratio = lr / target_r
dark_bg = (54, 67, 83)

for y in range(out_size):
    for x in range(out_size):
        dx = x + 0.5 - 256.0
        dy = y + 0.5 - 256.0
        dist = math.sqrt(dx * dx + dy * dy)
        
        if dist > target_r + 2.0:
            continue
            
        src_x = lcx + dx * ratio
        src_y = lcy + dy * ratio
        
        ix = int(round(src_x))
        iy = int(round(src_y))
        
        if 0 <= ix < large.width and 0 <= iy < large.height:
            p = large.getpixel((ix, iy))
            
            if dist <= target_r - 1.0:
                out_circle.putpixel((x, y), p)
            else:
                coverage = max(0.0, min(1.0, (target_r + 1.0 - dist) / 2.0))
                alpha = int(255 * coverage)
                out_circle.putpixel((x, y), (dark_bg[0], dark_bg[1], dark_bg[2], alpha))

# 1. Master PNGs
out_circle.save(os.path.join(public_dir, "icon-512.png"))

icon_192 = out_circle.resize((192, 192), Image.Resampling.LANCZOS)
icon_192.save(os.path.join(public_dir, "icon-192.png"))

favicon_48 = out_circle.resize((48, 48), Image.Resampling.LANCZOS)
favicon_48.save(os.path.join(public_dir, "favicon-48x48.png"))

favicon_32 = out_circle.resize((32, 32), Image.Resampling.LANCZOS)
favicon_32.save(os.path.join(public_dir, "favicon-32x32.png"))

favicon_16 = out_circle.resize((16, 16), Image.Resampling.LANCZOS)
favicon_16.save(os.path.join(public_dir, "favicon-16x16.png"))

# 2. Apple Touch Icon (180x180)
apple_icon = Image.new("RGBA", (180, 180), (0, 0, 0, 0))
circle_170 = out_circle.resize((170, 170), Image.Resampling.LANCZOS)
apple_icon.paste(circle_170, (5, 5), circle_170)
apple_icon.save(os.path.join(public_dir, "apple-touch-icon.png"))

# 3. Multi-resolution favicon.ico (16, 32, 48, 64)
ico_path = os.path.join(public_dir, "favicon.ico")
out_circle.save(ico_path, format="ICO", sizes=[(16, 16), (32, 32), (48, 48), (64, 64)])
print(f"Generated favicon.ico ({os.path.getsize(ico_path)} bytes)")

# 4. favicon.svg with embedded high-res PNG
with open(os.path.join(public_dir, "icon-512.png"), "rb") as f:
    b64_data = base64.b64encode(f.read()).decode("ascii")

svg_content = f'''<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 512 512" width="100%" height="100%">
  <image width="512" height="512" href="data:image/png;base64,{b64_data}"/>
</svg>
'''
with open(os.path.join(public_dir, "favicon.svg"), "w", encoding="utf-8") as f:
    f.write(svg_content)
print(f"Generated favicon.svg")
