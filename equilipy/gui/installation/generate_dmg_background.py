"""Generate the macOS DMG installer background image."""

from __future__ import annotations

from pathlib import Path

from PIL import Image, ImageDraw, ImageFont

ROOT = Path(__file__).resolve().parents[3]
OUTPUT = ROOT / "equilipy" / "gui" / "icons" / "dmg_background.png"

WIDTH = 760
HEIGHT = 480
SCALE = 3


def _font(
    size: int, *, bold: bool = False
) -> ImageFont.FreeTypeFont | ImageFont.ImageFont:
    """Return a system font with a conservative fallback."""
    candidates = [
        "/System/Library/Fonts/Supplemental/Arial Bold.ttf" if bold else "",
        "/System/Library/Fonts/Supplemental/Arial.ttf",
        "/System/Library/Fonts/SFNS.ttf",
    ]
    for candidate in candidates:
        if not candidate:
            continue
        try:
            return ImageFont.truetype(candidate, size * SCALE)
        except OSError:
            continue
    return ImageFont.load_default()


def _xy(value: float) -> int:
    """Scale a coordinate to the internal high-resolution canvas."""
    return round(value * SCALE)


def _draw_background(draw: ImageDraw.ImageDraw, image: Image.Image) -> None:
    """Draw the dark installation-panel background."""
    for y in range(image.height):
        for x in range(image.width):
            nx = x / image.width
            ny = y / image.height
            shade = nx * 0.10 + ny * 0.05
            image.putpixel(
                (x, y),
                (
                    max(15, int(24 - 9 * shade)),
                    max(25, int(36 - 9 * shade)),
                    max(23, int(34 - 9 * shade)),
                    255,
                ),
            )



def _draw_title(draw: ImageDraw.ImageDraw) -> None:
    """Draw the installer title."""
    title = "Drag Equilipy to Applications"
    font = _font(34, bold=True)
    bbox = draw.textbbox((0, 0), title, font=font)
    x = (WIDTH * SCALE - (bbox[2] - bbox[0])) / 2
    draw.text((x, _xy(72)), title, fill=(240, 246, 244, 255), font=font)


def _draw_label_pads(draw: ImageDraw.ImageDraw) -> None:
    """Draw light pads behind Finder labels so black labels stay readable."""
    fill = (219, 227, 225, 226)
    pads = [
        (195, 312, 164, 28),
        (565, 312, 220, 28),
    ]
    for cx, cy, width, height in pads:
        box = (
            _xy(cx - width / 2 + 30),
            _xy(cy - height / 2 + 2),
            _xy(cx + width / 2 - 30),
            _xy(cy + height / 2 - 2),
        )
        draw.rounded_rectangle(box, radius=_xy(14), fill=fill)


def _draw_arrow(draw: ImageDraw.ImageDraw) -> None:
    """Draw a clean drag direction arrow between the app and Applications."""
    center_x = 380
    center_y = 220
    start = (center_x - 50, center_y)
    end = (center_x + 30, center_y)
    head_tip = (center_x + 58, center_y)
    head_upper = (center_x + 30, center_y - 28)
    head_lower = (center_x + 30, center_y + 28)

    # Main arrow: purple chevron with rounded caps at every endpoint.
    arrow_color = (128, 96, 255, 255)
    shaft_width = _xy(10)
    radius = shaft_width // 2
    shaft = (
        _xy(start[0]),
        _xy(start[1]) - radius,
        _xy(end[0]),
        _xy(end[1]) + radius,
    )
    draw.rounded_rectangle(shaft, radius=radius, fill=arrow_color)
    draw.line(
        (
            _xy(head_upper[0]),
            _xy(head_upper[1]),
            _xy(head_tip[0]),
            _xy(head_tip[1]),
            _xy(head_lower[0]),
            _xy(head_lower[1]),
        ),
        fill=arrow_color,
        width=_xy(12),
        joint="curve",
    )

    # Round the chevron's two open ends (PIL's line() has flat caps only).
    cap_radius = _xy(12) // 2
    for cx, cy in (head_upper, head_lower):
        draw.ellipse(
            (
                _xy(cx) - cap_radius,
                _xy(cy) - cap_radius,
                _xy(cx) + cap_radius,
                _xy(cy) + cap_radius,
            ),
            fill=arrow_color,
        )


def main() -> int:
    """Generate the background image."""
    image = Image.new("RGBA", (WIDTH * SCALE, HEIGHT * SCALE), (23, 35, 33, 255))
    draw = ImageDraw.Draw(image)

    _draw_background(draw, image)
    _draw_title(draw)
    _draw_arrow(draw)
    _draw_label_pads(draw)

    image = image.resize((WIDTH, HEIGHT), Image.Resampling.LANCZOS)
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    image.save(OUTPUT)
    print(f"Wrote {OUTPUT} ({WIDTH} x {HEIGHT})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
