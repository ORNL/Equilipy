"""dmgbuild settings for the Equilipy macOS installer."""

from __future__ import annotations

from pathlib import Path

ROOT = Path.cwd()

application = "Equilipy"
format = "UDZO"
size = "1200M"

files = [str(ROOT / "pyinstaller-dist" / "Equilipy.app")]
symlinks = {"Applications": "/Applications"}
hide = [".background.png", ".VolumeIcon.icns"]

icon = str(ROOT / "equilipy" / "gui" / "icons" / "equilipy.icns")
background = str(ROOT / "equilipy" / "gui" / "icons" / "dmg_background.png")

default_view = "icon-view"
show_icon_preview = False
include_icon_view_settings = "auto"
arrange_by = None
grid_offset = (0, 0)
grid_spacing = 100
scroll_position = (0, 0)
label_pos = "bottom"
text_size = 13
icon_size = 144

window_rect = ((200, 120), (760, 480))
show_status_bar = False
show_tab_view = False
show_toolbar = False
show_pathbar = False
show_sidebar = False
sidebar_width = 0

icon_locations = {
    "Equilipy.app": (195, 220),
    "Applications": (565, 220),
    ".background.png": (5000, 5000),
    ".VolumeIcon.icns": (5100, 5000),
}
