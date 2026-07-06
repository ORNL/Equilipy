# GUI

The desktop application bundles calculations and database tools in one
window. Install:

```bash
pip install 'equilipy[gui]'
```

Then launch:

```bash
equilipy.gui
```

`equilipy-gui` (legacy alias) and `python -m equilipy.gui` launch the same
application.

## Standalone installers

Every release also ships self-contained GUI installers (no Python
required), built by CI and attached to the
[latest GitHub release](https://github.com/ORNL/Equilipy/releases/latest):

- macOS Apple Silicon, macOS 15: [Equilipy-macos-15-arm64.dmg](https://github.com/ORNL/Equilipy/releases/latest/download/Equilipy-macos-15-arm64.dmg)
- macOS Apple Silicon, macOS 14: [Equilipy-macos-14-arm64.dmg](https://github.com/ORNL/Equilipy/releases/latest/download/Equilipy-macos-14-arm64.dmg)
- Windows 64-bit: [Equilipy-windows-amd64.exe](https://github.com/ORNL/Equilipy/releases/latest/download/Equilipy-windows-amd64.exe)
- Linux 64-bit: [Equilipy-linux-x86_64](https://github.com/ORNL/Equilipy/releases/latest/download/Equilipy-linux-x86_64)
- Linux ARM 64-bit: [Equilipy-linux-arm64](https://github.com/ORNL/Equilipy/releases/latest/download/Equilipy-linux-arm64)

The header switches between the two workspaces:

- **[Calculation](calculation)** — run equilibrium and solidification
  calculations organized in sessions and modules.
- **[Database](database)** — inspect, edit, validate, and export
  thermochemical databases.

![Equilipy main window with the Calculation/Database workspace switch in the header](images/main_window.png)

```{toctree}
:maxdepth: 2

calculation
database
```
