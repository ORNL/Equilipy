# GUI Code Organization

The GUI should be organized by user-facing domain, not by vague extraction
roles. A developer should be able to guess where to add a feature from the file
name.

## Naming Rules

- Use domain names such as `equilibrium`, `solidification`, `functions`,
  `compounds`, `solutions`, `input_condition`, and `cp_editor`.
- Avoid broad names such as `helper`, `controller`, `manager`, or `utils` for new
  feature files. If temporary compatibility glue is needed during a refactor,
  make that explicit, for example `compat_exports.py`.
- Keep `main_window.py` as the application shell: window chrome, workspace
  switching, and top-level composition only.
- Keep files under 2000 lines when practical and under 3000 lines as a hard
  limit.
- If a feature is shared, name it for the concept it implements. For example,
  Cp/G/H/S editing shared by database Functions and Compounds belongs in
  `database/cp_editor.py`, not in a generic helper file.

## Calculation Package

Calculation files should follow the calculation workflow:

- `input_condition.py`: species, composition, units, balance species, phase
  selection inputs used across equilibrium and solidification modules.
- `module_input.py`: Qt forms for one calculation module.
- `session_overview.py`: calculation session tree and overview pages.
- `session_files.py`: saving and loading calculation sessions/modules.
- `session_serialization.py`: converting sessions/modules to and from portable
  payloads and generated scripts.
- `runner.py`: running calculations, progress, cancellation, and output routing.
- `results.py`: result table/window column selection and display helpers.

Future module-specific behavior should move into clear module files, such as
`equilibrium.py`, `solidification.py`, and `nucleation.py`, when the shared
module-input file becomes too broad.

## Database Package

Database files should follow the database editor workflow:

- `workspace.py`: database sidebar actions, import/export, validation orchestration,
  and render routing.
- `overview.py`: database overview and record summary views.
- `cp_editor.py`: Gibbs/Cp/H/S conversion, Cp range editors, and continuity checks.

As database editing grows, split `workspace.py` into domain files:

- `functions.py`
- `compounds.py`
- `solutions.py`
- `validation.py`
- `tree_actions.py`

The solution editor should use the same lazy strategy as Functions and Compounds:
load quickly, validate the selected/edited object locally, and run full database
validation only from the `Validate` button.
