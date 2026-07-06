# Code Standard

This document defines the Python and Fortran coding standards for Equilipy.

---

# Python Coding Standard

Follow the project configuration in `pyproject.toml` and the existing Equilipy
module structure. Prefer explicit scientific code, stable public APIs, and tests
that document numerical behavior.

## Table of Contents

1. [Guiding Principles](#guiding-principles)
2. [Required Tooling](#required-tooling)
3. [Python Version](#python-version)
4. [File Layout](#file-layout)
5. [Imports](#imports)
6. [Naming](#naming)
7. [Functions and Methods](#functions-and-methods)
8. [Classes](#classes)
9. [Type Hints](#type-hints)
10. [Docstrings](#docstrings)
11. [NumPy and Scientific Data](#numpy-and-scientific-data)
12. [Error Handling and Warnings](#error-handling-and-warnings)
13. [I/O and Serialization](#io-and-serialization)
14. [Command-Line Modules](#command-line-modules)
15. [Tests](#tests)
16. [Comments](#comments)
17. [Compatibility and Deprecation](#compatibility-and-deprecation)
18. [Public API Boundaries](#public-api-boundaries)
19. [Checklist for New Python Code](#checklist-for-new-python-code)

---

## Guiding Principles

- Prefer clear, explicit scientific code over clever compact code.
- Keep numerical data shapes, units, and dtypes visible in names, type hints,
  docstrings, or nearby validation.
- Follow the existing module structure and naming vocabulary before introducing
  new abstractions.
- Make public APIs stable and documented; keep implementation helpers private.
- Use automated formatting and linting as the source of truth.

---

## Required Tooling

Use Ruff for both linting and formatting. Match the project settings:

```toml
[tool.ruff]
line-length = 88
lint.select = [
    "F", # Pyflakes
    "B", # flake8-bugbear
    "I", # isort
    "E", # pycodestyle errors
    "D", # pydocstyle
]
lint.extend-ignore = ["D417", "D100"]

[tool.ruff.lint.pydocstyle]
convention = "numpy"
```

The standard formatter is `ruff-format`; do not use Black separately.

**Daily commands:**

```bash
ruff check . --fix --show-fixes
ruff format .
pytest
```

**Recommended pre-commit hooks:**

- Remove trailing whitespace and ensure final newline.
- Check YAML files; prevent accidentally added large files.
- Normalize line endings to LF.
- Run `ruff --fix --show-fixes` and `ruff-format`.

---

## Python Version

Target Python 3.10 or newer. Use modern annotation syntax:

| Old style | Modern style |
|---|---|
| `Optional[X]` | `X \| None` |
| `Union[A, B]` | `A \| B` |
| `List[str]`, `Dict[str, int]` | `list[str]`, `dict[str, int]` |
| `Tuple[int, ...]` | `tuple[int, ...]` |

Also use `typing.Literal`, `typing.TypedDict`, `typing.TypeAlias`, and
`typing.TypeVar` when they make contracts clearer.

Most nontrivial modules should include:

```python
from __future__ import annotations
```

Place it after the module docstring, before all other imports.

---

## File Layout

Order the contents of a Python source file as follows:

1. Module docstring.
2. Project copyright/license header comments, if required.
3. `from __future__ import annotations`.
4. Standard-library imports.
5. Third-party imports.
6. Local package imports.
7. Module-level constants, regular expressions, and `TypeAlias` definitions.
8. `TypedDict` definitions and small data containers.
9. Public classes and functions.
10. Private helpers.

**Package `__init__.py`:** export only the public API and declare `__all__`.
Add a `py.typed` marker file at the package root to signal PEP 561 type
support.

```python
# equilipy/__init__.py
from equilipy.core import PhaseAssemblage
from equilipy.version import __version__

__all__ = ["PhaseAssemblage", "__version__"]
```

**Minimal source module example:**

```python
"""Utilities for phase assemblage post-processing."""

from __future__ import annotations

from collections.abc import Sequence
from typing import Literal, TypedDict

import numpy as np
from numpy.typing import NDArray

from equilipy.structure import PhaseRecord


class PhaseResult(TypedDict):
    """Serializable phase result."""

    names: list[str]
    amounts: NDArray[np.double]
```

---

## Imports

Ruff/isort owns import ordering. Rules to follow:

- Use **absolute imports** for project modules.
- Prefer `collections.abc` for abstract types: `Callable`, `Iterator`,
  `Sequence`.
- Import `NDArray` from `numpy.typing` for NumPy array annotations.
- Use `import numpy as np`.
- Keep optional heavyweight imports **local to the function** that needs them.
- Use `typing.cast` or `assert isinstance(...)` to document type boundaries.
- **No wildcard imports.**

```python
from collections.abc import Sequence
from typing import Literal, cast

import numpy as np
from numpy.typing import NDArray

from equilipy.core.phase import Phase
```

---

## Naming

| Kind | Convention | Example |
|---|---|---|
| Modules / packages | `snake_case` | `phase_record.py` |
| Functions / methods | `snake_case` | `get_phase_amounts` |
| Variables / parameters | `snake_case` | `cutoff_energy` |
| Classes | `CapWords` | `PhaseAssemblage` |
| Custom exceptions | `CapWordsError` | `ConvergenceError` |
| Constants | `UPPER_SNAKE_CASE` | `MAX_ITERATIONS` |
| Private helpers | `_snake_case` | `_trim_cell` |
| Private attributes | `_snake_case` | `self._symprec` |
| Boolean parameters | `is_/has_/use_` prefix | `is_converged`, `use_cache` |

Scientific abbreviations may keep established domain capitalization when it
aids recognition (`NAC`, `QHA`, `VCA`, `SNF`). Do not invent mixed-case names
for ordinary variables.

---

## Functions and Methods

Annotate all public and nontrivial functions. Keep signatures explicit.

```python
def get_phase_amounts(
    names: Sequence[str],
    amounts: Sequence[float] | NDArray[np.double],
    cutoff: float = 1e-12,
) -> dict[str, float]:
    """Return phase amounts above cutoff."""
    values = np.array(amounts, dtype="double", order="C")
    return {name: float(v) for name, v in zip(names, values) if v > cutoff}
```

**Conventions:**

- Prefer keyword arguments for optional behavior; keep defaults simple and
  immutable.
- Return `None` explicitly only when `None` is meaningful to callers.
- Validate invalid inputs early; build long error messages in a `msg` variable.
- Use `"\n".join([...])` for multiline messages; f-strings for interpolation.
- Use `%` formatting only when matching an existing fixed-width numerical
  output style.
- Do not add comments that restate the code.
- Add a short comment before a non-obvious numerical transformation.

---

## Classes

Use classes for stateful domain objects with clear invariants. Keep the public
surface small.

```python
class PhaseAssemblage:
    """Represent a calculated phase assemblage."""

    def __init__(
        self,
        names: Sequence[str],
        amounts: Sequence[float] | NDArray[np.double],
    ) -> None:
        """Initialize phase assemblage."""
        self._names = list(names)
        self._amounts = np.array(amounts, dtype="double", order="C")
        self._check()

    @property
    def amounts(self) -> NDArray[np.double]:
        """Return phase amounts."""
        return self._amounts.copy()

    def _check(self) -> None:
        if len(self._names) != len(self._amounts):
            raise ValueError("names and amounts have different lengths.")
```

**Conventions:**

- Store state in private attributes; expose via `@property`.
- Return **copies** of mutable arrays when external mutation would break
  invariants.
- Keep validation in small private methods such as `_check`.
- Use `dataclasses.dataclass(frozen=True)` for small immutable value objects.
- Use `@dataclass` (mutable) for configuration / parameter containers.
- Use `TypedDict` for dictionaries passed across API or serialization
  boundaries; support optional fields with inheritance and `total=False`.
- Define custom exceptions in a dedicated `exceptions.py` module; subclass
  `RuntimeError` or `ValueError` with a one-sentence docstring.

**TypedDict with optional fields:**

```python
class DisplacementEntry(TypedDict):
    """Required displacement fields."""
    number: int
    displacement: NDArray[np.double]

class DisplacementEntryWithForces(DisplacementEntry, total=False):
    """Displacement entry optionally containing forces."""
    forces: NDArray[np.double]
    supercell_energy: float
```

**Abstract base classes:**

Use `abc.ABC` and `@abstractmethod` to define interfaces that multiple
implementations share.

```python
from abc import ABC, abstractmethod

class SolverBase(ABC):
    """Abstract base for equilibrium solvers."""

    @abstractmethod
    def solve(self, composition: NDArray[np.double]) -> PhaseAssemblage:
        """Solve for the equilibrium assemblage."""
```

**Exceptions module (`exceptions.py`):**

```python
class ConvergenceError(RuntimeError):
    """Exception when minimization fails to converge."""

class AssemblageNotFoundError(RuntimeError):
    """Exception when no valid phase assemblage is found."""
```

---

## Type Hints

Type hints are part of the coding style, not decorative comments.

- Annotate all public functions, methods, and properties.
- Annotate private helpers when inputs or outputs are not obvious.
- Use `Sequence[T]` for read-only inputs that may be list, tuple, or
  array-like.
- Use `NDArray[np.double]`, `NDArray[np.int64]`, or `NDArray[np.cdouble]` when
  the dtype matters.
- Use `Literal[...]` for constrained string **or integer** modes.
- Use `TypedDict` or a small dataclass over unstructured `dict`.
- Use `TypeVar` for generic utility functions.
- Use `Any` sparingly at boundaries to untyped libraries or legacy data.
- Use `cast(...)` only after runtime validation or when a library type is
  incomplete.
- Use `# type: ignore[<code>]` only at unavoidable library typing gaps; include
  the narrowest ignore code possible.

```python
_T = TypeVar("_T")

def first_or_default(items: Sequence[_T], default: _T) -> _T:
    """Return the first item or default if empty."""
    return items[0] if items else default


class ForceDataset(TypedDict):
    """Forces and displacements for one calculation set."""

    displacements: NDArray[np.double]
    forces: NDArray[np.double]
    energies: NDArray[np.double] | None
```

---

## Docstrings

Use **NumPy-style** docstrings (`convention = "numpy"` in Ruff).

Every public module, class, function, method, and property needs a docstring.
Private helpers need one when their purpose or return contract is not obvious.

**Short function:**

```python
def get_number_of_phases(amounts: NDArray[np.double]) -> int:
    """Return number of nonzero phases."""
    return int(np.count_nonzero(amounts))
```

**Public function with parameters:**

```python
def normalize_amounts(
    amounts: Sequence[float] | NDArray[np.double],
    total: float = 1.0,
) -> NDArray[np.double]:
    """Normalize phase amounts.

    Parameters
    ----------
    amounts : array_like
        Phase amounts.
    total : float, optional
        Target sum after normalization. Default is 1.0.

    Returns
    -------
    ndarray
        Normalized phase amounts. shape=(n_phase,), dtype='double'.

    """
    values = np.array(amounts, dtype="double", order="C")
    return values / values.sum() * total
```

**Class with attributes:**

```python
class PhaseGrid:
    """Grid of phase amounts.

    Attributes
    ----------
    points : ndarray
        Grid coordinates. shape=(n_point, n_component), dtype='double'.
    amounts : ndarray
        Phase amounts. shape=(n_point, n_phase), dtype='double'.

    """
```

**Conventions:**

- Start with a short imperative or descriptive sentence.
- Use `Parameters`, `Returns`, `Attributes`, `Raises`, `Notes`, and `Examples`
  sections when useful.
- Document **shapes, dtypes, and units** for all numerical data.
- Use `array_like` when the implementation accepts sequences and converts them.
- Keep deprecation notes in the parameter or method docstring; emit a
  `DeprecationWarning` in code.
- Do not document every obvious private implementation detail.

---

## NumPy and Scientific Data

- Convert user-provided numerical inputs with `np.array(...)` or
  `np.asarray(...)` near the boundary.
- Specify `dtype="double"`, `dtype="int64"`, or another exact dtype when the
  downstream calculation depends on it.
- Specify `order="C"` when C-contiguous memory is required.
- Keep shape checks close to the conversion point.
- Prefer vectorized NumPy operations over Python loops for numerical kernels.
- Use `np.rint`, `np.linalg`, `np.dot`, `np.c_`, and similar standard APIs
  directly when they express the calculation clearly.
- Avoid hidden unit conversions; name the unit in the variable, docstring, or
  a conversion helper.

```python
def get_distances(
    scaled_positions: Sequence[Sequence[float]] | NDArray[np.double],
    lattice: Sequence[Sequence[float]] | NDArray[np.double],
) -> NDArray[np.double]:
    """Return Cartesian distances from fractional coordinates."""
    positions = np.array(scaled_positions, dtype="double", order="C")
    cell = np.array(lattice, dtype="double", order="C")
    return np.linalg.norm(np.dot(positions, cell), axis=1)
```

**In tests:** use `np.testing.assert_allclose` for floating arrays and
`np.testing.assert_array_equal` for exact arrays.

---

## Error Handling and Warnings

Use explicit exceptions with actionable messages.

```python
if len(names) != len(amounts):
    raise ValueError("names and amounts have different lengths.")
```

For multiline or contextual errors:

```python
msg = "\n".join(
    [
        "Phase assemblage creation failed.",
        "Probably some phase amounts are negative.",
        str(amounts),
    ]
)
raise RuntimeError(msg)
```

**Conventions:**

| Situation | Exception / warning |
|---|---|
| Invalid argument value | `ValueError` |
| Failed calculation / inconsistent state | `RuntimeError` |
| Common domain failure callers may catch | Custom exception in `exceptions.py` |
| Chained exception | `raise NewError(msg) from exc` |
| Deprecated API | `warnings.warn(..., DeprecationWarning, stacklevel=2)` |
| Nonfatal scientific condition | `warnings.warn(..., UserWarning, stacklevel=2)` |

- Define custom exceptions in a dedicated `exceptions.py` module.
- Avoid bare `except`; always catch specific exception classes.

---

## I/O and Serialization

Keep parsing and formatting separated from numerical computation.

- Accept `str | os.PathLike` for filenames; use `pathlib.Path` in new code.
- Build text output as `list[str]` and join with `"\n".join(lines)` for
  testability.
- Validate required fields before constructing domain objects in parsers.
- Use `yaml.safe_load` for YAML input.
- Keep file-format quirks inside interface modules.

```python
def get_phase_yaml_lines(amounts: NDArray[np.double]) -> list[str]:
    """Return phase amounts as YAML lines."""
    lines = ["phase_amounts:"]
    for i, amount in enumerate(amounts):
        lines.append(f"- index: {i}")
        lines.append(f"  amount: {amount:.16f}")
    return lines
```

---

## Command-Line Modules

Use thin script entry points; put real logic in importable modules.

```python
from equilipy.cui.main_script import main


def run():
    """Run equilipy script."""
    main()
```

Command-line code may print status, warnings, and banners. Core library code
should return data or raise exceptions instead of printing.

---

## Tests

Use pytest. Test files live under `test/` or `tests/` and are named
`test_*.py`. Shared fixtures go in `conftest.py`.

```python
# conftest.py
@pytest.fixture(scope="session")
def phase_assemblage() -> PhaseAssemblage:
    """Return a prebuilt PhaseAssemblage for reuse across tests."""
    return PhaseAssemblage(names=["Al", "Fe"], amounts=[0.3, 0.7])
```

**Conventions:**

- Test functions start with `test_`.
- Use `@pytest.fixture(scope="session")` for expensive reusable objects.
- Use `@pytest.mark.parametrize` for input/output matrices.
- Use plain `assert` for scalars and exact strings.
- Use `pytest.approx` for scalar floating-point comparisons.
- Use `np.testing.assert_allclose` / `assert_array_equal` for arrays.
- Use `pytest.raises` to test expected errors.
- Keep regression data small and clearly named.

```python
@pytest.mark.parametrize("scale", [1.0, 2.0])
def test_normalize_amounts(scale: float):
    """Test normalize_amounts."""
    amounts = normalize_amounts([1.0, 3.0], total=scale)
    np.testing.assert_allclose(amounts, [0.25 * scale, 0.75 * scale])


def test_normalize_amounts_zero_sum():
    """Test normalize_amounts raises on zero-sum input."""
    with pytest.raises(ValueError, match="sum"):
        normalize_amounts([0.0, 0.0])
```

---

## Comments

Comments explain **why** a non-obvious step exists, never **what** it does.

```python
# Good: explains the physical reason for the operation
# Fold fractional coordinates back into the first periodic image before
# computing Cartesian distances.
diff -= np.rint(diff)

# Bad: just restates the code
# Subtract rounded diff from diff.
diff -= np.rint(diff)
```

---

## Compatibility and Deprecation

Keep compatibility wrappers small and explicit:

```python
def OldPhase(*args, **kwargs) -> "Phase":
    """Deprecated alias of Phase."""
    warnings.warn(
        "OldPhase is deprecated. Please use Phase instead.",
        DeprecationWarning,
        stacklevel=2,
    )
    return Phase(*args, **kwargs)
```

- Prefer warnings over silent behavior changes.
- Include the replacement API in the warning message.
- Keep deprecated wrappers out of new code paths.

---

## Public API Boundaries

- Public functions must have stable parameter names.
- Public return dictionaries should be `TypedDict` when feasible.
- Public classes should expose data through properties or documented methods.
- Avoid returning internal mutable state directly.
- Avoid changing default behavior without a deprecation period.

---

## Checklist for New Python Code

Before submitting Python changes:

- [ ] `ruff format .` passes cleanly.
- [ ] `ruff check . --fix --show-fixes` passes cleanly.
- [ ] Relevant `pytest` tests pass.
- [ ] Every public object has a NumPy-style docstring.
- [ ] Array shapes, dtypes, and units are documented where relevant.
- [ ] Optional dependencies are imported lazily or declared properly.
- [ ] New tests cover both successful behavior and important failure cases.
- [ ] No unrelated formatting churn is included.

---

# Fortran Coding Standard

Covers the free-form Fortran 90 code in `equilipy/fsrc`. New project code
follows the Thermochimica-style procedural structure used in `1_Modules`,
`2_Utils`, `3_SystemDef`, `4_GibbsModels`, `5_Submin`,
`6_ActiveSetLagrangian`, and `7_MultiPhaseMinimizer`. LAPACK and BLAS drivers
are linked from the platform numerical library during the f2py build; the old
Netlib source copy is archived under `archive/fsrcArchives/0Lapack` for
provenance only.

## Table of Contents

1. [Source File Policy](#source-file-policy)
2. [File Layout](#file-layout-1)
3. [Indentation and Formatting](#indentation-and-formatting)
4. [Keyword and Operator Style](#keyword-and-operator-style)
5. [Naming](#naming-1)
6. [Types and Declarations](#types-and-declarations)
7. [Modules and Shared State](#modules-and-shared-state)
8. [Control Flow](#control-flow)
9. [Numerical Code](#numerical-code)
10. [LAPACK and Linear Solves](#lapack-and-linear-solves)
11. [Allocation and Cleanup](#allocation-and-cleanup)
12. [Error Handling and Diagnostics](#error-handling-and-diagnostics)
13. [Comments and Documentation](#comments-and-documentation)
14. [Compatibility With Legacy Code](#compatibility-with-legacy-code)
15. [Checklist for New Fortran Code](#checklist-for-new-fortran-code)

---

## Source File Policy

| Source area | Style rule |
|---|---|
| `1_Modules`, `3_SystemDef`, `4_GibbsModels`, `5_Submin`, `6_ActiveSetLagrangian`, `7_MultiPhaseMinimizer` | Follow this standard. |
| `2_Utils` project utilities | Follow this standard; preserve third-party attribution and algorithm-specific structure. |
| `archive/fsrcArchives/0Lapack/*.f` | Archived Netlib LAPACK/BLAS source; do not restyle or compile. |

---

## File Layout

Use one primary module or subroutine per `.f90` file when practical. The file
name should match the primary routine or module name.

**Subroutine file order:**

1. Optional leading `!` separator lines.
2. `subroutine Name(arguments)`.
3. Doxygen-style comment block.
4. `USE` statements.
5. `implicit none`.
6. Dummy argument declarations.
7. Local variable declarations.
8. Initialization block.
9. Algorithm body.
10. Cleanup, `return`, and `end subroutine Name`.

```fortran
!
!
subroutine CompExample(iSolnPhaseIndex, lPhasePass)
!
    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    CompExample.f90
    !> \brief   Compute an example thermodynamic quantity.
    !> \author  Initials. Lastname
    !> \date    May 4, 2026
    !
    ! Purpose:
    ! ========
    !
    !> \details The purpose of this subroutine is to demonstrate the local
    !! Fortran layout.
    !
    ! Pertinent variables:
    ! ====================
    !
    !> \param[in]  iSolnPhaseIndex  Absolute index of the solution phase.
    !> \param[out] lPhasePass       Whether the phase passes the test.
    !
    !-------------------------------------------------------------------------------------------------------------
    USE ModuleThermo
    USE ModuleGEMSolver
!
    implicit none
!
    integer, intent(in)  :: iSolnPhaseIndex
    logical, intent(out) :: lPhasePass
    integer              :: i, iFirst, iLast
    real(8)              :: dTemp
!
    ! Initialize variables:
    lPhasePass = .FALSE.
    dTemp      = 0D0
!
    iFirst = nSpeciesPhase(iSolnPhaseIndex-1) + 1
    iLast  = nSpeciesPhase(iSolnPhaseIndex)
!
    do i = iFirst, iLast
        dTemp = dTemp + dMolFraction(i)
    end do
!
    if (dTemp > 1D0 - dTolerance(1)) lPhasePass = .TRUE.
!
    return
!
end subroutine CompExample
!
!
```

**Module file order:**

1. `module ModuleName`.
2. Doxygen-style module comment block.
3. `implicit none`.
4. `SAVE` when the module stores global state.
5. Parameters, scalar state, allocatable arrays, character data, logical data.
6. `end module ModuleName`.

---

## Indentation and Formatting

- Use **free-form `.f90`** for project code.
- Indent executable code and declarations with **4 spaces**; no tabs.
- Keep `USE` statements before `implicit none`.
- Keep declaration groups visually aligned when it improves scanning.
- Use `!` separator lines between declaration, initialization, algorithm, and
  cleanup blocks when matching the surrounding file.
- Prefer readable code over dense one-line expressions.
- Keep new executable lines below about 120 characters; use `&` continuation
  for long expressions. Align the continuation under the expression.

```fortran
dChemicalPotential(i) = dStdGibbsEnergy(i) &
                      + DLOG(DMAX1(dMolFraction(i), 1D-75)) &
                      + dPartialExcessGibbs(i)
```

Do not reformat whole legacy files just to satisfy these rules. Apply the
standard to new code and to touched blocks only.

---

## Keyword and Operator Style

| Token | Style | Examples |
|---|---|---|
| `USE` | uppercase | `USE ModuleThermo` |
| `implicit none` | lowercase | `implicit none` |
| Control keywords | lowercase | `if`, `then`, `else`, `end if`, `do`, `end do`, `select case`, `case`, `return`, `call` |
| `PARAMETER` | uppercase | `INTEGER, PARAMETER :: iTolNum = 15` |
| `SAVE` | uppercase | `SAVE` |
| `recursive` | lowercase | `recursive subroutine Qsort(...)` |
| Logical operators | uppercase | `.AND.`, `.OR.`, `.NOT.` |
| Logical constants | uppercase | `.TRUE.`, `.FALSE.` |
| Double literals | `D` exponent | `0D0`, `1D0`, `1D-8` |
| Legacy double intrinsics | uppercase | `DABS`, `DLOG`, `DMAX1`, `DMIN1`, `DFLOAT` |

Use legacy double intrinsics consistently with nearby code; do not mix them
with their generic equivalents in the same block.

---

## Naming

Use the existing Thermochimica-style prefixes. They make type and meaning
visible in numerical code.

| Kind | Convention | Examples |
|---|---|---|
| Modules | `Module` + `CapWords` | `ModuleThermo`, `ModuleGEMSolver` |
| Subroutines | `CapWords` | `CompThermoData`, `GEMNewton` |
| Integer counters / indices | `i`, `j`, `k`; descriptive `i...` | `iSolnIndex` |
| Integer size / count | `n...` | `nElements`, `nSpecies` |
| Reals (double precision) | `d...` | `dMolesPhase`, `dTolerance` |
| Characters | `c...` | `cSolnPhaseName` |
| Logicals | `l...` | `lConverged`, `lPhasePass` |
| Iteration counters | `iter...` | `iterGlobal`, `iterSub` |
| LAPACK status | `INFO` / `INFOLocal` | `INFO`, `INFOLocal` |
| Project error code | `INFOThermo` | `INFOThermo` |

Use established domain acronyms exactly as the surrounding code does: `GEM`,
`PEA`, `SUBL`, `SUBG`, `QKTO`, `RKMP`, `SRO`, `NAC`.

---

## Types and Declarations

Project thermodynamic code uses `real(8)` and `D` literals. Continue that
style for consistency.

```fortran
integer, intent(in)                :: iSolnPhaseIndex
integer                            :: i, j, iFirst, iLast
integer, dimension(:), allocatable :: iIndex
real(8), intent(out)               :: dDrivingForce
real(8), dimension(:), allocatable :: dMolFractionTemp
logical                            :: lPhasePass
character(25)                      :: cSolnPhaseNameTemp
```

Conventions:

- Always use `implicit none`.
- Add `intent(in)`, `intent(out)`, or `intent(inout)` to all dummy arguments
  in new and touched routines.
- Use `dimension(...)` declarations to match the project style.
- Use `allocatable` for arrays whose size depends on parsed data or system
  size.
- Use module parameters for shared constants.
- Keep local scalar declarations before local array declarations when
  practical.
- Prefer `character(n)` in new code. Preserve existing `character*n` syntax in
  untouched legacy declarations.
- Use `real(8)` consistently in active Fortran source. Do not add a local
  precision module just to alias double precision.

---

## Modules and Shared State

Global thermodynamic state is stored in modules under `1_Modules`.

- Put shared state in a module, not in ad hoc `COMMON` blocks.
- Keep module names descriptive and stable.
- Use `SAVE` in state modules to match existing modules.
- Group variables by type and rank: integers, reals, characters, logicals.
- Use `USE ModuleName, ONLY: ...` when importing a small, clear subset.
- Use plain `USE ModuleName` when a routine depends broadly on that module, as
  the solver routines commonly do.

Do not introduce new global module variables unless the data is genuinely
shared state. Prefer local variables for temporary work arrays.

---

## Control Flow

Use structured control flow.

- Use `select case` for solution phase types, unit names, and error-code
  dispatch.
- Use **named loops** for long loops, nested loops, and loops that use `exit`
  or `cycle`. End the loop with the same label.
- Use **named `if` blocks** for long or deeply nested conditionals where naming
  the block improves readability or allows `exit`. End the block with the same
  label.
- Use early `return` after setting `INFOThermo` for unrecoverable errors.
- Avoid new `goto` in project code. Preserve it only in translated or vendored
  algorithms such as `nnls`.

```fortran
LOOP_SolnPhaseSys: do i = 1, nSolnPhasesSys
    if (lSolnPhases(i)) cycle LOOP_SolnPhaseSys

    select case (cSolnPhaseType(i))
        case ('SUBG', 'SUBQ')
            call CompExcessGibbsEnergySUBG(i)
        case default
            INFOThermo = 17
            return
    end select
end do LOOP_SolnPhaseSys

IF_ParamPass: if (iParamPassCS(j) /= 0) then
    ! Process parameter only when it passes the filter.
    ...
end if IF_ParamPass
```

---

## Numerical Code

Make numerical assumptions visible.

- Initialize scalars and arrays before use: `i = 0`, `A = 0D0`,
  `lPhasePass = .FALSE.`.
- Use project tolerances (`dTolerance`, `dPEATol`, `dMaxPotentialTol`) instead
  of scattering new magic thresholds. When a new threshold is needed, name it
  and document why it is safe.
- Use array operations (`sum`, `MATMUL`, `dot_product`, `MINVAL`, `MAXVAL`)
  when they improve clarity.
- Use 1-based indexing and document derived ranges such as first and last
  species in a phase.
- Keep units clear in comments, especially where dimensionless Gibbs energies,
  J/mol, mole fractions, atom fractions, and gram fractions meet.

```fortran
m = nSpeciesPhase(k-1) + 1      ! First constituent in phase.
n = nSpeciesPhase(k)            ! Last  constituent in phase.
dTemp = sum(dStoichSpecies(m:n,j) * dMolesSpecies(m:n))
```

---

## LAPACK and Linear Solves

The code uses LAPACK-style status handling.

- Active builds link LAPACK/BLAS externally through `setup.py`; do not add new
  Netlib routines back under `equilipy/fsrc`.
- Name pivot arrays `IPIV`.
- Name status variables `INFO` or `INFOLocal`.
- Initialize solver inputs and status before the call.
- Check `INFO` immediately after the call.
- Set `INFOThermo` when a solver failure should propagate to the caller.
- Keep LAPACK call signatures explicit; do not hide important dimensions.

```fortran
INFO = 0
IPIV = 0
call DGESV(nElements, 1, A, nElements, IPIV, dMolesPhase, nElements, INFO)

if (INFO /= 0) then
    INFOThermo = 10
    return
end if
```

---

## Allocation and Cleanup

Memory management is explicit.

- Check `allocated(...)` before deallocating reusable arrays.
- Allocate work arrays after dimensions are known and close to first use.
- Initialize allocated arrays immediately.
- Deallocate local allocatable work arrays before returning.
- For cleanup routines, continue the project pattern of testing each module
  allocation and deallocating it explicitly.
- Use `STAT` when deallocation failure needs to set `INFOThermo`.

```fortran
if (allocated(dTempVec)) deallocate(dTempVec)
allocate(dTempVec(nSolnPhases))
dTempVec = 0D0
...
deallocate(dTempVec, STAT = INFO)
if (INFO /= 0) INFOThermo = 24
```

---

## Error Handling and Diagnostics

Use the project error-code mechanism.

- Use `INFOThermo` for project-level error propagation.
- Use local `INFO` for LAPACK and allocation/deallocation status.
- Do not add `stop` or `error stop` in library routines.
- Use `print *` only in user-facing output routines (`PrintResults`),
  diagnostic routines (`ThermoDEBUG`), or under `lDebugMode`.
- For NaN checks, preserve the local idiom `x /= x` unless replacing a whole
  block with a clearer checked helper.

```fortran
if (dTemperature /= dTemperature) then
    INFOThermo = 1
    return
end if
```

---

## Comments and Documentation

Use Doxygen-style comments for modules and nontrivial subroutines.

**Required header fields for new substantial routines:**

- `\file`, `\brief`, `\author`, `\date`
- `Purpose` section
- `Pertinent variables` section
- `\param[in]`, `\param[out]`, or `\param[in,out]` for each dummy argument

Use inline comments to explain physical meaning, algorithmic steps, units, or
non-obvious numerical safeguards.

```fortran
! Good: explains the physical reason
! Normalize per gram-atom before comparing phase potentials.
dTemp = dTemp / dSpeciesTotalAtoms(i)

! Bad: restates the code
! Divide dTemp by dSpeciesTotalAtoms.
dTemp = dTemp / dSpeciesTotalAtoms(i)
```

Remove temporary debug prints before committing, or guard them with
`lDebugMode`.

---

## Compatibility With Legacy Code

This codebase has a long numerical history. Do not churn stable legacy style
without a behavioral reason.

- Preserve archived fixed-form LAPACK files in `archive/fsrcArchives/0Lapack`
  for provenance only.
- Preserve third-party attribution in utility routines.
- Keep existing public subroutine and module names stable.
- Keep established error codes stable unless the caller contract is updated.
- When modernizing a touched block, keep the change local and verify
  numerical behavior.

---

## Checklist for New Fortran Code

Before submitting Fortran changes:

- [ ] The file is free-form `.f90` unless it is vendored fixed-form code.
- [ ] The primary file name matches the module or subroutine name.
- [ ] `USE` statements appear before `implicit none`.
- [ ] All dummy arguments have `intent(...)`.
- [ ] New real thermodynamic variables use `real(8)` and `D` literals.
- [ ] Arrays are initialized before use.
- [ ] Allocatable arrays are deallocated on all relevant paths.
- [ ] LAPACK `INFO` values and project `INFOThermo` errors are checked.
- [ ] Long or nested loops that use `exit` or `cycle` are named.
- [ ] New numerical tolerances are named or clearly documented.
- [ ] Debug printing is removed or guarded by `lDebugMode`.
