"""Active-component bundle projection helpers for minimizer development."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class ActiveComponentProjection:
    """Projection of one full-element candidate onto active-component bundles."""

    coefficients: np.ndarray
    reconstructed_stoichiometry: np.ndarray
    residual: np.ndarray
    residual_norm: float
    max_abs_residual: float
    rank: int
    compatible: bool
    nonnegative: bool


@dataclass(frozen=True)
class ActiveComponentBasis:
    """Active-component bundles expressed in full element space."""

    element_names: tuple[str, ...]
    component_names: tuple[str, ...]
    component_matrix: np.ndarray

    @classmethod
    def from_stoichiometries(
        cls,
        components: Mapping[str, Mapping[str, float] | Sequence[float]],
        element_names: Sequence[str],
    ) -> "ActiveComponentBasis":
        """Build an active-component basis from component stoichiometries.

        The matrix convention follows the Step 7 design note: columns are
        active-component bundles in full element space, so a candidate is
        compatible when ``nu = C @ a`` within tolerance.
        """
        normalized_elements = _normalize_names(element_names)
        component_names = tuple(str(name) for name in components)
        if not component_names:
            raise ValueError("At least one active component is required.")

        columns = [
            stoichiometry_vector(stoichiometry, normalized_elements)
            for stoichiometry in components.values()
        ]
        matrix = np.column_stack(columns)
        return cls(
            element_names=normalized_elements,
            component_names=component_names,
            component_matrix=matrix,
        )

    def project(
        self,
        candidate_stoichiometry: Mapping[str, float] | Sequence[float],
        *,
        require_nonnegative: bool = True,
        atol: float = 1e-10,
        rtol: float = 1e-10,
    ) -> ActiveComponentProjection:
        """Project a candidate full stoichiometry onto active components."""
        candidate = stoichiometry_vector(candidate_stoichiometry, self.element_names)
        matrix = np.asarray(self.component_matrix, dtype=float)
        if require_nonnegative:
            coefficients = _nnls(matrix, candidate, tolerance=max(atol, rtol))
        else:
            coefficients = np.linalg.lstsq(matrix, candidate, rcond=None)[0]

        reconstructed = matrix @ coefficients
        residual = candidate - reconstructed
        residual_norm = float(np.linalg.norm(residual))
        max_abs_residual = float(np.max(np.abs(residual))) if residual.size else 0.0
        scale = max(1.0, float(np.linalg.norm(candidate)))
        nonnegative = bool(np.all(coefficients >= -atol))
        compatible = bool(
            residual_norm <= atol + rtol * scale
            and max_abs_residual
            <= atol + rtol * max(1.0, float(np.max(np.abs(candidate))))
            and (nonnegative or not require_nonnegative)
        )
        return ActiveComponentProjection(
            coefficients=coefficients,
            reconstructed_stoichiometry=reconstructed,
            residual=residual,
            residual_norm=residual_norm,
            max_abs_residual=max_abs_residual,
            rank=int(np.linalg.matrix_rank(matrix)),
            compatible=compatible,
            nonnegative=nonnegative,
        )

    def driving_force(
        self,
        candidate_gibbs: float,
        coefficients: Sequence[float],
        component_potentials: Mapping[str, float] | Sequence[float],
    ) -> float:
        """Evaluate driving force against the active-component plane."""
        potentials = self._potential_vector(component_potentials)
        return float(
            candidate_gibbs
            - np.dot(np.asarray(coefficients, dtype=float), potentials)
        )

    def representative_element_potentials(
        self,
        component_potentials: Mapping[str, float] | Sequence[float],
    ) -> np.ndarray:
        """Return a minimum-norm full-element potential vector for diagnostics.

        This solves ``C.T @ mu_e = mu_active``.  For rank-deficient full element
        spaces this vector is only a gauge choice; driving-force decisions
        should use active-component coefficients directly.
        """
        potentials = self._potential_vector(component_potentials)
        return np.linalg.lstsq(self.component_matrix.T, potentials, rcond=None)[0]

    def _potential_vector(
        self,
        component_potentials: Mapping[str, float] | Sequence[float],
    ) -> np.ndarray:
        if isinstance(component_potentials, Mapping):
            return np.asarray(
                [float(component_potentials[name]) for name in self.component_names],
                dtype=float,
            )
        values = np.asarray(component_potentials, dtype=float)
        if values.shape != (len(self.component_names),):
            raise ValueError(
                "component potential vector length must match active components."
            )
        return values


def stoichiometry_vector(
    stoichiometry: Mapping[str, float] | Sequence[float],
    element_names: Sequence[str],
) -> np.ndarray:
    """Return a stoichiometry vector in the requested element order."""
    normalized_elements = _normalize_names(element_names)
    if isinstance(stoichiometry, Mapping):
        upper_map = {
            str(key).strip().upper(): float(value)
            for key, value in stoichiometry.items()
        }
        return np.asarray(
            [upper_map.get(element.upper(), 0.0) for element in normalized_elements],
            dtype=float,
        )

    values = np.asarray(stoichiometry, dtype=float)
    if values.shape != (len(normalized_elements),):
        raise ValueError("stoichiometry vector length must match element_names.")
    return values


def _normalize_names(names: Sequence[str]) -> tuple[str, ...]:
    normalized = tuple(str(name).strip() for name in names)
    if not normalized:
        raise ValueError("At least one element is required.")
    if any(not name for name in normalized):
        raise ValueError("Element and component names cannot be empty.")
    return normalized


def _nnls(matrix: np.ndarray, rhs: np.ndarray, *, tolerance: float) -> np.ndarray:
    """Solve a small nonnegative least-squares problem by active set."""
    a = np.asarray(matrix, dtype=float)
    b = np.asarray(rhs, dtype=float)
    if a.ndim != 2:
        raise ValueError("matrix must be two-dimensional.")
    if b.shape != (a.shape[0],):
        raise ValueError("rhs length must match matrix rows.")

    _, n_columns = a.shape
    x = np.zeros(n_columns, dtype=float)
    passive = np.zeros(n_columns, dtype=bool)
    gradient = a.T @ (b - a @ x)
    max_iter = max(3 * n_columns, 1)
    iterations = 0

    while np.any((~passive) & (gradient > tolerance)):
        if iterations >= max_iter:
            break
        trial_scores = np.where(~passive, gradient, -np.inf)
        passive[int(np.argmax(trial_scores))] = True

        while True:
            z = np.zeros(n_columns, dtype=float)
            if np.any(passive):
                z[passive] = np.linalg.lstsq(a[:, passive], b, rcond=None)[0]

            if np.all(z[passive] > tolerance):
                x = z
                break

            blockers = passive & (z <= tolerance)
            if not np.any(blockers):
                x = np.maximum(z, 0.0)
                break

            denom = x[blockers] - z[blockers]
            valid = denom > tolerance
            if not np.any(valid):
                x[blockers] = 0.0
                passive[blockers] = False
                break

            alpha = float(np.min(x[blockers][valid] / denom[valid]))
            x = x + alpha * (z - x)
            small = passive & (x <= tolerance)
            x[small] = 0.0
            passive[small] = False
            if not np.any(passive):
                break

        gradient = a.T @ (b - a @ x)
        iterations += 1

    x[x < tolerance] = 0.0
    return x
