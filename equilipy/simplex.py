"""Simplex grid generation helpers."""

from __future__ import annotations

import math

import numpy as np

from equilipy.equilifort import simplex_grid_f


def simplex_count(n_vertex: int, n_spacing: int):
    """
    Count the number of simplex grid points.

    Input
    -----
    n_vertex: The number of vertices (e.g. n_vertex = 3 for a ternary system )

    n_spacing: The number of spacings in one dimension. Note that 10 spacings
    result in 11 points
    (e.g. n_spacing = 10 for 0.1 grid spacing in the range of 0 to 1)

    Output
    ------
    The total number of points in the simplex grid.
    """
    return math.comb(n_spacing + n_vertex - 1, n_vertex - 1)


def simplex_grid(n_vertex: int, n_spacing: int):
    """
    Calculate simplex grid evenly spaced between 0 and 1.

    Input
    -----
    n_vertex: The number of vertices (e.g. n_vertex = 3 for a ternary system )

    n_spacing: The number of spacings in one dimension. Note that 10 spacings
    result in 11 points
    (e.g. n_spacing = 10 for 0.1 grid spacing in the range of 0 to 1)

    Output
    ------
    Grid matrix with size of (simplex_count,n_vertex)

    """
    n_points = simplex_count(n_vertex, n_spacing)
    return np.asarray(_simplex_grid(n_points, n_vertex, n_spacing))


def _simplex_grid(n_points, n_vertex, n_spacing):
    res = simplex_grid_f(n_points, n_vertex, n_spacing)
    return res


def simplex_grid_shift(simplex_grid, corner_vertex: int, shrink_factor: float):
    """
    Sift the simplex grid to the given corner_vertex with shrink_factor.

    Input
    -----
    simplex_grid  : Grid matrix with size of (simplex_count,n_vertex)
    corner_vertex : An integer indicating the vertex column in simplex_grid
    shrink_factor : A float ranging from 0 to 1

    Output
    ------
    Grid matrix re-evaluated based on corner_vertex and shrink_factor
    """
    Shift = np.zeros(simplex_grid.shape)
    Shift[:, corner_vertex] = 1 - shrink_factor
    return simplex_grid * shrink_factor + Shift
