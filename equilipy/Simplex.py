import numpy as np, math
from numba import njit

def simplex_count(nVertex:int,nSpacing:int):
    '''
    Count the number of simplex grid points.
    ----------------------------------------------------------------------------------------------------------------
    
    Input
    -----
    nVertex: The number of vertices (e.g. nVertex = 3 for a ternary system )
    
    nSpacing: The number of spacings in one dimension. Note that 10 spacings result in 11 points 
    (e.g. nSpacing = 10 for 0.1 grid spacing in the range of 0 to 1)

    Output
    ------
    The total number of points in the simplex grid.
    '''
    return math.comb(nSpacing+nVertex-1, nVertex-1)


def simplex_grid(nVertex:int,nSpacing:int):
    '''
    Calculate simplex grid evenly spaced between 0 and 1.
    ----------------------------------------------------------------------------------------------------------------

    Input
    -----
    nVertex: The number of vertices (e.g. nVertex = 3 for a ternary system )
    
    nSpacing: The number of spacings in one dimension. Note that 10 spacings result in 11 points 
    (e.g. nSpacing = 10 for 0.1 grid spacing in the range of 0 to 1)

    Output
    ------
    Grid matrix with size of (simplex_count,nVertex)
    
    '''

    L = simplex_count(nVertex, nSpacing)
    return np.asarray(_simplex_grid(L,nVertex,nSpacing))


@njit
def _simplex_grid(L,nVer,nSpacing):
    res = np.zeros((L, nVer))
    
    vertices=np.zeros(nVer)
    vertices[nVer-1] = nSpacing

    for i in range(nVer):
        res[0, i] = vertices[i]

    h = nVer

    for i in range(1, L):
        h -= 1

        val = vertices[h]
        vertices[h] = 0
        vertices[nVer-1] = val - 1
        vertices[h-1] += 1

        for j in range(nVer):
            res[i, j] = vertices[j]

        if val != 1:
            h = nVer
    return res/nSpacing


def simplex_grid_shift(SimplexGrid,CornerVertex:int,ShrinkFactor:float):
    '''
    Sift the simplex grid to the given CornerVertex with ShrinkFactor.
    ----------------------------------------------------------------------------------------------------------------

    Input
    -----
    SimplexGrid  : Grid matrix with size of (simplex_count,nVertex)
    CornerVertex : An integer indicating the vertex column in SimplexGrid
    ShrinkFactor : A float ranging from 0 to 1

    Output
    ------
    Grid matrix re-evaluated based on CornerVertex and ShrinkFactor
    '''
    
    Shift=np.zeros(SimplexGrid.shape)
    Shift[:,CornerVertex]=1-ShrinkFactor
    return SimplexGrid*ShrinkFactor+Shift