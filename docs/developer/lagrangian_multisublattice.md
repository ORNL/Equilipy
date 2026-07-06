# Lagrangian GEM For Multisublattice CEF Phases

This page records the target formulation for fixed-active-set Lagrangian Gibbs
energy minimization when one or more active solution phases use a multisublattice
CEF model (`SUBL`, `SUBLM`, or `SUBOM`).  The immediate motivation is that the
legacy GEM Newton system treats every CEF endmember as an independent species.
For phases such as `BCC_B2`, `FCC_4SL`, `SIGMA`, and sublattice liquids, this
creates a poor Newton direction because the true independent coordinates are
site fractions, not pseudo-endmember mole fractions.

## Variable Choice

For an active CEF phase \(p\), write the phase amount as \(N_p\) and the
independent site-fraction variables as \(\eta_{p,a}\).  Each
\(\eta_{p,a}\) is an exchange variable on one sublattice:

\[
\eta_{p,a} = y_{s,c}, \qquad
y_{s,r} = 1 - \sum_{c \ne r} y_{s,c},
\]

where \(r\) is the reference constituent on sublattice \(s\).  The reference is
chosen from the largest current site fraction so the reduced coordinate remains
well conditioned.

The phase Gibbs energy and phase composition are functions of the site
fractions:

\[
g_p = g_p(y_p),
\qquad
a_{p,e} = a_{p,e}(y_p).
\]

For a stoichiometric compound \(c\), the amount is \(N_c\), the Gibbs energy is
\(G_c\), and the composition vector is fixed as \(a_{c,e}\).

## Fixed-Assemblage Lagrangian

Using dimensionless internal units, define

\[
\mathcal{L}
= \sum_p N_p g_p(y_p)
  + \sum_c N_c G_c
  - \sum_e \lambda_e
    \left[
      \sum_p N_p a_{p,e}(y_p)
      + \sum_c N_c a_{c,e}
      - b_e
    \right].
\]

The fixed-assemblage residuals are:

\[
R^N_p = g_p(y_p) - \sum_e \lambda_e a_{p,e}(y_p),
\]

\[
R^\eta_{p,a}
= N_p
  \left[
    \frac{\partial g_p}{\partial \eta_{p,a}}
    - \sum_e \lambda_e
      \frac{\partial a_{p,e}}{\partial \eta_{p,a}}
  \right],
\]

\[
R^\lambda_e
= b_e
 - \sum_p N_p a_{p,e}(y_p)
 - \sum_c N_c a_{c,e}.
\]

For a compound:

\[
R^N_c = G_c - \sum_e \lambda_e a_{c,e}.
\]

## Newton/KKT Matrix

The mixed-coordinate Newton step solves for:

\[
\Delta z =
\left[
  \Delta N,\ \Delta \eta,\ \Delta \lambda
\right]^T.
\]

The symmetric KKT matrix is assembled from the second variation of
\(\mathcal{L}\):

\[
\begin{bmatrix}
0              & R^\eta     & -A_N^T \\
(R^\eta)^T     & N_p H_p    & -N_p A_\eta^T \\
-A_N           & -N_p A_\eta & 0
\end{bmatrix}
\begin{bmatrix}
\Delta N \\ \Delta \eta \\ \Delta \lambda
\end{bmatrix}
=
-
\begin{bmatrix}
R^N \\ R^\eta \\ R^\lambda
\end{bmatrix}.
\]

Here

\[
A_{N,p,e} = a_{p,e},
\qquad
A_{\eta,p,a,e} =
\frac{\partial a_{p,e}}{\partial \eta_{p,a}},
\]

and

\[
H_{p,a,b}
=
\frac{\partial^2 g_p}{\partial \eta_{p,a}\partial \eta_{p,b}}
- \sum_e \lambda_e
  \frac{\partial^2 a_{p,e}}{\partial \eta_{p,a}\partial \eta_{p,b}}.
\]

For ordinary CEF constituents, \(a_{p,e}(y_p)\) is linear in site fractions, so
the second composition derivative is zero and \(H_p\) is the analytical CEF
site Hessian from `CompHessianSUBL`.

## CEF Plus Other Phases

This formulation naturally couples CEF phases, ordinary solution phases, and
stoichiometric compounds through the same elemental-potential vector
\(\lambda\) and the same mass-balance rows.  Each active phase contributes one
phase-amount row.  CEF phases additionally contribute site-fraction exchange
rows.  Non-CEF solution models can keep the legacy local composition variables
until equivalent analytical gradients and Hessians exist for those models.

## Positivity

The Newton step is not allowed to make \(N_p\) or any site fraction negative.
The CEF line search therefore applies a feasibility cap:

\[
0 < N_p + \alpha \Delta N_p,
\qquad
0 < y_{s,c} + \alpha \Delta y_{s,c},
\]

including the reference constituent update

\[
\Delta y_{s,r} = -\sum_{c \ne r} \Delta y_{s,c}.
\]

This replaces the pseudo-endmember lambda correction for active CEF internal
coordinates.  A phase amount that wants to become negative should still be
reported as an active-set problem, not hidden by arbitrary species floors.

## Raw Phase-Amount Boundary Check

For CEF fixed-assemblage GEM, phase-amount positivity and site-fraction
positivity are different mathematical events.

A site-fraction boundary event is internal to an active solution phase:

\[
0 < y_{s,c} + \alpha \Delta y_{s,c}.
\]

The CEF line search owns this case because the phase remains active and only
its internal constitution must stay feasible.

A phase-amount boundary event changes the active assemblage:

\[
N_p^\mathrm{raw} = N_p + \Delta N_p < 0.
\]

This must be handled before the line search.  If the line search simply clips
or damps the step, the phase never reaches the active-set boundary and the
minimizer can report artificial stagnation even though the Newton direction has
already identified the phase that should leave the assemblage.

The production CEF path therefore audits the raw Newton phase direction
immediately after `GEMNewtonCEF` and before `GEMLineSearchCEF`:

```text
RunLagrangianGEM
  -> GEMNewtonCEF
  -> AuditCEFRawPhaseDirection
       if any N_p + Delta N_p < -dTolerance(8):
           mark PHASE_CHANGE_REASON_NEGATIVE_PHASE_AMOUNT
           record the most negative raw phase target
  -> RemoveRawNegativeCEFPhase
       compact iAssemblage, dMolesPhase, and dMolesSpecies
  -> CompChemicalPotential
  -> CompFunctionNorm
  -> retry RunLagrangianGEM iteration on the reduced active set
```

This is deliberately **not** a PEA replacement.  PEA remains responsible for
phase discovery, phase addition, and broad active-set replacement.  The raw
CEF phase-removal rule only removes a phase that is already active and whose
full Newton phase-amount direction crosses the zero boundary.  This follows
the practical split:

- PEA is good at finding a new candidate active set.
- Lagrangian GEM is more reliable at recognizing when one currently active
  phase amount should become inactive.
- Line-search positivity control should protect CEF site fractions, not mask a
  phase-removal event.

Postprocess heat-capacity solves keep the fixed assemblage and therefore skip
this removal audit.

## Implementation Rule

For active CEF phases, Lagrangian GEM should use the same thermodynamic scalar,
gradient, and Hessian family:

- `CompGradientSUBL` for \(g_p\) and \(\partial g_p/\partial y\)
- `CompHessianSUBL` for \(\partial^2 g_p/\partial y^2\)
- CEF site-fraction residuals in `CompFunctionNorm`
- CEF site-fraction line search for positivity-preserving updates

This keeps Leveling, PEA/subminimization, and final fixed-assemblage
Lagrangian minimization on the same CEF thermodynamic surface.

## Current Status

The mixed-coordinate CEF Lagrangian path is enabled from `InitGEMSolver` and
uses `GEMNewtonCEF`, `GEMLineSearchCEF`, and the raw phase-amount boundary
audit described above.

Validation on Jun. 25, 2026:

- The two extracted Example03 Al-Cu-Mg-Si NaN rows converged after adding the
  raw CEF phase-removal audit.
- A full Example03 sweep over `AlCuMg`, `AlCuSi`, `AlMgSi`, `CuMgSi`, and
  `AlCuMgSi` produced no nonfinite `G/H/S/Cp` values and no negative `Cp`
  values.
- The same sweep still found high-Cp rows above 1000 J/K in `AlMgSi` and
  `CuMgSi`; those are separate Cp/phase-boundary validation cases and should
  not be hidden by this active-set fix.
