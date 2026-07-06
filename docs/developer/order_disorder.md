# Order/Disorder CEF And Magnetic Terms

This page documents Equilipy's current `SUBOM` implementation for ordered
compound-energy-formalism phases with a `DIS_PART` or `DISORDER_PART`
disordered partner. The intent is to separate the thermodynamic model from the
minimizer: fixed-site `G/H/S/Cp` must be correct before leveling, GEM, or
grid-seeding behavior is debugged.

## References

The implementation follows these model references:

- Hillert (2001), *The compound energy formalism*, for the CEF sublattice and
  endmember framework.
- Hillert and Jarl (1978), *A model for alloying effects in ferromagnetic
  metals*, for magnetic Gibbs, entropy, enthalpy, and heat-capacity terms.
- Zheng et al. (2017), *Thermodynamic assessment of the Al-C-Fe system*, for
  ordered/disordered bcc, fcc, and kappa partitioning.
- Noori and Hallstedt (2021), *Thermodynamic modelling of the Al-Co-Fe
  system*, for bcc/B2 order-disorder usage in an Al-Co-Fe database.

## Phase Scope

The Fortran runtime distinguishes:

```text
SUBL   CEF phase without magnetic terms
SUBLM  CEF phase with magnetic terms
SUBOM  ordered CEF phase with a disordered partner
```

`BCC_A2` and `A2_BCC` are aliases for the same disordered bcc physics.
Likewise, `FCC_A1` and `A1_FCC` are aliases for the same disordered fcc
physics. The writer may preserve helper aliases or canonicalize them, but the
runtime must evaluate one order/disorder relationship.

## Gibbs Energy

Let `y[s,k]` be the ordered phase site fraction of constituent `k` on ordered
sublattice `s`. Ordered sublattices map onto disordered sublattice groups `d`.
For each disordered group, collapse the ordered site fractions with ordered
site ratios `a[s]`:

```text
A[d]      = sum(s -> d, a[s])
x[d,k]   = (1 / A[d]) * sum(s -> d, a[s] * y[s,k])
y0[s,k]  = x[d,k] for every ordered sublattice s mapped to d
```

The ordered Gibbs energy is evaluated as:

```text
G_SUBOM(y) = G_ord_raw(y) + G_dis(x(y)) - G_ord_raw(y0(x(y)))
```

At the random state `y = y0`:

```text
G_SUBOM(y0) = G_dis(x)
```

so the ordering contribution vanishes. The bracket
`G_dis(x) - G_ord_raw(y0)` is not expected to vanish by itself.

## Analytical Chemical Potentials

Order/disorder chemical potentials are analytical. Numerical finite
differences are used only as development diagnostics.

For a scalar CEF property `Q(y)` represented by endmember partial values
`q_i`, the partial-molar transform for ordered endmember `i` is:

```text
Qbar_i = Q + sum_s a[s] * (R[s,c_i] - sum_c y[s,c] * R[s,c])
```

where `c_i` is the constituent on sublattice `s` for endmember `i`, and
`R[s,c]` is the relative conditional derivative:

```text
R[s,c] = <q_i | endmember i has constituent c on sublattice s> - Q
```

For the order/disorder correction this transform is applied to:

```text
Delta Q = Q_dis(x) - Q_ord_raw(y0)
```

and mapped back through `s -> d`:

```text
Delta Qbar_i += (a[s] / A[d]) *
    (Rcorr[d,c_i] - sum_c y[s,c] * Rcorr[d,c])
```

The runtime adds these analytical corrections to:

```text
dPartialExcessGibbs
dPartialEnthalpyXS
dPartialEntropyXS
dPartialHeatCapacityXS
```

## Magnetic Variable Partition

Hillert-Jarl magnetic Gibbs energy is nonlinear in Curie/Neel temperature
`TC` and magnetic moment `B` (`BM`, `BMAG`, or `BMAGN` in TDB spelling). For
`SUBOM`, Equilipy partitions these variables first:

```text
TC_part(y) = TC_ord(y) + TC_dis(x) - TC_ord(y0)
B_part(y)  = B_ord(y)  + B_dis(x)  - B_ord(y0)
```

Then the Hillert-Jarl function is evaluated once from `TC_part` and `B_part`.
The retired value-level expression

```text
Gmag_ord(y) + Gmag_dis(x) - Gmag_ord(y0)
```

is not equivalent because `Gmag(TC, B)` is nonlinear.

The magnetic partials mirror the existing `.dat`-compatible magnetic routine.
For endmember `i`, with the same auxiliary variables used in
`CompGibbsMagneticSoln.f90`:

```text
P1 = (B_i - B) / (1 + B)
P2 = ln(1 + B)
P3 = -(TC - TC_i)
```

the analytical magnetic partials are:

```text
G_i^mag  = F*P1 + P2*(P3*Qp + F)

H_i^mag  = -(TC*P1*Fp + P2*(TC*Fp + P3*(Qp + Qdp)))

S_i^mag  = -P1*(F + TC*Fp)
           -P2*(F + TC*Fp + P3*(2*Qp + Qdp))

Cp_i^mag = -P1*TC*(2*Fp + Fdp)
           -P2*(TC*(2*Fp + Fdp)
                + P3*(2*Qp + 4*Qdp + Qtp))
```

The final magnetic correction applied to each ordered endmember is:

```text
Delta X_i^mag =
    (X_i^mag,part - X^mag,part)
  - (X_i^mag,ord_raw - X^mag,ord_raw)
```

for `X` in `{G, H, S, Cp}`. The wrapper already carries the raw ordered
magnetic partial; this correction replaces it with the partitioned magnetic
contribution without double counting.

## Code Map

The main implementation lives in:

```text
equilipy/fsrc/4_GibbsModels/CompExcessGibbsEnergySUBL.f90
  ApplyOrderDisorderCorrectionSUBL
  EvaluateMagneticScalarSUBL
  EvaluateMagneticPartialSUBL
  CompMagneticVariablePartialsSUBL
```

Supporting files:

```text
equilipy/fsrc/4_GibbsModels/CompGibbsMagneticSoln.f90
equilipy/fsrc/4_GibbsModels/CompMagneticTemperatureMoment.f90
equilipy/database_ir/runtime.py
equilipy/database_ir/tdb_writer.py
```

## Validation

The focused regression suite validates the same shared runtime from both DAT
and TDB inputs:

```bash
source /Users/69e/miniforge3/etc/profile.d/conda.sh
conda activate equilipy
zsh bdist_wheel.sh
pytest tests/tdb_bcc_a2_magnetic_test.py \
       tests/tdb_runtime_f_option_test.py \
       tests/tdb_order_disorder_alias_test.py \
       tests/database_ir_test01.py -q
```

The BCC_A2 rows validate the disordered magnetic phase against the Al-Fe DAT
reference. The BCC_B2 rows validate DAT, direct TDB, and TDB roundtrip against
TC-Python 2025b variable-partition reference values for `G`, `S`, and `Cp`.

Before changing the minimizer, any new order/disorder work should use the
fixed-site workflow: obtain TC-Python site fractions, plug the same fractions
into Equilipy fixed evaluation, compare `G`, then `H/S`, then `Cp`.

The fixed-site and forced-phase order/disorder thermodynamics are considered
the reference layer. If ordered phases such as `B2_BCC` or `FCC_4SL` cause slow
or unstable equilibrium calculations, treat that as a minimizer initialization,
phase-selection, or subminimization problem first. Do not change the `SUBOM`
thermodynamic equations unless a fixed-site `G/H/S/Cp` comparison fails.
