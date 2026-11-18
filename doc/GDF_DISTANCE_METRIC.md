# The GDF data distance

This note describes the **total** distance printed by
[`bin/gdf_distance`](../bin/gdf_distance), in exact real arithmetic. It is a
metric on the logical data described below, not on HDF5 files as byte strings.
For example, file layout, grid partitioning, unrelated attributes, boundary
types, and particle properties outside the supported list are not compared.
Two files differing only in such information represent the same logical data.

To generate a PDF, run `make gdf-distance-pdf` from the repository root.
This requires Pandoc and XeLaTeX with its standard LaTeX packages. The output
is `doc/GDF_DISTANCE_METRIC.pdf`, which is ignored by Git and rebuilt when
this source changes. If XeLaTeX reports a missing `xelatex.fmt`, initialize
the user format with `fmtutil-user --byfmt xelatex`.

## Domain and definition

Consider valid Piernik v2 GDF outputs accepted by the script, with finite,
real-valued time, grid fields and supported particle quantities. Their
grid levels are contiguous starting at zero, with refinement factor two when
there is more than one level. Each grid has the same canonical field names,
non-overlapping in-bounds cells on its level, unambiguous reported field-group
names, and the required data. Particle IDs are unique integers, including
integer values stored in floating-point datasets; each particle has all
quantities listed in `pkeys`.
Old and new CRESP field names are identified by `canonical_field_name`.

The structure of a data set is its dimensionality, geometry, domain edges and
dimensions, base-grid dimensions, canonical field names, covered-cell mask at
each level, and set of particle IDs. Different grid partitions of the same
covered cells are equivalent. If structures differ, the distance is **1**.

Within one structure, the logical values are the simulation time, all field
values on covered cells at all levels, and the supported particle values
indexed by ID. The script groups components deterministically into vectors, or
compares components individually with `--novec`; the chosen mode must be
fixed when considering a metric. For each value or vector it uses

\[
  \delta(x,y) =
  \begin{cases}
  0 & x=y=0,\\
  \dfrac{\lVert x-y\rVert_2}{\lVert x\rVert_2+\lVert y\rVert_2}
    & \text{otherwise}.
  \end{cases}
\]

For scalar time, this gives 1 for opposite signs. The total for matching
structures is \(D(x,y)=1-\prod_i(1-\delta_i(x_i,y_i))\), over the time,
field groups and particle groups. Other zero-valued report entries do not
change the result.

The optional level arguments compare projections containing only the selected
grid level, while retaining time and particles. The level number is part of
the projection's structure. The metric statement then applies to those
projections, not to the full input files. Three-file mode reports individual
component distances rather than this total.

## Proof

First, \(\delta\) is a metric on any finite-dimensional real Euclidean space.
It is non-negative and symmetric, is at most 1 by the ordinary triangle
inequality, and vanishes exactly when its arguments agree. To prove its
triangle inequality, put

\[
 a=\lVert x\rVert,\quad b=\lVert y\rVert,\quad c=\lVert z\rVert,\quad
 u=\lVert x-y\rVert,\quad v=\lVert y-z\rVert,\quad w=\lVert x-z\rVert.
\]

We need \(w/(a+c)\le u/(a+b)+v/(b+c)\), with zero denominators
interpreted as in the definition. Assume \(a\le c\), exchanging \(x,z\) if
necessary. The zero cases follow directly from the definition, so consider
positive \(a,b,c\). We have \(w\le u+v\), and also the Ptolemy inequality

\[
 bw\le cu+av.
\]

For the latter, apply the triangle inequality to the inverted points
\(p\mapsto p/\lVert p\rVert^2\), using
\(\lVert p/\lVert p\rVert^2-q/\lVert q\rVert^2\rVert
=\lVert p-q\rVert/(\lVert p\rVert\lVert q\rVert)\).

If \(b\le a\), both \((a+c)/(a+b)\) and \((a+c)/(b+c)\) are at
least 1, and the claim follows from \(w\le u+v\). If \(b>a\), take a
convex combination of the ordinary and Ptolemy bounds, with weights
\(a/(a+b)\) and \(b/(a+b)\):

\[
 w\le \frac{a+c}{a+b}u+\frac{2a}{a+b}v
   \le \frac{a+c}{a+b}u+\frac{a+c}{b+c}v.
\]

The last inequality follows from \((c-a)(b-a)\ge0\). Dividing by
\(a+c\) proves the triangle inequality. This includes scalar time as the
one-dimensional case.

Next, let \(d_i\) be any finite family of metrics bounded by 1. Their
probabilistic sum \(D=1-\prod_i(1-d_i)\) is a metric. Non-negativity,
symmetry, and the zero condition are immediate. For the triangle inequality,
write \(p_i=d_i(x,y)\) and \(q_i=d_i(y,z)\). Since
\(d_i(x,z)\le\min(1,p_i+q_i)\), it suffices to show
\(F(\min(1,p+q))\le F(p)+F(q)\), where
\(F(r)=1-\prod_i(1-r_i)\). To see this, take independent uniform random
variables \(U_i\in[0,1]\), and set
\(A_i=\{U_i<p_i\}\), \(B_i=\{U_i>1-q_i\}\).
Then \(P(A_i\cup B_i)=\min(1,p_i+q_i)\). Independence across indices and
the union bound give

\[
 F(\min(1,p+q))
 =P\!\left(\bigcup_i(A_i\cup B_i)\right)
 \le P\!\left(\bigcup_i A_i\right)+P\!\left(\bigcup_i B_i\right)
 =F(p)+F(q).
\]

Within a fixed structure, field and particle groups form a fixed collection
of coordinates, so the preceding result proves that the total is a metric on
their logical values. For different structures the distance is 1, while
every within-structure distance is at most 1. A triangle with endpoints in
different structures has at least one cross-structure edge of length 1; if
the endpoints have the same structure but the intermediate point does not,
both other edges have length 1. Thus the triangle inequality also holds
across structures. Distinct structures have nonzero distance, completing the
proof.

## Numerical and input limitations

The proof is for exact arithmetic, not a guarantee of exact metric axioms for
floating-point output or the printed rounded values. Grid and particle values
are loaded into float64 arrays. The implementation scales
norms to avoid common overflow and underflow; it reports an error if a
nonzero data difference underflows to zero. Rounding can still perturb
symmetry or the triangle inequality near equality. Inputs with NaN, infinity,
malformed or incomplete data are outside the stated domain; not every such
input is fully validated when a structure mismatch returns early.

Structure mismatches now return 1, rather than a graded score based on shared
fields or particles. This is an intentional change from the previous
comparison and may change regression results.
