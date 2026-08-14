# Method

Full description of what the tool computes, and why each choice was made.
For a summary, see the [README](../README.md); for every command-line flag, see [OPTIONS.md](OPTIONS.md).

## Collective variables

Collective variables (CVs) are coordinates that describe the macroscopic state of a system. They
reduce a trajectory of millions of atomic coordinates to the few degrees of freedom that matter:
principal components, a distance between two atoms, an angle, a dihedral, or any descriptor you can
compute per frame.

The tool takes two of them. The input is any two-column text file — frame index, CV value — so the
CV can come from PCA, from a geometric measurement, or from anything else.

## Step 1 — Density

The two CV columns give one point per frame in a 2D space. A Gaussian kernel density estimate
turns those points into a smooth probability density. For $n$ points $\{x_i\}$:

$$\hat{f}(x) = \frac{1}{n h} \sum_{i=1}^{n} K\!\left( \frac{x - x_i}{h} \right)$$

$K$ is the Gaussian kernel and $h$ the bandwidth, which controls smoothness. Without a bandwidth
given, `scipy.stats.gaussian_kde` uses Scott's rule. The estimate is evaluated on a regular grid
(`--grid`, default 100×100).

A **single** kernel is fitted over the whole trajectory. Parallelism (`--n_jobs`) splits only the
*evaluation* of that one kernel, so the number of cores never changes the result.

> Version 1.x fitted a separate kernel inside each parallel chunk. Since chunks were contiguous in
> frame order, each covered a different region of CV space and was normalised to itself, and the
> answer depended on `cpu_count()`. See [CHANGELOG.md](../CHANGELOG.md) for the measured impact.

## Step 2 — Free energy

Boltzmann inversion converts density into free energy:

$$G(x,y) = -k_B T \ln \hat{f}(x,y), \qquad G \leftarrow G - \min G$$

The shift puts the global minimum at zero, so every energy the tool reports is relative to the
deepest point of the landscape. $k_B$ is set by `--kb` (default 8.314e-3 kJ/mol·K) and $T$ by
`--temperature` (default 300 K).

## Step 3 — Sampling mask

A KDE returns a number everywhere, including regions that contain no frame at all. Those numbers
are extrapolation, not measurement — on the bundled trajectory they reach ~1700 kJ/mol.

A grid cell counts as sampled only when the estimate expects at least `--mask_min_count` frames
inside it:

$$\hat{f} \cdot A_{\text{cell}} \cdot N \ \geq\ \texttt{mask\_min\_count}$$

Everything else is excluded from the analysis and from the colour scale. On the bundled trajectory
this keeps 1173 of 10000 grid cells.

Because the criterion is a threshold on $\hat{f}$, and $G$ is a monotonic function of $\hat{f}$, the
mask boundary is exactly a contour of constant free energy.

This is why `--max_energy` defaults to `none`: with the extrapolated region gone, there is no
runaway range left to clip. Version 1.x hard-coded a 25 kJ/mol ceiling, which flattened every real
barrier above it and turned most of the figure into a saturated plateau.

## Step 4 — Steepest descent

Each cell looks at its 8 neighbours — the full 3×3 block minus the centre — and points at the one
with the **lowest** energy. If no neighbour is lower, the cell points at itself: it is a local
minimum.

Following those pointers to their fixed point tells which minimum each cell drains into. That
minimum defines its **basin of attraction**. Since a step only accepts a strictly lower neighbour,
no cycle can form and the iteration always terminates; the implementation reaches the fixed point by
pointer doubling, in $O(\log)$ passes.

### Why eight neighbours

Eight is the complete Moore neighbourhood in two dimensions — $3^d - 1$ for $d = 2$.

With only the four axis directions, descent cannot move diagonally, so a diagonal valley makes it
zigzag and stall: a cell can be lower than its four axis neighbours while a diagonal one is lower
still, which registers as a minimum that does not exist in the underlying continuous surface. On the
bundled trajectory:

| neighbourhood | cells | raw local minima |
|---|---|---|
| 4 (von Neumann) | 4 | 25 |
| **8 (Moore, r=1)** | **8** | **12** |
| 24 (r=2) | 24 | 9 |
| 48 (r=3) | 48 | 9 |
| 120 (r=5) | 120 | 7 |

Thirteen of the 25 minima found with four neighbours are artefacts of the grid.

A **larger radius is not a generalisation** — it is a different, worse algorithm. With radius > 1 a
cell no longer follows the gradient; it jumps to the best cell within the radius. If the radius spans
a ridge, the cell jumps straight over the barrier and the two basins merge. That silently does what
`--basin_min_depth` does explicitly, except without measuring the barrier and without reporting it.

|  | mechanism | reportable |
|---|---|---|
| large radius | merges by spatial proximity | no |
| `--basin_min_depth` | merges by saddle height, in kJ/mol | yes, with the full spectrum printed |

The knobs that legitimately control smoothness are `--kde_bandwidth`, `--grid` and
`--basin_min_depth`, because all three change $G$ or the reported depths, and are therefore visible
in the output.

The genuine limitation of eight neighbours is that flow directions are quantised to 45° — the same
D8 problem known from terrain hydrology, which runs this exact algorithm on elevation maps. The
remedies there are D∞ or multiple-flow-direction, not a bigger radius.

## Step 5 — Saddles

Two neighbouring basins meet along a ridge. The **saddle** between them is the lowest point of that
ridge — the cheapest crossing:

$$\text{saddle}(A,B) = \min_{\substack{x \in A,\ y \in B \\ x,\,y \text{ adjacent}}} \max\big(G(x),\,G(y)\big)$$

## Step 6 — Depth

$$\text{depth}(A) = \min_B \text{saddle}(A,B)\ -\ \min_{x \in A} G(x)$$

How far a basin must climb to escape. This is its topological persistence, and it is what decides
whether a dip is a state or a wrinkle. Unlike an energy cut-off, it does not depend on the grid
resolution.

A basin with no neighbour has no saddle, so its depth is undefined — reported as `-`. It sits on its
own island of sampled space.

## Step 7 — Merging

Basins shallower than `--basin_min_depth` are absorbed, one at a time, by the neighbour across their
lowest saddle, recomputing after each merge.

The cut defaults to automatic: sort all depths and place it at the **largest step** in that list.
Real structure and estimator ripple usually sit on opposite sides of a wide gap. The automatic cut
never exceeds $k_B T$, because a barrier taller than the thermal scale is real by definition — below
it, thermal agitation crosses freely and the two sides are not distinct states.

What it chose, and what it merged, is printed on every run:

```
Watershed found 12 local minima; --basin_min_depth 1.01 kJ/mol (auto) merged 7 of them.
Barrier depths (deepest first): 2.10, 1.81, 1.59, 1.43, 1.41, 0.61, 0.47, 0.27, 0.17, ...
Largest step is 1.41 -> 0.61 kJ/mol: --basin_min_depth between them keeps 6 basins.
```

On this trajectory a fixed $k_B T = 2.49$ would have merged all twelve minima into a single basin
holding 98.6% of the frames. The automatic cut at 1.01 keeps six, including three sub-basins holding
5201, 9263 and 36447 frames.

## Step 8 — Frames and representatives

Each frame falls in a grid cell and inherits its basin id. The **representative** of a basin is the
frame of lowest free energy inside it — an actual conformation from the trajectory, ready to
extract, not a grid node.

Basins are numbered by increasing energy, so basin 1 is always the global minimum.

## A worked example, in one dimension

The same procedure on a small 1D profile, where every step can be checked by hand:

```
index:    0     1     2     3     4     5     6     7     8     9    10
G:       6.0   3.0   0.0   2.0   4.5   3.8   4.6   2.2   0.9   3.0   6.5
                      ^                 ^                 ^
                    min A             min B             min C
descent:  →     →     •     ←     ←     •     →     →     •     ←     ←
```

Cell 6 has neighbours 5 (3.8) and 7 (2.2), so it descends to 7, not to 5. That is why basin B ends
up holding a single cell.

```
basins:   A = {0,1,2,3,4}      B = {5}      C = {6,7,8,9,10}

saddles:  A|B at (4,5):  max(4.5, 3.8) = 4.5
          B|C at (5,6):  max(3.8, 4.6) = 4.6

depths:   A: 4.5 - 0.0 = 4.5
          C: 4.6 - 0.9 = 3.7
          B: 4.5 - 3.8 = 0.7      (its lowest saddle is the one shared with A)

spectrum: 4.5    3.7  |  0.7
                      ↑ largest step (3.0) → cut at 2.2
```

B is shallower than the cut, so it is absorbed by A, across its lowest saddle. Two states remain:
`A = {0..5}` and `C = {6..10}`. B was a 0.7 kJ/mol kink on the way up, not a state.

## Determinism

Nothing in the pipeline is random or iterative. There is no `k` to choose, no seed, no
initialisation. The same input gives byte-identical output on any machine and for any `--n_jobs`,
which the test suite enforces.

## Why not DBSCAN, k-means or a Gaussian mixture

| method | uses the energy surface | gives barrier heights | assigns every frame | parameter |
|---|---|---|---|---|
| **watershed + persistence** | yes | yes (saddle) | yes | kJ/mol, physical |
| connected components | yes | no | no | kJ/mol |
| DBSCAN | no | no | no (noise label) | `eps`, `minPts` |
| HDBSCAN | no | hierarchy only | no | `min_cluster_size` |
| k-means / GMM | no | no | yes | `k`, not deterministic |
| MSM + PCCA+ | kinetics | yes (timescales) | yes | lag time, states |

DBSCAN groups by local density — which is exactly the quantity the free energy is the logarithm of.
Running it on the points is therefore a density cut-off in disguise, and it throws away what Step 2
already produced: no barrier heights, no ΔG, parameters with no physical unit (`eps` would mix
degrees with ångströms), and the low-density points — the transition states — discarded as noise.

HDBSCAN prunes its hierarchy by persistence, which is the same idea as `--basin_min_depth`, but
applied to density rather than to energy.

Density-based and hierarchical methods do win in high dimension, where no grid can be built:
watershed needs a discretised surface, which is practical only for one to three CVs. To cluster
conformations by RMSD in a hundred dimensions, they are the right tools.

The real upgrade over what is implemented here is not DBSCAN but a **Markov state model with
PCCA+**, which groups by kinetics — what interconverts quickly is one state — and yields transition
timescales alongside ΔG. That is a substantially larger piece of work and would add a dependency.

## Known limitations

- **No error bars.** MD frames are strongly autocorrelated, so the effective sample size is far
  below the frame count. Block averaging or a block bootstrap would be needed to put uncertainties
  on ΔG and on barrier heights.
- **No Jacobian correction.** For curvilinear CVs (radial distance, angle, RMSD) the phase-space
  volume element biases the profile.
- **Cost.** `gaussian_kde` evaluation is $O(N^2)$; around $2 \times 10^5$ frames it becomes
  impractical, and a histogram-based or FFT-based estimator would be needed.
- **Two CVs.** The grid, the descent and the saddle search all assume $d = 2$.
