# Changelog

All notable changes to this project are documented here.
This project adheres to [Semantic Versioning](https://semver.org/).

## [2.0.1] — 2026-08-15

No change to the analysis or its results. Same code, same numbers as 2.0.0.

### Added

- Releases publish to PyPI through Trusted Publishing. PyPI authenticates the workflow directly
  over OIDC, so there is no API token to create, store, rotate or leak. Publishing a GitHub release
  now builds, tests and uploads on its own; nothing has to be run by hand.
- The release workflow refuses to publish when the tag and the version in `pyproject.toml`
  disagree, and runs the suite and `twine check` before uploading — a PyPI version can never be
  replaced once it is up.

### Changed

- README trimmed. The project page on PyPI carries whatever README shipped with the version that
  built it, and cannot be edited afterwards, so bringing it in line with the repository takes a
  release of its own — which is what this one is.

## [2.0.0] — 2026-08-14

**This release changes numerical results.** Any analysis produced with version 1.x should be
re-run. See "Fixed" below for the magnitude of the difference.

### Fixed

#### Parallel KDE was fitted per chunk (critical, affects all published numbers)

`calculate_and_save_free_energy` split the trajectory into `multiprocessing.cpu_count()`
contiguous blocks and fitted a **separate** `gaussian_kde` inside each block, then concatenated
the resulting densities as if they were comparable.

Because `np.array_split` cuts along frame order, each block covered a different region of CV
space, and each block's density was normalised to itself. A conformation that was rare over the
whole trajectory but common inside one block received a spuriously low free energy. The bandwidth
(Scott's rule) was also recomputed per block, and the number of blocks was the machine's CPU
count — so **the same input produced different results on different machines**.

Measured on the bundled example trajectory (51,643 frames), against a single global kernel:

| chunks (= CPU count) | max &#124;ΔG − ΔG_ref&#124; | frames ≤ 1 kJ/mol | Jaccard vs. correct set |
|---|---|---|---|
| reference (1 kernel) | 0.000 kJ/mol | 12,003 | 1.000 |
| 2 | 11.4 kJ/mol | 12,197 | 0.506 |
| 4 | 12.7 kJ/mol | 6,702 | 0.081 |
| 8 | 17.4 kJ/mol | 1,981 | 0.146 |
| 16 | 17.9 kJ/mol | 1,719 | 0.123 |

On an 8-core machine the tool selected 1,981 frames instead of 12,003, and roughly 85% of the
selected frames were wrong.

A single `gaussian_kde` is now fitted over the whole dataset (`get_kernel`), and parallelism
applies only to *evaluating* that fitted kernel (`evaluate_density`). Results are now bit-identical
for any value of `--n_jobs`, and a regression test enforces this.

#### The grid surface and the per-frame energies used different kernels

The contour/3D plots and the `.tsv` export each fitted their own KDE, so the plotted landscape and
the exported per-frame energies described different surfaces in the same run. Both paths now share
one kernel via `get_kernel`.

#### `calculate_free_energy` returned a stale cache

The method returned `self.cached_results` whenever it existed, ignoring its `data` argument, the
bandwidth, the grid size and the axis limits. Reusing an instance for a second dataset silently
returned the first result. The cache is now keyed on all of those inputs.

#### No validation of the input files

CV files with different lengths, mismatched frame indices, fewer than two columns, or NaN/inf
values were accepted and failed later with an obscure error — or, worse, produced a landscape that
paired CV1 and CV2 across *different* conformations. `load_data` now validates all of these and
raises an actionable message. When the frame counter restarts at different rows in the two files
(concatenated trajectories of unequal length) the error points at the offending rows.

#### `log(0)` warnings

Densities are clipped to the smallest positive float before the logarithm, removing
`RuntimeWarning: divide by zero` and the resulting `-inf` cells.

### Changed

- **The 25 kJ/mol ceiling is gone.** v1.x hard-coded `np.clip(G, 0, 25)`, silently flattening
  every barrier above 25 kJ/mol. It is now `--max_energy`, defaulting to `none` — masking the
  unsampled region removed the reason the cap existed. When set, it applies **only to the plotted
  surface**; the per-frame energies in the `.tsv` are never truncated.
- **Grid resolution is configurable** via `--grid` (default 100, previously hard-coded).
- **`--n_jobs`** controls the number of parallel workers. It no longer affects results.
- `calculate_and_save_free_energy` accepts a `filename` argument instead of hard-coding
  `discrete_values_energy_frames.tsv`.
- Frame indices are read once in `load_data` and exposed as `self.frames`, instead of re-reading
  the CV1 file from disk.
- Unhandled exceptions now print a traceback instead of only the message.
- `-h` / `--help` works as the first argument.
- Options given without their value now report a clear error instead of raising `IndexError`.
- `requirements.txt` uses lower bounds and matches `pyproject.toml` (it previously pinned exact
  old versions that contradicted the package metadata). See "Supported Python" below.

### Added

#### Basin clustering replaces point-by-point marking

`--energy T` previously coloured **individual grid nodes** whose energy fell in each
`--discretize` band. Those nodes are cells of the KDE grid, not sampled conformations, and two
nodes in the same energy band could sit in regions with no path between them — so the marked
points were not representative structures and did not correspond to the frames listed in the
`.tsv`.

The surface is now segmented into **basins of attraction** by watershed. Every grid cell follows
the steepest downhill neighbour until it reaches a local minimum; the minimum it lands on defines
its basin. This is the physical definition of a metastable state, it partitions the whole sampled
region, it needs no threshold, and it leaves no frame unassigned. Deterministic — no `k` to pick,
no random seed — and built on numpy and scipy, both already dependencies.

Adjacent basins separated by a barrier shallower than `--basin_min_depth` are merged. The
criterion is topological persistence: a basin's depth is the height of the lowest saddle to a
neighbour, above its own minimum. Unlike a threshold, persistence does not depend on the grid
resolution.

The cut defaults to automatic, placed at the largest step in the sorted depths — real structure
and estimator ripple usually fall on opposite sides of a wide gap — and never above `kB·T`, since
a barrier taller than the thermal scale is real by definition. Every run prints how many local
minima were found, the full depth spectrum, the cut that was applied and how many basins it
merged, so the choice is never invisible. On the bundled trajectory the spectrum is
`2.10 1.81 1.59 1.43 1.41 | 0.61 0.47 0.27 …`: a fixed `kB·T = 2.49` would have merged all twelve
minima into a single basin holding 98.6% of the frames, while the automatic cut at 1.01 keeps six.

For each basin the tool reports the position and depth of its minimum, how many frames it holds,
its population fraction, and a **representative frame**: the actual sampled conformation with the
lowest free energy inside that basin, ready to extract from the trajectory. Basins are numbered by
increasing energy, so basin 1 is always the global minimum.

With the watershed method `--energy` no longer changes the partition — it only filters which
basins are reported (those whose minimum lies below it). The old behaviour is still available as
`--basin_method connected`.

New outputs:

- `basins.tsv` — one row per basin: `cluster`, `cv1_min`, `cv2_min`, `G_min_kJ_mol`,
  `depth_kJ_mol`, `n_frames`, `population_fraction`, `rep_frame`, `rep_cv1`, `rep_cv2`,
  `rep_energy_kJ_mol`, `area_cells`.
- A summary table printed to the terminal.
- A `cluster` column in `discrete_values_energy_frames.tsv`, so every frame carries the id of its
  basin (0 = outside the sampled region).
- The 2D and 3D landscapes draw the watershed ridges and mark each minimum with a numbered
  callout; the legend gives ΔG, depth, population and representative frame.

New options: `--basin_method` (`watershed`, default, or `connected`), `--basin_min_depth`,
`--basin_connectivity` (`connected` only) and `--basin_min_frames`.

`--discretize` still works and keeps the old per-node rendering, now documented as legacy.

#### Unsampled regions are no longer painted as measured energy

The KDE was evaluated over the full bounding box and every cell was drawn, including regions
containing not a single frame. Extrapolated values there reached ~1700 kJ/mol, which the old
hard-coded 25 kJ/mol ceiling flattened into a saturated plateau covering most of the canvas — the
plateau, not the landscape, dominated the figure.

A grid cell now counts as sampled only when the KDE expects at least `--mask_min_count` frames
inside it (default 1). The rest is left blank. With that, the ceiling is no longer needed and
`--max_energy` defaults to `none`; the colour scale spans the range that was actually measured.

#### Figure quality

- Axes zoom to the sampled region instead of the full data rectangle, so the basins are no longer
  slivers in a corner.
- The area beyond the sampled region takes the top colour of the palette instead of staying blank,
  so the figure reads as one continuous field rather than a landscape cut out of a white page. Its
  edge is marked by a dashed white contour, which is where measurement actually stops.
- Perceptually uniform, colourblind-safe colormap by default: `viridis_r`, so minima are bright
  and the unsampled background recedes instead of glaring. `--cmap` takes any matplotlib name, and
  `classic` restores the 1.x blue-to-red map.
- Contour levels span only the useful range, in round steps, with thin light lines instead of a
  black hairline mesh.
- Colorbar sits against the axes and the legend moves below the plot; neither steals plot area.
- Basin labels are drawn as callouts with leader lines, and their placement is solved rather than
  fixed: each label is tried in several directions and distances from its marker and keeps the
  first position clear of every label already placed and of every other marker. Fixed offsets
  cycled by index collided as soon as two basins landed near each other on screen. The search
  works in screen pixels, since two basins can be far apart in CV units and touching in the
  figure.
- Figures are saved at 200 dpi (`--dpi`).
- The 3D surface is cut a few cells beyond the sampled region rather than exactly at its edge.
  Cutting at the edge sliced the walls mid-rise: the sampled region is a diagonal ribbon and the
  grid is orthogonal, so every row ended at a different height and the silhouette came out as a
  sawtooth of fins.
- The sampled region no longer carries a thick white border. The basin outline, the mask boundary
  and the outermost watershed ridge all fell on nearly the same curve and stacked into a halo, which
  made the landscape read as a sticker pasted onto the background. Ridges are now drawn only between
  two real basins, and the sampling boundary is a faint dashed line.
- The GIF no longer flashes: the colorbar is drawn on every frame, not only the first and last.
- `discrete_values_energy_frames.tsv` is always written. It previously required `--energy`, even
  though the per-frame table is the main quantitative product and the basin column exists
  regardless of any display cut-off.
- The 1D free energy profiles no longer floor empty histogram bins at a density of `1e-10`, which
  drew them as ~57 kJ/mol spikes indistinguishable from measured barriers. Empty bins are now a
  break in the curve — the same principle as the sampling mask, applied in one dimension.
- The CV histograms are plotted on the collective variable's own axis instead of a rescaled 0–1
  one, so they can be read directly against the free energy profiles. Their legend no longer
  covers the second panel's title.
- Basin markers and labels on the 3D plot are placed by projecting each minimum into the figure
  plane and drawing there, with a leader line. Drawn in 3D they disappeared behind the surface
  whenever a basin sat on the far side of it — on the bundled trajectory basin 4 was invisible —
  and `zorder` cannot fix it, because mplot3d sorts by depth. Each basin keeps the marker shape
  that identifies it in the legend.
- README rewritten: the method in eight steps, then the landscape in 2D and 3D. The full
  derivation moved to `docs/METHOD.md`, all six figures with their axes to `docs/FIGURES.md`, and
  the complete flag reference to `docs/OPTIONS.md`. Figures are referenced by relative path, so
  they render on any branch, and the exact command that reproduces every one of them from
  `inputs/` is given and verified byte-identical.

### Supported Python

Minimum raised to **3.10** (3.9 reached end of life in October 2025) and the package is now tested
through **3.14**. Dependency floors move to the first releases that ship wheels for modern
interpreters, so a clean install on a recent Python never falls back to building an old one:

| | before | now |
|---|---|---|
| Python | `>=3.8`, declared but untested above 3.10 | `>=3.10`, tested on 3.10–3.14 |
| numpy | `1.23.5` pinned in requirements | `>=1.26` |
| scipy | `1.10.1` pinned | `>=1.11` |
| matplotlib | `3.7.4` pinned | `>=3.8` |
| imageio | `2.34.0` pinned | `>=2.34` |
| joblib | `1.3.2` pinned | `>=1.3` |

Verified against numpy 2.5, scipy 1.18 and matplotlib 3.11 on Python 3.14, with no deprecation
warnings. A GitHub Actions workflow runs the suite on every supported version and, separately,
installs a built wheel on 3.14 and reproduces the bundled example end to end.

#### Other additions

- `free_energy_landscape/__init__.py`, so `from free_energy_landscape import FreeEnergyLandscape`
  works and the package is no longer an implicit namespace package.
- `tests/` — 59 pytest tests, including an analytical anchor: a two-Gaussian mixture with known
  weights must reproduce ΔG = −k_B·T·ln(w₁/w₂), both on the grid and between basin minima.
- `.gitignore` (covers `__pycache__/`, build artefacts, and `.pypirc`/`*.token`).
- The 3D landscape is now saved to `3D_landscape.png`. Previously it was only shown on screen,
  so Figure 5 of the README could not be reproduced by running the tool.
- This changelog.

### Data

The bundled example trajectory in `inputs/` was misaligned: both files are concatenated
trajectories whose frame counter restarted at different rows, so from row 23977 onward CV1 and CV2
were paired one frame apart. `proj2Out.txt` carried a duplicated row at the end of its first block
(frame 23977 repeating the value of frame 23976) and `proj1Out.txt` a surplus final frame. Both
rows were removed, leaving 51,642 rows with matching frame indices. The effect on the landscape
was up to 8.9 kJ/mol (RMS 1.2 kJ/mol) over the sampled region; the global minimum did not move.
The figures in `outputs/` were regenerated from the corrected input.

### Removed

- `calculate_density_for_chunk`. It only existed to implement the per-chunk refit described above.

## [1.0.2]

Initial published releases.
