# Figures

All six come from one command, run on the two files in `inputs/`:

```bash
free_energy_landscape inputs/proj1Out.txt inputs/proj2Out.txt --names "CV1 (Angle)" "CV2 (Distance)"
```

CV1 is an angle and CV2 a distance in ångströms, used here only as an example — any pair of
collective variables works. See [OPTIONS.md](OPTIONS.md) for the flags that change them.

## 1. Free energy profile of each CV

![Free energy profile](../outputs/1_Combined_Free_Energy_Profile_Normalized.png)

Boltzmann inversion of the 1D histogram of each collective variable, on a shared energy axis.

- **x** — the collective variable, in its own units
- **y** — free energy in kJ/mol, relative to the lowest bin of that variable

Minima are the preferred values of that CV; the rise between them is the barrier along that
coordinate alone. A break in the curve is a range where no frame was observed — an empty histogram
bin carries no energy, so none is drawn.

Read this together with the 2D landscape: a barrier that appears here may vanish there, because
projecting onto one coordinate can hide a path that goes around it in the other.

## 2. Distribution of each CV

![Histograms](../outputs/2_histograms_normalized_side_by_side.png)

The raw counts that Figure 1 turns into energy, on the same axes.

- **x** — the collective variable, in its own units
- **y** — frequency as a percentage of the tallest bin across both variables

Useful as a sampling check: a peak that rests on a handful of frames will show up here as a thin,
isolated bar, and any energy read off it should be treated with care.

## 3. CVs against frame

![CV by frame](../outputs/3_cv_by_frame_combined_normalized.png)

Both collective variables over the trajectory, each rescaled to 0–100 so they share one axis.

- **x** — position in the trajectory
- **y** — each CV rescaled to its own range

Shows when transitions happen and whether sampling is converged: a landscape built from a
trajectory that never leaves one basin describes that basin, not the system. A step discontinuity
marks the join between concatenated trajectories.

## 4. Free energy landscape

![2D landscape](../outputs/4_Free_energy_landscape.png)

The two collective variables combined — the main result.

- **x, y** — the two collective variables
- **colour** — free energy in kJ/mol, bright for low, dark for high
- **dashed line** — edge of the sampled region; beyond it no energy is claimed
- **white lines** — watershed ridges separating basins
- **numbered callouts** — each basin minimum, with depth, population and representative frame in
  the legend

This is where transition states, stable conformations and the paths between them are read off.

## 5. Landscape in three dimensions

![3D landscape](../outputs/5_3D_landscape.png)

The same surface as relief.

- **x, y** — the two collective variables
- **z** and **colour** — free energy in kJ/mol

The depth of each basin and the height of the barriers between them become directly comparable,
which is harder to judge from colour alone. The surface stops at the sampled region.

## 6. Rotating view

![3D animation](../outputs/6_energy_landscape_3D.gif)

The 3D surface across a full rotation, so features hidden behind a ridge from one angle become
visible from another. Controlled by `--gif_angles`, `--gif_elevation` and `--gif_duration`.

## Interpretation

Together these give a layered view: Figures 1–3 describe each coordinate separately and check that
the sampling supports any conclusion, Figures 4–6 show how the two combine. Reading them side by
side is what connects a change in a structural parameter to the stability of the states it
produces — the basis for predicting reaction pathways, designing ligands, and engineering materials.
