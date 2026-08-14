# Command-line reference

```bash
free_energy_landscape path/to/cv1.txt path/to/cv2.txt [options]
```

Every option is optional. `-h` / `--help` prints this list.

## Physical constants

| option | type | default | meaning |
|---|---|---|---|
| `--temperature` | float | `300` | Simulation temperature, in Kelvin |
| `--kb` | float | `8.314e-3` | Boltzmann constant, in kJ/(mol·K) |

## Density estimate

| option | type | default | meaning |
|---|---|---|---|
| `--kde_bandwidth` | float | `none` | Bandwidth of the Gaussian kernel. Default uses Scott's rule |
| `--grid` | int | `100` | Grid resolution per axis for the landscape |
| `--mask_min_count` | float | `1` | A grid cell counts as sampled when the estimate expects at least this many frames inside it. The rest is left out instead of extrapolated |
| `--n_jobs` | int | all CPUs | Parallel workers for kernel evaluation. Results are identical for any value |

## Basins

| option | type | default | meaning |
|---|---|---|---|
| `--basin_method` | str | `watershed` | `watershed` assigns every frame to the minimum it descends to; `connected` groups only the cells below `--energy` |
| `--basin_min_depth` | float | auto | Merge basins separated by a barrier shallower than this, in kJ/mol. Automatic cut goes at the largest step in the depth spectrum, capped at $k_B T$. `0` keeps every local minimum |
| `--basin_min_frames` | int | `1` | Discard basins holding fewer sampled frames than this, as estimator artefacts |
| `--basin_connectivity` | int | `2` | `connected` only. `1` = 4 neighbours, `2` = 8 |
| `--energy` | float | `none` | Energy threshold, in kJ/mol. Under `watershed` it only filters which basins are reported — those whose minimum lies below it. Under `connected` it defines the grouping |

## Plots

| option | type | default | meaning |
|---|---|---|---|
| `--names` | str str | `CV1 CV2` | Axis labels for the two collective variables |
| `--cmap` | str | `viridis_r` | Any matplotlib colormap name, or `classic` for the 1.x blue-to-red map |
| `--dpi` | int | `200` | Resolution of the saved figures |
| `--max_energy` | float | `none` | Upper cap for the plotted surface, in kJ/mol. Affects only the colour scale, never the `.tsv` values |
| `--bins_energy_histogram` | int | `100` | Bins for the 1D energy histograms |
| `--xlim_inf`, `--xlim_sup` | float | `none` | Fix the x-axis range. By default the axes zoom to the sampled region |
| `--ylim_inf`, `--ylim_sup` | float | `none` | Fix the y-axis range |

## 3D animation

| option | type | default | meaning |
|---|---|---|---|
| `--gif_angles` | int | `10` | Number of frames around the rotation |
| `--gif_elevation` | int | `10` | Elevation angle of the camera |
| `--gif_duration` | float | `0.1` | Seconds per frame |

## Legacy

| option | type | default | meaning |
|---|---|---|---|
| `--discretize` | float | `none` | Restores the v1 rendering, colouring individual grid nodes by energy band. Those nodes are cells of the KDE grid, not conformations, and two nodes in the same band can sit in regions with no path between them. Kept for backwards compatibility only |

## Input format

Each CV file is a two-column text file: the frame index in column 0 and the CV value in column 1.

```
0 83.81493927422635
1 82.63629337277925
2 79.48816566948845
```

Both files must describe the **same** trajectory — same number of rows, and the same frame index on
every row. Otherwise CV1 of one conformation would be paired with CV2 of another, and the landscape
would describe a molecule that never existed. This is checked at load time and reported with the
offending row numbers.

Concatenated trajectories are fine, as long as the frame counter restarts at the same row in both
files. When it does not, the error names the rows where each file restarts.

## Output files

Written to the current directory.

| file | contents |
|---|---|
| `Combined_Free_Energy_Profile_Normalized.png` | 1D free energy profile of each CV |
| `histograms_normalized_side_by_side.png` | Distribution of each CV |
| `cv_by_frame_combined_normalized.png` | Both CVs against frame number |
| `Free_energy_landscape_with_discrete_points.png` | 2D landscape with basins |
| `3D_landscape.png` | 3D surface |
| `energy_landscape_3D.gif` | Rotating 3D surface |
| `basins.tsv` | One row per basin |
| `discrete_values_energy_frames.tsv` | Per-frame energy and basin id |

### `basins.tsv`

`cluster`, `cv1_min`, `cv2_min`, `G_min_kJ_mol`, `depth_kJ_mol`, `n_frames`,
`population_fraction`, `rep_frame`, `rep_cv1`, `rep_cv2`, `rep_energy_kJ_mol`, `area_cells`

`depth_kJ_mol` is `inf` when the basin has no neighbour to be separated from.

### `discrete_values_energy_frames.tsv`

`frame`, `cv1`, `cv2`, `energy`, `cluster` — sorted by energy. `cluster = 0` means the frame lies
outside the sampled region. With `--energy`, only frames below the threshold are written.

## Upgrading from 1.x

Version 2.0.0 fixes a bug in the parallel KDE that made results depend on the number of CPU cores
and could misidentify most low-energy frames. Results from 1.x differ and should be re-computed.
See [CHANGELOG.md](../CHANGELOG.md) for the measured impact.
