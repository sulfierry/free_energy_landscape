[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.10850229.svg)](https://doi.org/10.5281/zenodo.10850229)

# Free Energy Landscape Analysis

This tool is designed to analyze and visualize the free energy landscape from molecular dynamics simulations. It calculates the free energy from collective variable (CV) data across simulation frames and generates insightful visualizations to help understand the system's behavior over time.


## Features

- **Free Energy vs. Collective Variable (CV) Visualization**: This tool expertly maps the free energy landscape across varying values of collective variables (CVs), such as distances or angles, facilitating a deeper understanding of the system's energetics. It highlights stable states and potential energy barriers, crucial for identifying key transition states and conformational changes.

- **CV Value vs. Frame Insights**: This feature charts the trajectory of collective variables (CVs) across simulation frames, uncovering the temporal dynamics of the system. It provides a clear depiction of how CV values evolve, revealing underlying patterns and fluctuations that characterize the molecular system's behavior over time.

- **Histogram Analysis of CV Distributions**: By generating detailed histograms of CV values, this tool underscores the frequency and distribution of collective variable states within the simulation. It accentuates prevalent states, offering a statistical perspective on the system's most favored conformations.

- **Free Energy Landscape (CV1 vs. CV2)**: This visualization combines two collective variables, offering a two-dimensional free energy landscape that encapsulates the interplay between different CVs. It enables a multifaceted analysis of how combined CVs contribute to the system's stability and transitions, enriching the understanding of complex molecular mechanisms.

- **3D Free Energy Landscape**: This feature constructs a three-dimensional representation of the free energy landscape, complemented by an animated GIF. This immersive visualization affords a dynamic and comprehensive view of the energy landscape, enhancing the interpretation of the system's energetics from multiple angles and dimensions.

- **Low Energy Points Identification**: This functionality ($--energy$ and $--discretize$) identifies and marks points falling below specific energy ($KJ/mol$) thresholds on the landscape. This is instrumental in pinpointing regions associated with stable conformations and significant for understanding the energetics that govern molecular stability and transitions.

## Required Libraries

This project relies on several key Python libraries for numerical computations, image processing, and plotting capabilities. Ensure you have the following libraries installed along with their specified versions to guarantee compatibility and proper functionality of the scripts:

- **NumPy** (1.23.5): A fundamental package for scientific computing with Python, providing support for large, multi-dimensional arrays and matrices, along with a collection of mathematical functions to operate on these arrays.
- **ImageIO** (2.34.0): A library for reading and writing a wide range of image, video, scientific, and volumetric data formats. It is used in this project for creating and manipulating images and GIFs.
- **Matplotlib** (3.7.4): A comprehensive library for creating static, animated, and interactive visualizations in Python. It is used for plotting the free energy landscapes.
- **SciPy** (1.10.1): An open-source Python library used for scientific computing and technical computing. It contains modules for optimization, linear algebra, integration, interpolation, special functions, FFT, signal and image processing, and more.
- **Joblib** (1.3.2): A versatile library optimized for fast and easy saving and loading of Python objects, making it particularly useful for large data that does not fit into memory. In the context of this project, Joblib is instrumental in facilitating efficient parallel computations and optimizing performance when processing and analyzing large datasets or performing computationally intensive tasks like kernel density estimations and free energy calculations across multiple processors.

To install these libraries, you can use the following command:

```bash
python3.10 -m venv env
```
```bash
pip install numpy==1.23.5 imageio==2.34.0 matplotlib==3.7.4 scipy==1.10.1 joblib==1.3.2
```
These dependencies also are listed in the `requirements.txt` file. To install them, run the following command in your terminal:

```bash
pip install -r requirements.txt
```
To install the free energy landscape:

```bash
pip install free-energy-landscape
```

## Usage

Ensure your data file is in the correct two-column format. Run the script with the path to your data files and optional arguments as needed:

```bash
free_energy_landscape path/to/cv1_data.txt path/to/cv2_data.txt
```
Optional Arguments:

```bash
   --temperature           [int]       Simulation temperature in Kelvin (default: 300K)
   --kb                    [float]     Boltzmann constant in kJ/(mol·K) (default: 8.314e-3)
   --energy                [float]     Energy threshold (KJ/mol). With the watershed method it
                                       only filters which basins are reported (those whose
                                       minimum lies below it) (default: None)
   --basin_method          [str]       'watershed' (default) assigns every frame to the minimum
                                       it descends to; 'connected' groups only the cells
                                       below --energy
   --basin_min_depth       [float]     Merge basins separated by a barrier shallower than this
                                       (KJ/mol). Default: kB*T, the thermal scale
   --basin_connectivity    [int]       'connected' only. 1 = 4 neighbours, 2 = 8 (default: 2)
   --basin_min_frames      [int]       Discard basins holding fewer sampled frames than this,
                                       as KDE artefacts (default: 1)
   --mask_min_count        [float]     A grid cell counts as sampled when the KDE expects at
                                       least this many frames inside it; the rest is left
                                       blank instead of extrapolated (default: 1)
   --cmap                  [str]       Matplotlib colormap, or 'classic' for the 1.x
                                       blue-to-red map (default: viridis_r)
   --dpi                   [int]       Resolution of the saved figures (default: 200)
   --discretize            [float]     LEGACY. Colours individual grid nodes by energy band
                                       instead of grouping them into basins (default: None)
   --bins_energy_histogram [int]       Bins for energy histogram (default: 100)
   --kde_bandwidth         [float]     Bandwidth for kernel density estimation (default: None)
   --max_energy            [float]     Upper cap (KJ/mol) for the plotted energy surface
                                       (default: none). Affects only the colour scale of the
                                       plots, never the .tsv values.
   --grid                  [int]       Grid resolution per axis for the landscape (default: 100)
   --n_jobs                [int]       Parallel workers for KDE evaluation (default: all CPUs).
                                       Results are identical for any value.
   --names                 [str] [str] Names for the collective variables (default: CV1, CV2)
   --gif_angles            [int]       Angles for 3D GIF rotation (default: 10)
   --gif_elevation         [int]       Elevation angle for the 3D GIF (default: 10)
   --gif_duration          [float]     Duration per frame in the GIF in seconds (default: 0.1)
```

### Reproducing the example figures

Every figure in this README comes from the two files in `inputs/`, produced by a single command.
Clone the repository and run, from its root:

```bash
free_energy_landscape inputs/proj1Out.txt inputs/proj2Out.txt --names "CV1 (Angle)" "CV2 (Distance)"
```

Everything else is default. The run is deterministic: the same input gives byte-identical output on
any machine and for any `--n_jobs`. It writes into the current directory:

| generated file | figure in this README |
|---|---|
| `Combined_Free_Energy_Profile_Normalized.png` | Figure 1 |
| `histograms_normalized_side_by_side.png` | Figure 2 |
| `cv_by_frame_combined_normalized.png` | Figure 3 |
| `Free_energy_landscape_with_discrete_points.png` | Figure 4 |
| `3D_landscape.png` | Figure 5 |
| `energy_landscape_3D.gif` | Figure 6 |
| `basins.tsv` | basin table below |
| `discrete_values_energy_frames.tsv` | per-frame energies and basin ids |

The copies under `outputs/` are the same files, renamed with a leading number so they sort in
reading order.

### Input format

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
files.

> **Upgrading from 1.x:** version 2.0.0 fixes a bug in the parallel KDE that made results depend on
> the number of CPU cores and could misidentify most low-energy frames. Results from 1.x differ and
> should be re-computed. See [CHANGELOG.md](CHANGELOG.md) for the measured impact.

## How the analysis works

Eight steps, from raw CV values to a table of metastable states. Nothing here is random or
iterative: the same input always gives the same answer.

### Step 1 — Density

The two CV columns give one point per frame in a 2D space. A Gaussian kernel density estimate
(`scipy.stats.gaussian_kde`) turns those points into a smooth probability density $\rho$, evaluated
on a regular grid (`--grid`, default 100×100).

A single kernel is fitted over the whole trajectory. Parallelism (`--n_jobs`) only splits the
*evaluation* of that one kernel, so the number of cores never changes the result.

### Step 2 — Free energy

Boltzmann inversion converts density into free energy, then shifts the global minimum to zero:

$$G(x,y) = -k_B T \ln \rho(x,y), \qquad G \leftarrow G - \min G$$

Every energy the tool reports is relative to the deepest point of the landscape.

### Step 3 — Sampling mask

A KDE returns a number everywhere, including regions that contain no frame at all. Those numbers
are extrapolation, not measurement — on the bundled trajectory they reach ~1700 kJ/mol.

A grid cell counts as sampled only when the estimate expects at least `--mask_min_count` frames
inside it:

$$\rho \cdot A_{\text{cell}} \cdot N \ \geq\ \text{mask\\_min\\_count}$$

Everything else is excluded from the analysis and from the colour scale. This is why `--max_energy`
now defaults to `none`: there is no longer a runaway range to clip.

### Step 4 — Steepest descent

Each cell looks at its 8 neighbours — the full 3×3 block minus the centre — and points at the one
with the **lowest** energy. If no neighbour is lower, the cell points at itself: it is a local
minimum.

Following those pointers to their fixed point tells which minimum each cell drains into. That
minimum is its **basin of attraction**. Since a step only accepts a strictly lower neighbour, no
cycle can form and the iteration always terminates.

Eight rather than four neighbours matters: with only the four axis directions, descent cannot move
diagonally, so a diagonal valley makes it zigzag and stall. On the bundled trajectory that produces
25 local minima instead of 12 — thirteen of them artefacts of the grid, not of the molecule.

### Step 5 — Saddles

Two neighbouring basins meet along a ridge. The **saddle** between them is the lowest point of that
ridge — the cheapest way across:

$$\text{saddle}(A,B) = \min_{\substack{x \in A,\ y \in B \\ x,\,y \text{ adjacent}}} \max\big(G(x),\,G(y)\big)$$

### Step 6 — Depth

$$\text{depth}(A) = \min_B \text{saddle}(A,B)\ -\ \min_{x \in A} G(x)$$

How far a basin has to climb to escape. This is its topological persistence, and it is what decides
whether a dip is a state or a wrinkle. Unlike an energy cut-off, it does not depend on the grid
resolution.

### Step 7 — Merging

Basins shallower than `--basin_min_depth` are absorbed, one at a time, by the neighbour across
their lowest saddle.

The cut defaults to automatic: sort all depths and place it at the **largest step** in that list.
Real structure and estimator ripple usually sit on opposite sides of a wide gap. The automatic cut
never exceeds $k_B T$, because a barrier taller than the thermal scale is real by definition.

What it chose, and what it merged, is printed on every run:

```
Watershed found 12 local minima; --basin_min_depth 1.01 kJ/mol (auto) merged 7 of them.
Barrier depths (deepest first): 2.10, 1.81, 1.59, 1.43, 1.41, 0.61, 0.47, 0.27, 0.17, ...
Largest step is 1.41 -> 0.61 kJ/mol: --basin_min_depth between them keeps 6 basins.
```

Pass `--basin_min_depth` to override; `0` keeps every local minimum.

### Step 8 — Frames and representatives

Each frame falls in a grid cell and inherits its basin id. The **representative** of a basin is the
frame of lowest free energy inside it — an actual conformation from the trajectory, ready to
extract, not a grid node.

Basins are numbered by increasing energy, so basin 1 is always the global minimum.

### A worked example, in one dimension

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

### Reading the output

```
  #        CV1        CV2       dG    depth   frames    pop%  rep frame
-----------------------------------------------------------------------
  1     58.181      7.578     0.00     3.82    36447   70.58      10711
  2     40.777      5.786     2.41     1.41     5201   10.07      13367
  3     88.408     11.761     4.30     1.59     9263   17.94       4364
  4    114.971     22.219    13.77        -      164    0.32      23523
  5    109.475     32.975    13.92     2.10       98    0.19      23786
  6    116.803     31.780    14.59     1.43       43    0.08      23704
```

- **CV1, CV2** — where the basin's minimum sits.
- **dG** — its free energy, relative to the global minimum.
- **depth** — how far it must climb to reach a neighbour. A dash means it has no neighbour: it sits
  on its own island of sampled space, so no saddle exists.
- **frames**, **pop%** — how much of the trajectory drains into it.
- **rep frame** — the conformation to extract as its representative.

Written to `basins.tsv`, with the basin id of every frame added as a `cluster` column in
`discrete_values_energy_frames.tsv` (`0` = outside the sampled region).

### Other grouping modes

`--energy` does not change the partition under watershed; it only filters which basins are
reported, keeping those whose minimum lies below it.

> `--basin_method connected` restores the earlier grouping: connected components of $G \le$
> `--energy`. It depends on the threshold chosen and leaves every frame above it unassigned.
>
> `--discretize` restores the original v1 rendering, which coloured individual **grid nodes** by
> energy band. Those nodes are cells of the KDE grid, not conformations, and two nodes in the same
> band could sit in regions with no path between them. Kept only for backwards compatibility.

Density-based clustering such as DBSCAN is deliberately **not** used here. It groups by local
density, which is the quantity the free energy is already the logarithm of, so it re-derives less
than what Step 2 gives for free: no barrier heights, no ΔG, parameters (`eps`, `minPts`) with no
physical unit, and low-density points — the transition states — discarded as noise.

## Theoretical Background

The free energy landscape is a conceptual and computational tool used to understand the energetics and dynamics of molecular systems. It is particularly useful for interpreting molecular dynamics simulations, where it condenses a trajectory of millions of coordinates into the energetics of the processes the system actually undergoes: protein folding, chemical reactions, and phase transitions.

### Collective Variables (CVs)

Collective Variables (CVs) are a set of coordinates that describe the macroscopic state of a system. They are used to reduce the complexity of molecular systems by focusing on the relevant degrees of freedom. Examples include the analysis of principal components (PCA), the distance between two atoms, angles, dihedrals, and more complex descriptors.

Following PCA, the use of distance and angle as CVs in this study serves as an example to illustrate the tool's capabilities. However, it's important to note that the input can be any file containing two columns: the first for frames and the second for the value of the collective variable. This flexibility allows the tool to be applicable to a wide range of studies involving different types of collective variables.

These representations allow for a simplified description of the system's state, facilitating the study of its behavior and properties through the manipulation of a reduced set of variables rather than the full set of atomic coordinates.

### Kernel Density Estimation (KDE)

This statistical method is crucial for estimating the probability density function (PDF) of a dataset without assuming any predefined distribution shape. In the context of molecular dynamics and simulations, it allows us to visualize and analyze the distribution of free energy across different states or configurations defined by collective variables (CVs).

The `gaussian_kde` leverages a Gaussian (normal) distribution, placing it at each point in the dataset and summing these distributions to approximate the overall data's PDF. This technique is adept at capturing the underlying structure of the data, providing a smooth, continuous representation of the free energy landscape. The smoothness of the KDE is controlled by the bandwidth parameter, which determines the width of the Gaussian kernels used.

The KDE for a set of ($n$) points (${x_i}$) can be mathematically represented as: 

$$\hat{f}(x) = \frac{1}{n \cdot h} \sum_{i=1}^{n} K\left( \frac{x - x_i}{h} \right)$$ 

where $\hat{f}(x)$ is the estimated density at point ($x$), ($K$) is the kernel function (e.g., Gaussian), and ($h$) is the bandwidth, a parameter that controls the smoothness of the density estimate. 

The conversion from the estimated density to free energy is typically done using the relation: 

$$G(x) = -k_B T \ln(\hat{f}(x))$$ 

where $G(x)$ is the free energy at point $x$, $k_B$ is the Boltzmann constant, and $T$ is the temperature. 

The Gaussian KDE method provides a sophisticated approach to model the complex free energy landscapes encountered in molecular dynamics studies. It enables researchers to visualize the distribution of energy states without the constraints of parametric models, offering insights into molecular stability, transitions, and the energetics of molecular interactions.

## Visualizations Generated

The `FreeEnergyLandscape` class generates three key visualizations to aid in the analysis of the molecular system's free energy landscape. Each figure provides unique insights into the system's thermodynamic and kinetic behavior.

### Figure 1:  Normalized Free Energy Profile Comparison (CV1 and CV2)

This updated figure now integrates the normalized free energy profiles for both CV1 (Angle) and CV2 (Distance) in a single, unified visualization. This enhancement allows for a direct comparison between the two collective variables, illuminating their unique energy landscapes and pinpointing low-energy regions across specified energy threshold ranges.

![Alt text da image](https://github.com/sulfierry/free_energy_landscape/blob/main/outputs/1_Combined_Free_Energy_Profile_Normalized.png?raw=true)

- **X-Axis**: Displays both CV1 and CV2, with CV1 representing an angle and CV2 denoting distance in Ångströms, facilitating a comprehensive view of the system's energetics.
- **Y-Axis**: Represents the normalized free energy, derived from the probability distributions of CV1 and CV2 using the Boltzmann equation.
- **Visualization**: A comparative plot of free energy versus both CV1 and CV2, highlighting key energy barriers and minima corresponding to the molecular system's different conformational states.

### Figure 2: Histograms of Normalized CV1 and CV2 Values Side by Side

![Alt text da image](https://github.com/sulfierry/free_energy_landscape/blob/main/outputs/2_histograms_normalized_side_by_side.png?raw=true)

- **Visualization**: Two histograms placed side by side, one for CV1 (Angle) and the other for CV2 (Distance), each showing the normalized value distribution. This format allows for an easy comparison of the two variables' behaviors and predominant states within the simulation.


### Figure 3: CV Value vs. Frame for CV1 and CV2

This figure innovatively plots the values of CV1 (Angle) and CV2 (Distance) against simulation frames in the same image. It offers insights into how the values of these collective variables evolve over time, shedding light on the dynamic aspects of the system's behavior.

![Alt text da image](https://github.com/sulfierry/free_energy_landscape/blob/main/outputs/3_cv_by_frame_combined_normalized.png?raw=true)


- **X-Axis**: Represents the simulation frames, providing a temporal dimension to the analysis.
- **Y-Axis**: Displays values of CV1 and CV2, allowing for a direct observation of how each collective variable changes throughout the simulation.
- **Visualization**: A combined plot that tracks the fluctuation of both CV1 and CV2 over the course of the simulation, highlighting trends, transitions, and stability regions within the molecular system.

### Figure 4: 2D Free Energy Landscape (CV1 vs. CV2)

The third figure combines both CV1 and CV2 to produce a two-dimensional free energy landscape, offering a comprehensive view of the system's energetics over the considered collective variables.

![Alt text da image](https://github.com/sulfierry/free_energy_landscape/blob/main/outputs/4_Free_energy_landscape.png?raw=true)

- **X-Axis**: CV1, representing an angle, with the same labeling as in Figure 1 (left).
- **Y-Axis**: CV2, representing a distance, with the same labeling as in Figure 1 (right).
- **Visualization**: A heatmap showing the free energy levels across the CV1 and CV2 space. The color gradient represents different energy levels, with cooler colors indicating low-energy regions (minima) and warmer colors highlighting high-energy barriers. This visualization is crucial for identifying transition states, stable conformations, and understanding the molecular system's behavior under various conditions.

### Figure 5: 3D Free Energy Landscape

This figure provides a three-dimensional visualization of the free energy landscape, combining both collective variables, CV1 and CV2, along with the calculated free energy values to offer a dynamic perspective on the system's energetics.

![Alt Text](https://github.com/sulfierry/free_energy_landscape/blob/main/outputs/5_3D_landscape.png?raw=true)


- **X-Axis**: CV1, representing an angle, with the same labeling as in Figure 4.
- **Y-Axis**: CV2, representing a distance, with the same labeling as in Figure 4.
- **Z-Axis**: Shows the free energy in kJ/mol, calculated from the probability distribution of CV1 and CV2 using the Boltzmann equation.

**Visualization**: A 3D surface plot illustrates the variations in free energy across the space defined by CV1 and CV2. The color gradient enhances the visualization of energy barriers and minima, aiding in the identification of stable conformations and transition states. The 3D and animated visualizations of the free energy landscape extend the analysis beyond two dimensions, offering a richer and more nuanced understanding of the system's energetics. By exploring the landscape in three dimensions, researchers can better identify and analyze the regions of interest, such as low-energy conformations and high-energy transition states, which are vital for deciphering the molecular mechanisms underlying biological processes and material behaviors.

## Interpretation

Together, these figures provide a multi-faceted view of the molecular system's free energy landscape. Analyzing these visualizations helps in understanding how variations in critical structural parameters (angles and distances) influence the stability and dynamics of the system, which is vital for predicting reaction pathways, designing drugs, and engineering materials with desired properties.
