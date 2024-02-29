[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.10689690.svg)](https://doi.org/10.5281/zenodo.10689690)


# Free Energy Landscape Analysis

This tool is designed to analyze and visualize the free energy landscape from molecular dynamics simulations. It calculates the free energy from collective variable (CV) data across simulation frames and generates insightful visualizations to help understand the system's behavior over time.


## Features

- **Free Energy vs. Collective Variable (CV) Visualization**: This tool expertly maps the free energy landscape across varying values of collective variables (CVs), such as distances or angles, facilitating a deeper understanding of the system's energetics. It highlights stable states and potential energy barriers, crucial for identifying key transition states and conformational changes.

- **CV Value vs. Frame Insights**: This feature charts the trajectory of collective variables (CVs) across simulation frames, uncovering the temporal dynamics of the system. It provides a clear depiction of how CV values evolve, revealing underlying patterns and fluctuations that characterize the molecular system's behavior over time.

- **Histogram Analysis of CV Distributions**: By generating detailed histograms of CV values, this tool underscores the frequency and distribution of collective variable states within the simulation. It accentuates prevalent states, offering a statistical perspective on the system's most favored conformations.

- **Free Energy Landscape (CV1 vs. CV2)**: This visualization combines two collective variables, offering a two-dimensional free energy landscape that encapsulates the interplay between different CVs. It enables a multifaceted analysis of how combined CVs contribute to the system's stability and transitions, enriching the understanding of complex molecular mechanisms.

- **3D Free Energy Landscape**: This feature constructs a three-dimensional representation of the free energy landscape, complemented by an animated GIF. This immersive visualization affords a dynamic and comprehensive view of the energy landscape, enhancing the interpretation of the system's energetics from multiple angles and dimensions.

- **Low Energy Points Identification**: This functionality ($--energy$) identifies and marks points falling below specific energy ($KJ/mol$) thresholds on the landscape. This is instrumental in pinpointing regions associated with stable conformations and significant for understanding the energetics that govern molecular stability and transitions.

## Required Libraries

This project relies on several key Python libraries for numerical computations, image processing, and plotting capabilities. Ensure you have the following libraries installed along with their specified versions to guarantee compatibility and proper functionality of the scripts:

- **NumPy** (1.23.5): A fundamental package for scientific computing with Python, providing support for large, multi-dimensional arrays and matrices, along with a collection of mathematical functions to operate on these arrays.
- **ImageIO** (2.34.0): A library for reading and writing a wide range of image, video, scientific, and volumetric data formats. It is used in this project for creating and manipulating images and GIFs.
- **Matplotlib** (3.7.4): A comprehensive library for creating static, animated, and interactive visualizations in Python. It is used for plotting the free energy landscapes.
- **SciPy** (1.10.1): An open-source Python library used for scientific computing and technical computing. It contains modules for optimization, linear algebra, integration, interpolation, special functions, FFT, signal and image processing, and more.

To install these libraries, you can use the following command:

```bash
pip install numpy==1.23.5 imageio==2.34.0 matplotlib==3.7.4 scipy==1.10.1
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
   --energy                [float]     Energy (KJ/mol), single value (default: None)
   --discretize            [float]     Discrete value (KJ/mol) for energy (default: None)
   --bins_energy_histogram [int]       Bins for energy histogram (default: 100)
   --kde_bandwidth         [float]     Bandwidth for kernel density estimation (default: None)
   --names                 [str] [str] Names for the collective variables (default: CV1, CV2)
   --gif_angles            [int]       Angles for 3D GIF rotation (default: 10)
   --gif_elevation         [int]       Elevation angle for the 3D GIF (default: 10)
   --gif_duration          [float]     Duration per frame in the GIF in seconds (default: 0.1)
```

Example with optional arguments:

```bash
free_energy_landscape path/to/cv1_data.txt path/to/cv2_data.txt --names "CV1 (Angle)" "CV2 (Distance)" --energy 3
```


## Theoretical Background

The free energy landscape is a conceptual and computational tool used to understand the energetics and dynamics of molecular systems. It is particularly useful in the study of complex processes such as protein folding, chemical reactions, and phase transitions.

### Collective Variables (CVs)

Collective Variables (CVs) are a set of coordinates that describe the macroscopic state of a system. They are used to reduce the complexity of molecular systems by focusing on the relevant degrees of freedom. Examples include the analysis of principal components (PCA), the distance between two atoms, angles, dihedrals, and more complex descriptors. Below are examples of commonly used CVs and their mathematical formulations.

The inclusion of Principal Component Analysis (PCA) as an example of CVs is crucial for understanding the dimensional reduction in complex systems. PCA identifies the directions (principal components) along which the variance in the data is maximized. In the context of CVs, the first two principal components, denoted as $PC_1$ and $PC_2$, serve as collective variables that capture the most significant modes of variation within the system. Mathematically, PCA transforms the original data into a new set of variables, the principal components, which are linear combinations of the original variables with weights given by the eigenvectors of the data's covariance matrix.

PCA is an algebraic procedure that transforms a number of (possibly) correlated variables into a (smaller) number of uncorrelated variables called principal components. This transformation is defined in such a way that the first principal component has the highest possible variance (it accounts for as much of the variability in the data as possible), and each succeeding component in turn has the highest variance possible under the constraint that it is orthogonal to (i.e., uncorrelated with) the preceding components. The components are orthogonal because they are the eigenvectors of the data's covariance matrix $\mathbf{C}$, which is a symmetric matrix.

So, $\mathbf{X}$ is a $n \times p$ data matrix representing $n$ observations with $p$ variables, the goal of PCA is to transform $\mathbf{X}$ into a new coordinate system through a linear transformation that maps $\mathbf{X}$ onto a new set of variables, the principal components (PCs), which are linear combinations of the original variables. This is achieved by eigendecomposition of the covariance matrix $\mathbf{C} = \frac{1}{n-1} \mathbf{X}^T\mathbf{X}$, or more commonly through singular value decomposition (SVD) of $\mathbf{X}$.

Geometrically, PCA seeks the line (in the case of a single principal component) or hyperplane (in the case of multiple principal components) that best captures the distribution of the data in a high-dimensional space. The first principal component is the direction along which the projection of the data points has the largest variance. This direction corresponds to the eigenvector of the covariance matrix associated with the largest eigenvalue. The second principal component is orthogonal to the first and represents the direction of the next highest variance, corresponding to the eigenvector associated with the second largest eigenvalue, and so on.

1. **Principal Component Analysis (PCA)**: Can be expressed as:
$$PC_k = \sum_{i=1}^{p} a_{ki}X_i$$
where $PC_k$ is the $k_{th}$ principal component, $X_i$ are the original variables, and $a_{ki}$ are the coefficients (loadings) for the $k_{th}$ principal component, given by the $k_{th}$ eigenvector of the covariance matrix $\mathbf{C}$. The transformation can be represented in matrix form as:
$$\mathbf{PC} = \mathbf{X}\mathbf{A}$$
where $\mathbf{PC}$ is the matrix of principal components, $\mathbf{X}$ is the original data matrix (centered or standardized, if necessary), and $\mathbf{A}$ is the matrix whose columns are the eigenvectors of $\mathbf{C}$.

This algebraic and geometric perspective highlights the essence of PCA as a method for identifying the directions (in the space of the original variables) that maximize the variance of the projected data, thereby reducing its dimensionality while preserving as much of the original data's variation as possible.

2. **Angle**: The angle $\theta$ formed by three atoms $i$, $j$, and $k$, where $j$ is the vertex, can be calculated using the dot product:

   First, determine the vectors $\vec{r}_ji$ and $\vec{r}_jk$:
   $$\vec{r}_ji = \vec{r}_i - \vec{r}_j$$
   $$\vec{r}_jk = \vec{r}_k - \vec{r}_j$$

   Then, calculate the angle $\theta$ as:
   $$\cos(\theta) = \frac{\vec{r}_ji \cdot \vec{r}_jk}{|\vec{r}_ji| |\vec{r}_jk|}$$
   $$\theta = \arccos\left(\frac{\vec{r}_ji \cdot \vec{r}_jk}{|\vec{r}_ji| |\vec{r}_jk|}\right)$$

   - $\vec{r}_ji$ and $\vec{r}_jk$ are vectors pointing from atom $j$ to atoms $i$ and $k$, respectively.
   - $\cdot$ denotes the dot product between the vectors $\vec{r}_ji$ and $\vec{r}_jk$.
   - $|\vec{r}_ji|$ and $|\vec{r}_jk|$ represent the magnitudes of the vectors $\vec{r}_ji$ and $\vec{r}_jk$, respectively.
   - $\arccos$ is the inverse cosine function, used to find the angle $\theta$ from the cosine value.
  

3. **Distance**: The distance $d$ between two atoms $i$ and $j$ with positions $\vec{r}_i$ and $\vec{r}_j$ can be calculated using the Euclidean distance formula:

   First, determine the position vectors of atoms $i$ and $j$:
   $$\vec{r}_i = (x_i, y_i, z_i)$$
   $$\vec{r}_j = (x_j, y_j, z_j)$$

   Then, calculate the distance $d$ as:
   $$d = |\vec{r}_i - \vec{r}_j| = \sqrt{(x_i - x_j)^2 + (y_i - y_j)^2 + (z_i - z_j)^2}$$

   - $\vec{r}_i$ and $\vec{r}_j$ are the position vectors of atoms $i$ and $j$, respectively.
   - $(x_i, y_i, z_i)$ and $(x_j, y_j, z_j)$ denote the Cartesian coordinates of atoms $i$ and $j$.
   - $|\vec{r}_i - \vec{r}_j|$ represents the magnitude of the vector difference between $\vec{r}_i$ and $\vec{r}_j$, giving the direct distance between the two atoms.

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

### Boltzmann distribution

The Boltzmann distribution relates the energy of a state to its probability in a canonical ensemble, this provides a fundamental link between the microstate probabilities of a thermodynamic system and its macroscopic properties. It expresses the probability $P$ of a system being in a state with a certain free energy $\Delta G$ at a specific temperature $T$:

$$P \propto e^{-\frac{\Delta G}{k_B T}}$$

where:
- $P$ is the probability of observing the state,
- $\Delta G$ represents the free energy difference of that state relative to a chosen reference state,
- $k_B$ is the Boltzmann constant, a fundamental physical constant that relates energy scales to temperature,
- $T$ is the absolute temperature of the system, measured in Kelvin.

### Normalization of Probability

For the concept of probability to be meaningful in the context of statistical mechanics, it's imperative that the probabilities of all conceivable states sum to unity. This requirement ensures that the predicted behavior encompasses all possible configurations of the system. Consider a state that represents the maximum free energy, denoted as $\Delta G_{\text{max}}$, its corresponding maximum probability, $P_{\text{max}}$, can be analogously expressed through the Boltzmann factor.

### From Proportionality to Quantitative Relationship

To quantitatively relate $\Delta G$ with $P$ and $P_{\text{max}}$, we proceed as follows:

1. Express $P$ as a function of $\Delta G$:

   $$P = e^{-\frac{\Delta G}{k_B T}}$$

2. Similarly, express $P_{\text{max}}$ as a function of $\Delta G_{\text{max}}$:

   $$P_{\text{max}} = e^{-\frac{\Delta G_{\text{max}}}{k_B T}}$$

Taking the natural logarithm on both sides of these expressions yields:

$$\ln(P) = -\frac{\Delta G}{k_B T}$$

$$\ln(P_{\text{max}}) = -\frac{\Delta G_{\text{max}}}{k_B T}$$

Subtracting the equation for $\ln(P_{\text{max}})$ from the equation for $\ln(P)$ removes the dependency on the reference state $\Delta G_{\text{max}}$, leading to:

$$\ln(P) - \ln(P_{\text{max}}) = -\frac{\Delta G}{k_B T} + \frac{\Delta G_{\text{max}}}{k_B T}$$

Given the interest in determining the free energy difference $\Delta G$ relative to the state with maximum probability $P_{\text{max}}$, we can rearrange this relationship to solve explicitly for $\Delta G$:

$$\Delta G = -k_B T (\ln(P) - \ln(P_{\text{max}}))$$

The Gibbs free energy provides a measure of the amount of "useful work" that can be obtained from a thermodynamic system in a process at constant temperature and pressure. Here, $\Delta G$ is expressed as a function of the difference in the logarithms of the probabilities of finding the system in any state versus the state of highest probability (or lowest free energy), multiplied by the system's absolute temperature and the Boltzmann constant. This relationship shows how the free energy difference between two states can be calculated from their relative probabilities, providing a direct bridge between statistical thermodynamics and experimental observations or computational simulations.

This relationship allows us to convert a histogram of CV values into a free energy surface and this formulation also adjusts the free energy values so that the minimum energy associated with the highest probability $P_{\text{max}}$ is set to 0. This expression quantitatively links the probability distribution of states within a thermodynamic system to their respective free energy differences, providing a foundation for analyzing the system's behavior at the molecular level and this approach is useful for highlighting the relative differences in free energy between different states or conformations in an MD simulation, facilitating the identification of free energy minima and the relative comparison between different states. By setting the minimum of free energy to 0, this create a clear reference point to evaluate the relative stability of other states compared to the most stable state.

## Visualizations Generated

The `FreeEnergyLandscape` class generates three key visualizations to aid in the analysis of the molecular system's free energy landscape. Each figure provides unique insights into the system's thermodynamic and kinetic behavior.

### Figure 1:  Normalized Free Energy Profile Comparison (CV1 and CV2)

This updated figure now integrates the normalized free energy profiles for both CV1 (Angle) and CV2 (Distance) in a single, unified visualization. This enhancement allows for a direct comparison between the two collective variables, illuminating their unique energy landscapes and pinpointing low-energy regions across specified energy threshold ranges.

![Alt text da image](https://github.com/sulfierry/MolecularModeling/blob/main/FreeEnergyLandscape/1_Combined_Free_Energy_Profile_Normalized.png)

- **X-Axis**: Displays both CV1 and CV2, with CV1 representing an angle and CV2 denoting distance in Ångströms, facilitating a comprehensive view of the system's energetics.
- **Y-Axis**: Represents the normalized free energy, derived from the probability distributions of CV1 and CV2 using the Boltzmann equation.
- **Visualization**: A comparative plot of free energy versus both CV1 and CV2, highlighting key energy barriers and minima corresponding to the molecular system's different conformational states.

### Figure 2: Histograms of Normalized CV1 and CV2 Values Side by Side

![Alt text da image](https://github.com/sulfierry/MolecularModeling/blob/main/FreeEnergyLandscape/2_histograms_normalized_side_by_side.png)

- **Visualization**: Two histograms placed side by side, one for CV1 (Angle) and the other for CV2 (Distance), each showing the normalized value distribution. This format allows for an easy comparison of the two variables' behaviors and predominant states within the simulation.


### Figure 3: CV Value vs. Frame for CV1 and CV2

This figure innovatively plots the values of CV1 (Angle) and CV2 (Distance) against simulation frames in the same image. It offers insights into how the values of these collective variables evolve over time, shedding light on the dynamic aspects of the system's behavior.

![Alt text da image](https://github.com/sulfierry/MolecularModeling/blob/main/FreeEnergyLandscape/3_cv_by_frame_combined_normalized.png)


- **X-Axis**: Represents the simulation frames, providing a temporal dimension to the analysis.
- **Y-Axis**: Displays values of CV1 and CV2, allowing for a direct observation of how each collective variable changes throughout the simulation.
- **Visualization**: A combined plot that tracks the fluctuation of both CV1 and CV2 over the course of the simulation, highlighting trends, transitions, and stability regions within the molecular system.

### Figure 4: 2D Free Energy Landscape (CV1 vs. CV2)

The third figure combines both CV1 and CV2 to produce a two-dimensional free energy landscape, offering a comprehensive view of the system's energetics over the considered collective variables.

![Alt text da image](https://github.com/sulfierry/MolecularModeling/blob/main/FreeEnergyLandscape/fel_angle_distance.png)

- **X-Axis**: CV1, representing an angle, with the same labeling as in Figure 1 (left).
- **Y-Axis**: CV2, representing a distance, with the same labeling as in Figure 1 (right).
- **Visualization**: A heatmap showing the free energy levels across the CV1 and CV2 space. The color gradient represents different energy levels, with cooler colors indicating low-energy regions (minima) and warmer colors highlighting high-energy barriers. This visualization is crucial for identifying transition states, stable conformations, and understanding the molecular system's behavior under various conditions.

### Figure 5: 3D Free Energy Landscape

This figure provides a three-dimensional visualization of the free energy landscape, combining both collective variables, CV1 and CV2, along with the calculated free energy values to offer a dynamic perspective on the system's energetics.

![Alt Text](https://github.com/sulfierry/MolecularModeling/blob/main/FreeEnergyLandscape/energy_landscape_3D.gif)


- **X-Axis**: CV1, representing an angle, with the same labeling as in Figure 4.
- **Y-Axis**: CV2, representing a distance, with the same labeling as in Figure 4.
- **Z-Axis**: Shows the free energy in kJ/mol, calculated from the probability distribution of CV1 and CV2 using the Boltzmann equation.

**Visualization**: A 3D surface plot illustrates the variations in free energy across the space defined by CV1 and CV2. The color gradient enhances the visualization of energy barriers and minima, aiding in the identification of stable conformations and transition states. The 3D and animated visualizations of the free energy landscape extend the analysis beyond two dimensions, offering a richer and more nuanced understanding of the system's energetics. By exploring the landscape in three dimensions, researchers can better identify and analyze the regions of interest, such as low-energy conformations and high-energy transition states, which are vital for deciphering the molecular mechanisms underlying biological processes and material behaviors.

## Interpretation

Together, these figures provide a multi-faceted view of the molecular system's free energy landscape. Analyzing these visualizations helps in understanding how variations in critical structural parameters (angles and distances) influence the stability and dynamics of the system, which is vital for predicting reaction pathways, designing drugs, and engineering materials with desired properties.


# Comparative Analysis of Free Energy Calculation Methods

This document provides a detailed mathematical comparison of free energy calculation and interpolation methods between the custom Python script `freeEnergyLandscape.py` and the method of GROMACS tools `gmx_sham` and `gmx_wham`.

# GMX SHAM Mechanism Description

The `gmx_sham` tool is part of the GROMACS suite, designed to analyze free energy landscapes from simulation data. The core functionality revolves around creating multidimensional histograms from the provided data and calculating the free energy landscape. Here's a detailed breakdown of the process, including the improved explanation of "accumulated probability of all data points":

## Data Preparation

1. **Extremes Determination**: For each eigenvector, the code calculates minimum and maximum values across all data points, adding a small buffer to ensure all data is encompassed (`find_extremes` function).

2. **Normalization and Volume Correction**: Data undergo normalization and volume correction if necessary, based on user-specified dimensions. This step adjusts the probability calculations for spatial dimensions, correcting for volume effects that increase with distance in 2D and 3D spaces, as seen in the `correct_for_volume_effects` method.

## Histogram Creation and Free Energy Calculation

1. **Binning**: Data is sorted into multidimensional bins, each representing specific intervals of the eigenvectors (`binning_data` function). Points are allocated to bins based on their eigenvector values.

2. **Accumulated Probability Calculation**: The accumulated probability ($P_{acum}(bin)$) for each bin is calculated as the number of data points within the bin ($(N_{bin})$) divided by the total number of data points ($(N_{total})$), reflecting the proportion of occurrences for each bin relative to the entire dataset.

$$P_{acum}(bin) = \frac{N_{bin}}{N_{total}}$$

3. **Free Energy Calculation**: Free energy ($G$) for each bin is determined using the inverse Boltzmann relation, where the accumulated probability informs the energy level:

$$G(bin) = -kT \ln(P_{acum}(bin))$$

   Here, ($k$) is Boltzmann's constant, and ($T$) is the system's temperature. This step is executed within the `calculate_free_energy` method.

4. **Probability Normalization and Energy Adjustment**: Finally, the probability in each bin is normalized, and the free energy values are adjusted so the minimum free energy across the landscape is set to zero, facilitating easier interpretation and visualization.


# GMX WHAM Mechanism Description

The `gmx_wham` tool within the GROMACS suite utilizes the Weighted Histogram Analysis Method (WHAM) to derive free energy landscapes from multiple simulations. This sophisticated statistical approach combines data from various histograms to estimate the system's free energy landscape accurately. Additionally, `gmx_wham` employs linear interpolation for handling tabulated potentials. Below is a detailed explanation of WHAM's process, incorporating the clarification regarding interpolation methods:

## Data Preparation and Histogram Creation

1. **Histogram Generation**: Each simulation dataset contributes a histogram based on a specific reaction coordinate, representing the frequency distribution of system states across the sampled parameter space.

2. **Bias Correction**: To effectively explore specific regions of the parameter space, simulations often apply a bias. WHAM corrects these biases across all histograms, ensuring equitable contributions to the final analysis.

## Statistical Combination and Free Energy Calculation

WHAM utilizes a statistical method to combine histograms and derive the free energy landscape:

1. **Combining Histograms**: The method combines biased histograms using iteratively adjusted weights. Each histogram's weight reflects its contribution to the overall free energy calculation, considering the applied bias during simulation.

2. **Iterative Solution**: WHAM solves a set of self-consistent equations to find the weights that maximize the likelihood of the combined histogram data:

   - Let $N_i(j)$ be the number of counts in the $j^{th}$ bin of the $i^{th}$ histogram, with $n_i$ total observations and a biasing energy $U_i(j)$.
   - The unbiased probability distribution $P(j)$ for the $j^{th}$ bin is estimated by:

$$P(j) = \frac{\sum_{i} N_i(j)}{\sum_{i} n_i \exp\left[-\beta (U_i(j) - F_i)\right]}$$

Here, $\beta = 1/(k_BT)$, $k_B$ is Boltzmann's constant, $T$ is the temperature, and $F_i$ is the iteratively adjusted free energy of the $i^{th}$ simulation.

3. **Free Energy Landscape**: The free energy $G$ for each state is calculated from the probability distribution $P(j)$:

$$G(j) = -k_BT \ln(P(j))$$

   This reveals the energetically favorable states and barriers between them.

## Linear Interpolation for Tabulated Potentials

In addition to the statistical combination of histograms, `gmx_wham` uses linear interpolation within the `tabulated_pot` function for tabulated potentials. This method estimates potential energy values at intermediate distances, providing a continuous potential energy landscape.

However, this linear interpolation method is specifically applied to the scenario of dealing with tabulated potentials and does not directly influence the primary WHAM algorithm's statistical combination of histograms for free energy calculation. The WHAM methodology itself does not inherently use linear interpolation as part of its core algorithm for combining histograms or calculating the free energy landscape. Instead, WHAM relies on a statistical approach to optimally combine data from multiple biased simulations to reconstruct the unbiased free energy profile.


## Free Energy Landscape Calculation with Python Script

This document outlines the use of Kernel Density Estimation (KDE) and Boltzmann Inversion by the `freeEnergyLandscape.py` script to estimate the free energy landscape of a system from its collective variables.

### Kernel Density Estimation (KDE) Method

KDE smooths the distribution of collective variables to enhance the estimation of probability densities. It employs Scott's rule for bandwidth selection, ensuring an appropriately smooth density estimate for the data's scale.

$$
\text{Bandwidth (Scott's Rule)} = \sigma \cdot n^{-1/5}
$$

- $\sigma$ is the standard deviation of the sample.
- $n$ is the sample size.

KDE for a set of points is given by:

$$
\hat{f}(x) = \frac{1}{n \cdot h} \sum_{i=1}^{n} K\left( \frac{x - x_i}{h} \right)
$$

- $n$ is the number of points
- $K$ is the kernel function, typically Gaussian.
- $h$ is the bandwidth.

### Boltzmann Inversion

The free energy $G(x)$ of a state is calculated from the probability density $\hat{f}(x)$, obtained via KDE:

$$
G(x) = -k_BT \ln(\hat{f}(x))
$$

- $k_B$ is the Boltzmann constant.
- $T$ is the temperature.

### Probability Normalization

To ensure accurate energy calculation, the probability densities are normalized over the entire range of possible values, ensuring they sum to one across the configurational space:

$$
\int_{-\infty}^{\infty} \hat{f}(x) \, dx = 1
$$

### Interpolation

The smoothly interpolated probability density function provided by KDE enables calculating the free energy landscape over a continuous range of collective variable values, highlighting energetically favorable states and barriers between them.

## Conclusion

In conclusion, while `gmx_sham` and `gmx_wham` are powerful tools for free energy analysis within the GROMACS environment, `freeEnergyLandscape.py` offers a more user-friendly and customizable approach, particularly beneficial for those seeking immediate visualization and flexible data analysis options.
