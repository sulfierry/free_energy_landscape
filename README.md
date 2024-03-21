[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.10689690.svg)](https://doi.org/10.5281/zenodo.10689690)


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
free_energy_landscape path/to/cv1_data.txt path/to/cv2_data.txt --names "CV1 (Angle)" "CV2 (Distance)" --energy 3.0 --discretize 1.0
```
In this previous bash example we will display in the plots all points that exhibit energy ($<=3 KJ/mol$) for every $1.0 KJ/mol$. And the axis title of each plot will be "CV1 (Angle)" and "CV2 (Distance)" respectively.

## Theoretical Background

The free energy landscape is a conceptual and computational tool used to understand the energetics and dynamics of molecular systems. It is particularly useful in the study of complex processes such as protein folding, chemical reactions, and phase transitions.

### Collective Variables (CVs)

Collective Variables (CVs) are a set of coordinates that describe the macroscopic state of a system. They are used to reduce the complexity of molecular systems by focusing on the relevant degrees of freedom. Examples include the analysis of principal components (PCA), the distance between two atoms, angles, dihedrals, and more complex descriptors. Below are examples of commonly used CVs and their mathematical formulations.

Geometrically, PCA seeks the line (in the case of a single principal component) or hyperplane (in the case of multiple principal components) that best captures the distribution of the data in a high-dimensional space. The first principal component is the direction along which the projection of the data points has the largest variance. This direction corresponds to the eigenvector of the covariance matrix associated with the largest eigenvalue. The second principal component is orthogonal to the first and represents the direction of the next highest variance, corresponding to the eigenvector associated with the second largest eigenvalue, and so on.

1. **Principal Component Analysis (PCA)**: Can be expressed as:
$$PC_k = \sum_{i=1}^{p} a_{ki}X_i$$
where $PC_k$ is the $k_{th}$ principal component, $X_i$ are the original variables, and $a_{ki}$ are the coefficients (loadings) for the $k_{th}$ principal component, given by the $k_{th}$ eigenvector of the covariance matrix $\mathbf{C}$. The transformation can be represented in matrix form as:
$$\mathbf{PC} = \mathbf{X}\mathbf{A}$$
where $\mathbf{PC}$ is the matrix of principal components, $\mathbf{X}$ is the original data matrix (centered or standardized, if necessary), and $\mathbf{A}$ is the matrix whose columns are the eigenvectors of $\mathbf{C}$.

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

## Visualizations Generated

The `FreeEnergyLandscape` class generates three key visualizations to aid in the analysis of the molecular system's free energy landscape. Each figure provides unique insights into the system's thermodynamic and kinetic behavior.

### Figure 1:  Normalized Free Energy Profile Comparison (CV1 and CV2)

This updated figure now integrates the normalized free energy profiles for both CV1 (Angle) and CV2 (Distance) in a single, unified visualization. This enhancement allows for a direct comparison between the two collective variables, illuminating their unique energy landscapes and pinpointing low-energy regions across specified energy threshold ranges.

![Alt text da image](https://github.com/sulfierry/free_energy_landscape/blob/main/1_Combined_Free_Energy_Profile_Normalized.png?raw=true)

- **X-Axis**: Displays both CV1 and CV2, with CV1 representing an angle and CV2 denoting distance in Ångströms, facilitating a comprehensive view of the system's energetics.
- **Y-Axis**: Represents the normalized free energy, derived from the probability distributions of CV1 and CV2 using the Boltzmann equation.
- **Visualization**: A comparative plot of free energy versus both CV1 and CV2, highlighting key energy barriers and minima corresponding to the molecular system's different conformational states.

### Figure 2: Histograms of Normalized CV1 and CV2 Values Side by Side

![Alt text da image](https://github.com/sulfierry/free_energy_landscape/blob/main/2_histograms_normalized_side_by_side.png?raw=true)

- **Visualization**: Two histograms placed side by side, one for CV1 (Angle) and the other for CV2 (Distance), each showing the normalized value distribution. This format allows for an easy comparison of the two variables' behaviors and predominant states within the simulation.


### Figure 3: CV Value vs. Frame for CV1 and CV2

This figure innovatively plots the values of CV1 (Angle) and CV2 (Distance) against simulation frames in the same image. It offers insights into how the values of these collective variables evolve over time, shedding light on the dynamic aspects of the system's behavior.

![Alt text da image](https://github.com/sulfierry/free_energy_landscape/blob/main/3_cv_by_frame_combined_normalized.png?raw=true)


- **X-Axis**: Represents the simulation frames, providing a temporal dimension to the analysis.
- **Y-Axis**: Displays values of CV1 and CV2, allowing for a direct observation of how each collective variable changes throughout the simulation.
- **Visualization**: A combined plot that tracks the fluctuation of both CV1 and CV2 over the course of the simulation, highlighting trends, transitions, and stability regions within the molecular system.

### Figure 4: 2D Free Energy Landscape (CV1 vs. CV2)

The third figure combines both CV1 and CV2 to produce a two-dimensional free energy landscape, offering a comprehensive view of the system's energetics over the considered collective variables.

![Alt text da image](https://github.com/sulfierry/free_energy_landscape/blob/main/4_Free_energy_landscape.png?raw=true)

- **X-Axis**: CV1, representing an angle, with the same labeling as in Figure 1 (left).
- **Y-Axis**: CV2, representing a distance, with the same labeling as in Figure 1 (right).
- **Visualization**: A heatmap showing the free energy levels across the CV1 and CV2 space. The color gradient represents different energy levels, with cooler colors indicating low-energy regions (minima) and warmer colors highlighting high-energy barriers. This visualization is crucial for identifying transition states, stable conformations, and understanding the molecular system's behavior under various conditions.

### Figure 5: 3D Free Energy Landscape

This figure provides a three-dimensional visualization of the free energy landscape, combining both collective variables, CV1 and CV2, along with the calculated free energy values to offer a dynamic perspective on the system's energetics.

![Alt Text](https://github.com/sulfierry/free_energy_landscape/blob/main/5_3D_landscape.png?raw=true)


- **X-Axis**: CV1, representing an angle, with the same labeling as in Figure 4.
- **Y-Axis**: CV2, representing a distance, with the same labeling as in Figure 4.
- **Z-Axis**: Shows the free energy in kJ/mol, calculated from the probability distribution of CV1 and CV2 using the Boltzmann equation.

**Visualization**: A 3D surface plot illustrates the variations in free energy across the space defined by CV1 and CV2. The color gradient enhances the visualization of energy barriers and minima, aiding in the identification of stable conformations and transition states. The 3D and animated visualizations of the free energy landscape extend the analysis beyond two dimensions, offering a richer and more nuanced understanding of the system's energetics. By exploring the landscape in three dimensions, researchers can better identify and analyze the regions of interest, such as low-energy conformations and high-energy transition states, which are vital for deciphering the molecular mechanisms underlying biological processes and material behaviors.

## Interpretation

Together, these figures provide a multi-faceted view of the molecular system's free energy landscape. Analyzing these visualizations helps in understanding how variations in critical structural parameters (angles and distances) influence the stability and dynamics of the system, which is vital for predicting reaction pathways, designing drugs, and engineering materials with desired properties.
