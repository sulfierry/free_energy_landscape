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

## Usage

Ensure your data file is in the correct two-column format. Run the script with the path to your data files and optional arguments as needed:

```bash
python freeEnergyLandscape.py path/to/cv1_data.txt path/to/cv2_data.txt
```
```bash
Optional Arguments:

   --temperature           [int]       Simulation temperature in Kelvin (default: 300K)
   --kb                    [float]     Boltzmann constant in kJ/(mol·K) (default: 8.314e-3)
   --energy                [int]       Energy, single value (default: None)
   --bins_energy_histogram [int]       Bins for energy histogram (default: 100)
   --kde_bandwidth         [float]     Bandwidth for kernel density estimation (default: None)
   --names                 [str] [str] Names for the collective variables (default: CV1, CV2)
   --gif_angles            [int]       Angles for 3D GIF rotation (default: 10)
   --gif_elevation         [int]       Elevation angle for the 3D GIF (default: 10)
   --gif_duration          [float]     Duration per frame in the GIF in seconds (default: 0.1)
```
```bash
python freeEnergyLandscape.py path/to/cv1_data.txt path/to/cv2_data.txt  --names CV1_angle CV2_distance --energy 3
```

## Theoretical Background

The free energy landscape is a conceptual and computational tool used to understand the energetics and dynamics of molecular systems. It is particularly useful in the study of complex processes such as protein folding, chemical reactions, and phase transitions.

### Collective Variables (CVs)


Collective Variables (CVs) are a set of coordinates that describe the macroscopic state of a system. They are used to reduce the complexity of molecular systems by focusing on the relevant degrees of freedom. Examples include the distance between two atoms, angles, dihedrals, and more complex descriptors. 

The use of distance and angle as CVs in this study serves merely as an example to illustrate the tool's capabilities. However, it's important to note that the input can be any file containing two columns: the first for frames and the second for the value of the collective variable. This flexibility allows the tool to be applicable to a wide range of studies involving different types of collective variables

CVs in biomolecular systems are mathematically represented as functions of the atomic coordinates. These functions are designed to capture the essential features of the system's configuration that are relevant to its macroscopic properties or behaviors. Below are examples of commonly used CVs and their mathematical formulations:

1. **Distance**: The distance $d$ between two atoms $i$ and $j$ with positions $\vec{r}_i$ and $\vec{r}_j$ can be calculated using the Euclidean distance formula:

   First, determine the position vectors of atoms $i$ and $j$:
   $$\vec{r}_i = (x_i, y_i, z_i)$$
   $$\vec{r}_j = (x_j, y_j, z_j)$$

   Then, calculate the distance $d$ as:
   $$d = |\vec{r}_i - \vec{r}_j| = \sqrt{(x_i - x_j)^2 + (y_i - y_j)^2 + (z_i - z_j)^2}$$

   - $\vec{r}_i$ and $\vec{r}_j$ are the position vectors of atoms $i$ and $j$, respectively.
   - $(x_i, y_i, z_i)$ and $(x_j, y_j, z_j)$ denote the Cartesian coordinates of atoms $i$ and $j$.
   - $|\vec{r}_i - \vec{r}_j|$ represents the magnitude of the vector difference between $\vec{r}_i$ and $\vec{r}_j$, giving the direct distance between the two atoms.

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

These representations allow for a simplified description of the system's state, facilitating the study of its behavior and properties through the manipulation of a reduced set of variables rather than the full set of atomic coordinates.

### Kernel Density Estimation (KDE)

This statistical method is crucial for estimating the probability density function (PDF) of a dataset without assuming any predefined distribution shape. In the context of molecular dynamics and simulations, it allows us to visualize and analyze the distribution of free energy across different states or configurations defined by collective variables (CVs).

The `gaussian_kde` leverages a Gaussian (normal) distribution, placing it at each point in the dataset and summing these distributions to approximate the overall data's PDF. This technique is adept at capturing the underlying structure of the data, providing a smooth, continuous representation of the free energy landscape. The smoothness of the KDE is controlled by the bandwidth parameter, which determines the width of the Gaussian kernels used.

Mathematically, the density estimation at a point $`x`$ is calculated as follows:

$$\ f(x) = \frac{1}{n \cdot h} \sum_{i=1}^{n} K(u) \$$   

where:
- $K(u)$ represents the Gaussian kernel function

$$ K(u) = \frac{1}{\sqrt{2\pi}} e^{-\frac{1}{2}u^2} $$

- $u$ is the normalized distance

$$ u = \frac{x - x_i}{h} $$

- $n$ is the number of data points
- $h$ is the bandwidth
- $x$ is where the density is being estimated
- $x_i$ are the data points
  
note that:
- $u^2$ is the square of this normalized distance, which serves to weight the contribution of each data point $x_i$ to the density estimate at $x$, based on how far $x_i$ is from $x$, adjusted by the bandwidth $h$.
- $e^{-\frac{1}{2}u^2}$ decreases rapidly as $u$ increases, meaning that points further away from $x$ will have less influence on the density estimate at $x$.
- $\frac{1}{\sqrt{2\pi}}$ is a normalization term that ensures the Gaussian kernel function integrates to $1$, keeping it as a valid probability distribution.

In the context of KDE the variable $u$ is utilized within the Gaussian kernel function $K$, which is a probability density function. Thus, $u$ is crucial in determining how the distance between $x$ and the data points $x_i$ affects the density estimate at $x$, with the bandwidth $h$ controlling the sensitivity of this influence.

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
