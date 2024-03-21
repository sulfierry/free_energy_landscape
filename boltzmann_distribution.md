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
