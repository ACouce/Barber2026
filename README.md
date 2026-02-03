# Barber2026
Simulation model of mutation-biased adaptation:

The file mutT_scan_mu.R contains an example of the basic code used to simulate how strong mutation biases can drive mutator populations along divergent adaptive pathways across a range of population sizes in a haploid asexual system. All relevant parameters, variables, and processes are commented throughout the script. To run the model, simply execute the source file in any R interpreter (e.g., RStudio on Windows, macOS, or Linux; bash console in Linux) as follows:

source("mutT_scan_mu.R")

The code outputs several figures: one showing the genotype–fitness map used in all simulations, and a variable number of plots illustrating how the strength of mutation bias scales with population size across different values of the interrogated parameter (here, the basal mutation rate). One plot is generated per population size and mutation rate combination. Typical execution time is approximately 20–25 minutes on a 2.50 GHz Intel(R) Core(TM) i7-11850H CPU.

REFERENCE
Barber, J., & Couce, A. (2025). Mutation-biased adaptation is consequential even in large bacterial populations. bioRxiv. https://doi.org/10.1101/2025.06.16.655099
