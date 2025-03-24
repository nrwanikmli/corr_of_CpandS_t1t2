# corr_of_CpandS_t1t2
Scripts for generating graphs and getting result to evaluate the prediction potential of the first and second general temperature indices for thermodynamic properties(entropy & heat capacity) of Benzenoid Hydrocarbons

Each file's brief description:

-octave-
1) combined_corr_21b2015.m: This script analyzes the predictive potential of general temperature indices (first & second general temperature indices) for entropy and heat capacity by finding optimal correlation coefficients and plotting correlation curves. This script generates four correlation curve plots, two for each thermodynamic properties.

2) good_alpha_21b2015.m: This script generates 4 plots. It calculates the best value of alpha and highlights the "good" alpha-intervals on the graphs. The intervals and correlation values are also printed in the console

3) scatter_21b2015.m: This script generates 4 scatter plots, each plot includes the regression line between the data points. The script computes the correlation coefficient, finds the best value of alpha that maximizes the correlation, and prints the results for the slope, intercept, and standard error.

4) golden_section_search_max.m: This script must be placed in the working directory as it is referenced by the other three scripts. It is the Octave 9.2 implementation of the Golden Section Search algorithm, manually translated from Python code. The function uses the Golden Section Search method to find the interval containing the maximum value of a function within a specified range, with a defined tolerance for accuracy

-Rstudio-
5) multcorr_t1.R: Generates visualization of the relationship between a first general temperature index and the two thermodynamic properties (S & Cp) and matrix plot visualising the variable distribution and their bivariate relationship
6) multcorr_t2.R: Generates visualization of the relationship between a second general temperature index and the two thermodynamic properties (S & Cp) and matrix plot visualising the variable distribution and their bivariate relationship
