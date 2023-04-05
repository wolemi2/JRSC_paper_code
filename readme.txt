#Title: Multivariate sensitivity analysis for a large-scale climate impact and adaptation model
Oyebamiji, O., Nemeth, C., Harrison, P., Dunford, R., & Cojocaru, G. (2022). Multivariate sensitivity analysis for a large-scale climate impact and adaptation model. arXiv preprint arXiv:2201.09681.

#This main folder contains R scripts and data for producing results in this paper
(i) The folder named as "Rcode" contains all the R codes for fitting and plotting results in the paper.
(ii) The folder named as "data_in" contains all the data for the analysis.
(iii) The folder named as "data_out" contains all the model outputs and results.
(iv) 


##Some fixed parameters
p <- 24 #number of input variables
m <- 18 #no of output 
N <- 2000 # no of original samples
T <- 3 ## time slices
sc <- 4 # no of scenarios considered
TS <- 15 # 
ns <- 600 #number of subsamples
nr <- r <- 6 #number of region

#Data used for the analysis
1. Ysen0.txt 15*2000*18*6 array of output data. Note: 4th dimension [6 <- 5 regional averages + 1 EU-wide averages].
2. Xsen0.txt 15*2000*24 array of input data.
3. sce_range.txt  # 17*24 matrix of scenario ranges.
4. scenario.txt # 15*2000*18*4 # scenario-based data.
5. lonlat.txt - longitude and latitude data for plotting regional and spatial maps.
6. IAP_regions.mat - IAP regional data.


#R routines
1. config.r - Routine for setting various configuration parameters eg number of iterations and thin for MCMC
2. main.r - The main R codes calling other routines. The main path/directory has to be set from here.
3. EU_wide_MSGP.r - The main code for fitting the EU_wide MSGP emulator
4. regional_MSGP.r - This code is for fitting the regional MSGP emulators.
5. sparsity_MSGP.r - This code is for fitting the different levels of sparsity to the EU_wide MSGP emulator.
6. fit_matern.r - This code fit the Matern covariance using the Wendland function.
7. regional_sensitivity.r - This code performs the regional-based sensitivity indices.
8. scenario_sensitivity.r - This code performs the scenario-based sensitivity indices.
9. sensitivity_func.r - Contains miscellaneous functions for sensitivity analysis.
10. mcmc_func.r - This function is called by other routines (regional_MSGP and sparsity_MSGP) to fit the MSGP models and perform the crossvalidation.
11. lik.r - This file contains the likelihood functions.
12. metropolis.r - This file contains the metropolis function.
13. functions_libraries.r - This routine contain useful functions and R libraries.
14. pred.r - This file contains the predictive function for the MSGP models.
15. plot.r - The R code for plotting Figures 1, 8, 9, 10, 11, 12, 13, 14, and 15.
16. diagnostic.r - Routine for plotting Diagnostics in Figures 6, 7, 16 and 17.
17. test_cases.r - for fitting three test cases and plotting Figures 2, 3 and 4.
18. dicekrig_mult.r - This code fit Standard GP using DiceKriging library and multivariate regression models to sample of EU_wide data.


#Examples on how to fit the models
#######to fit EU_wide MSGP models and generate plots
system.time(source("../JRSC_paper_code/Rcode/EU_wide_MSGP.r",echo=TRUE))
####to fit regional models
system.time(source("../JRSC_paper_code/Rcode/regional_MSGP.r",echo=TRUE))

#The following R libraries are required
gplots
Rcpp
sensitivity
data.table
parallel
doParallel
RcppArmadillo
adaptMCMC
abind
RColorBrewer
ggplot2
TSP
coda
Matrix

Dr. O.K. Oyebamiji
Flood and Water Management
HR Wallingford
Howbery Park, Wallingford
Oxfordshire, OX10 8FZ
United Kingdom
E-mails: wolemi2@yahoo.com; o.oyebamiji@hrwallingford.com
 