This folder contains the scipts used for the bootstrap analysis of the Affect Valuation Theory data. 


sem ram apply AVT3g bootstrap.R does the actual resampling and model fitting. Each round of resampling 
uses seeds taken from the random_seed_bootstrap.txt file. The seeds that were used were then dumped in
the used_seeds.txt file. The results from each model fit to the resampled data are found in the Results
directory. 

An analysis of the results are done in the avtboot_posthoc_analysis.R script.

