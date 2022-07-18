# pathogen-transmission-contact-network
Code for "Pathogen transmission modes determine contact network structure, altering other pathogen characteristics"

The full respository for networks used in this study can be found in the Animal Social Network Repository (ASNR) https://github.com/bansallab/asnr.
"Network_summary_master_file" is a summary of the networks and metrics used inlcuded in this study.

The folder "232_Network_Subsample" is a random sample of networks from each study included in our analysis, pulled from the ASNR with a max of 15 networks per study.

"GLMM_Sript.R" is the R script use to run 1000 GLMMs on our data, pulling a different subsample of networks from the ASNR each time. It will generate the "Model_List_1000_runs.rds" needed in the "pull_from_posterior.R" script, as well as the pvalues used to infer signficant differences in network structure among networks representative of each transmission mode.

"pull_from_posterior.R" is an R script that will use the effect sizes generated in the 1000 GLMMs to generate predictions for each network metric across transmission modes. (Used to generate Figure 1 in the manuscript). It will generate "Model_Predictions.csv"

"FindTc_3scenarios.py" will calculat the critical transmissibility (Tc) for each network in the 232_Network_subsample folder. It will generate the Tc values using a disease simulation, as well as through two control scenarios that control the effects of homogeneous and heterogeneous degree. It will generate the "Tc_3Scenarios.csv file", and is used to generate Figure 2 in the manuscript. "important_functions.py" are functions used in the "FindTc_3scenarios.py" script.
