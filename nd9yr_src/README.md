# Data & code repository for study on metacognitive insight around COVID-19 knowledge and its relation to health protective behaviours

This repository contains all data and code to run all main analyses and generate the figures in the main text.

All the data that was collected in this study is included in the `data` folder; in particular the following files are included here exactly in the format they were provided by the agency that collected them (YouGOv)

- `UniOfEssex_Results_210414.xls`: averages (weighted) for each question
- `UoEssex_Results_210414_Client.csv`: detailed data
- `UoEssex_Results_210414_Codebook.xlsx`: codebook providing explanations for all codes used in the detailed data file.

All other files in the `data` folder do not contain actual data but only provide the labels and recoded categories that are used in the analyses and plotting.

All these analyses use MCMC sampling (either in Jags or Stan) and the results consists of sampels from the posterior distribution, which occupy quite a lot of memory space and are therefore not provided here. 

To reproduce analyses and results, run the scripts as follow in R:


    # Bayesian multilevel estiamtion of metacognitive efficiency:
    source('prepare_data.R')
    source('run_jags_model.R')      # (very slow)
    
    # Ordinal regressions
    source('recode_variables_for_regression.R')
    source('run_ordinal_regressions.R') # also slow-ish
    
    # Figures 
    source('make_fig1_main_text.R') # to make fig 1 as in main text
    source('generate_plots_regression_results.R') # including fig. 2 and 3 in main text

    # Other/Misc
    source('ordinal_pseudo_Rsquared.R') # as the name say, compute pseudo R-squared for ordinal models
    source('supplemental_analysis_confidence_bias.R') # Supplemental analysis that     compare un-signed confidence across COVID-19 and science items

See the Supplemental Material (`SupplementalInformation.pdf`) to see the original computing environment in which these analyses were run.


If you have questions about the code and data, please do get in touch at: matteo.lisi [at] rhul.ac.uk

