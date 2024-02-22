This folder contains the following files:

1.revision_code_yieldgaps_clean.R : R script estimating the yield gap between organic and conventional farming due to differences in pesticides use;
Notes: Please be aware that this script requires access to farm-level FADN data for its execution, which under the EU GDPR regulation cannot 
be shared by any circumstances to third parties. Nevertheless, to maximize transparency we include the R-script used to calculate the yield gap. 
The R-script consists of two main functions: One function calculates the coefficients of the covariates influencing yields and a second one 
converts the coefficient to a yield gap. The script aligns with "Yield gaps between organic and conventional farming" under the method section.

2.mixed_effect_model.R : R script running the mixed effect model capturing the impact of natural pest control (NPC) on yield gaps;
Notes: Please ensure you have the necessary data (input_npc_yield_gap_data.xlsx) saved in the same directory as this R script. This will allow you
to execute the script without any problems and obtain the same results as presented in output_mixed_effect.xlsx.

3.input_npc_yield_gap_data.xlsx : Excel file containing the data on LF-NPC and yield crops inserted to run mixed_effect_model.R;

4.output_mixed_effect.xlsx : Excel file containing the output from mixed_effect_model.R;

5.results_NPC_income.xlsx: Excel file containg the data on region LF-NPC and income (baseline and future simulation).
Notes: Kindly note that the income has been estimated by the CAPRI modelling system. The model is open-source. 
See: https://www.capri-model.org/doku.php?id=capri:install.

6. variation_FADN_regions.xlsx: Excel file containing the LF-NPC and yield gap (per crop) variation within FADN regions.
The variation is captured  using various metrics: minimum, maximum, 25th and 75th percentiles, mean, median, and standard deviation.
For Landscape Feature Natural Pest Control (LF-NPC) data, statistics are calculated across all agricultural pixels within each region.
Yield gap data is presented per crop, with dedicated tabs for each individual crop.
