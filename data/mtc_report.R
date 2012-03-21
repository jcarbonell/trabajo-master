########################### MTC

# input data
variant_files <- c(
  "../data/mtc_bfast_gatk_filtered_test_5_02_09_2011/family_variant_db_casefam_mgp41.txt",
  "../data/mtc_bfast_gatk_filtered_test_5_02_09_2011/family_variant_db_casefam_mgp43.txt",
  "../data/mtc_bfast_gatk_filtered_test_5_02_09_2011/family_variant_db_casefam_mgp44.txt"
  )

family_names <- c("MGP41","MGP43","MGP44")
outdir <- "../mtc_results_5_02_09_2011"
dominant_model <- TRUE

global.disease_name <- "Medullary thyroid cancer"
global.disease_phenotype <- "Dominant"
global.study_description <- "This is a Case/Control study with 3 MTC families and 10 healthy controls"
global.phenotype_selection <- "variant must be in all family cases (heterozygous or homozygous) and they cannot be present in healthy family members and healthy controls"
global.effect_selection <- "non synonymous, stop codon gain/loss, splicing site and ncRNA"

########################### EXECUTION

source("disease_report.R")
source("interactome_analysis.R")
source("do_report.R")