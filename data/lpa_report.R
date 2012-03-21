########################### LPA 2 replicates

# input data
variant_files <- c(
  "../data/lymphoma_2_replicates/family_variant_db_lpa1.txt",
  "../data/lymphoma_2_replicates/family_variant_db_lpa2.txt",
  "../data/lymphoma_2_replicates/family_variant_db_lpa3.txt",
  "../data/lymphoma_2_replicates/family_variant_db_lpa4.txt",
  "../data/lymphoma_2_replicates/family_variant_db_lpa5.txt"
)

family_names <- c("LPA1","LPA2","LPA3","LPA4","LPA5")
outdir <- "../lpa_results_09_11_2011"
dominant_model <- TRUE

global.disease_name <- "Acute Promyelocytic Leukemia"
global.disease_phenotype <- "Dominant"
global.study_description <- "This is a Case/Control study with APL families and 40 healthy controls"
global.phenotype_selection <- "variant must be in all family cases (heterozygous or homozygous) and they cannot be present in healthy family members and healthy controls"
global.effect_selection <- "non synonymous, stop codon gain/loss, splicing site and ncRNA"

########################### EXECUTION

source("disease_report.R")
source("interactome_analysis.R")
source("do_report.R")