
# input data
variant_files <- c(
  "/mnt/data2/ngs/autismo_test__08_11_2011/family_variant_db_casefam.txt"
)
family_names <- c("grupo_afectos")
outdir <- "../autismo_results_08_11_2011"
dominant_model <- FALSE

global.disease_name <- "Idiopathic Autism"
global.disease_phenotype <- "Recessive"
global.study_description <- "This is a Case/Control study with 29 cases and 40 healthy controls"
global.phenotype_selection <- "variant must be homozygous in at least 30% of cases and they cannot be homozygous in any controls"
global.effect_selection <- "non synonymous, stop codon gain/loss, splicing site and ncRNA"

########################### EXECUTION

source("disease_report.R")
source("interactome_analysis.R")
source("do_report.R")