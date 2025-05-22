## ------------------- Parse Snakemake Object ------------------- ##
# Check if the "snakemake" object exists
# This snippet is run at the beginning of a snakemake run to setup the env
# Helps to load the workspace if the script is run independently or debugging
if (exists("snakemake")) {
    INPUT <- snakemake@input
    OUTPUT <- snakemake@output
    WILDCARDS <- snakemake@wildcards
    THREADS <- snakemake@threads

    # setup logger if log file is provided
    if (length(snakemake@log) > 0)
        sink(
            file = snakemake@log[[1]],
            append = FALSE,
            type = c("output", "message"),
            split = TRUE
        )

    # Assuming that this script is named after the rule
    # Saves the workspace to "resources/"buildTreatmentResponseExperiment"
    file.path("resources", paste0(snakemake@rule, ".RData")) |>
        save.image()
} else {
    # If the snakemake object does not exist, load the workspace
    file.path("resources", "buildTreatmentResponseExperiment.RData") |>
        load()
}
####################################################################################
# 0.2 read in metadata
# --------------------
rawdt <- data.table::fread(
    INPUT$preprocessed_raw,
    header = TRUE,
    sep = ",",
    fill = TRUE
)

print(rawdt)


print(unique(rawdt$tech_rep))