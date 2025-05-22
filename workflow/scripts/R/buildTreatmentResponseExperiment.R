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

# -- Create a TREDataMapper
print("Creating a TREDataMapper")
treDataMapper <- CoreGx::TREDataMapper(rawdata = rawdt)

CoreGx::rowDataMap(treDataMapper) <- list(
    id_columns = c("treatmentid", "dose", "tech_rep"),
    mapped_columns = c()
)

CoreGx::colDataMap(treDataMapper) <- list(
    id_columns = c("sampleid", "culture_media"),
    mapped_columns = c()
)

CoreGx::assayMap(treDataMapper) <- list(
    raw = list(
        id_columns = c("treatmentid", "dose", "tech_rep", "sampleid", "culture_media"),
        mapped_columns = c("viability")
    )
)



print("Running CoreGx::metaConstruct")
(tre <- CoreGx::metaConstruct(treDataMapper))



print("Endoaggregating to create mono_viability Assay")
tre_sens <- tre |>
    CoreGx::endoaggregate(
        assay="raw",
        target="sensitivity",  # create a new assay named mono_viability
        mean_viability=pmin(100, mean(viability)), # pmin takes the minimum of two vectors element-wise # idk why the old script did this step...
        by=c("treatmentid", "dose", "sampleid"),
        nthread=THREADS
    )

saveRDS(
    tre_sens,
    file = OUTPUT$tre,
)
