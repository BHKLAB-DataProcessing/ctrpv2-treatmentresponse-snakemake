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

###############################################################################
# Load INPUT
###############################################################################
sampleMetadata <- data.table::fread(INPUT$sampleMetadata)
treatmentMetadata <- data.table::fread(INPUT$treatmentMetadata)
preprocessed_raw <- data.table::fread(INPUT$preprocessed_raw)
preprocessed_profiles <- data.table::fread(INPUT$preprocessed_profiles)

###############################################################################
# Main Script
###############################################################################

# merge preprocessed raw data with sample and treatment metadata
# raw$sampleid == sampleMetadata$CCLE.sampleid
# raw$treatmentid == treatmentMetadata$CCLE.treatmentid

message("Merging preprocessed raw data with sample and treatment metadata...")
preprocessed_raw <- merge(
    preprocessed_raw,
    sampleMetadata,
    by.x = "sampleid",
    by.y = "CCLE.sampleid",
    all.x = TRUE
)
preprocessed_raw <- merge(
    preprocessed_raw,
    treatmentMetadata,
    by.x = "treatmentid",
    by.y = "CCLE.treatmentid",
    all.x = TRUE
)
# remove any duplicate rows
preprocessed_raw <- preprocessed_raw[!duplicated(preprocessed_raw), ]

# CONSTRUCT THE TREDataMapper OBJECT
message("Constructing the treatment response experiment object...")
tdm <- CoreGx::TREDataMapper(rawdata = preprocessed_raw)

# mapped columns are the columns in treatmentMetadata that are not in the
# id_columns
CoreGx::rowDataMap(tdm) <- list(
    id_columns = c("treatmentid", "dose"),
    mapped_columns = intersect(
        names(treatmentMetadata),
        names(preprocessed_raw)
    )
)

CoreGx::colDataMap(tdm) <- list(
    id_columns = c("sampleid"),
    mapped_columns = intersect(
        names(sampleMetadata),
        names(preprocessed_raw)
    )
)

CoreGx::assayMap(tdm) <- list(
    sensitivity = list(
        id_columns = c("treatmentid", "sampleid", "dose"),
        mapped_columns = c("viability")
    )
)

(ccle_tre <- CoreGx::metaConstruct(tdm))
ccle_tre$sensitivity


######################################################
# Compute on the sensitivity assay
message("Computing the profiles assay...")
tre_fit <- ccle_tre |>
    CoreGx::endoaggregate(
        {
            # the entire code block is evaluated for each group in our group by
            # 1. fit a log logistic curve over the dose range
            fit <- PharmacoGx::logLogisticRegression(
                dose,
                viability,
                viability_as_pct = TRUE
            )
            # 2. compute curve summary metrics
            ic50 <- PharmacoGx::computeIC50(dose, Hill_fit = fit)
            aac <- PharmacoGx::computeAUC(dose, Hill_fit = fit)
            # 3. assemble the results into a list, each item will become a
            #   column in the target assay.
            list(
                HS = fit[["HS"]],
                E_inf = fit[["E_inf"]],
                EC50 = fit[["EC50"]],
                Rsq = as.numeric(unlist(attributes(fit))),
                aac_recomputed = aac,
                ic50_recomputed = ic50
            )
        },
        assay = "sensitivity",
        target = "recomputed_profiles",
        enlist = FALSE, # this option enables the use of a code block for aggregation
        by = c("treatmentid", "sampleid"),
        nthread = THREADS # parallelize over multiple cores to speed up the computation
    )

###################################################################
# Add the published profiles
message("Adding the published profiles...")
CoreGx::assay(tre_fit, "published_profiles") <- preprocessed_profiles
######################################################
# Add metadata
CoreGx::metadata(tre_fit) <- list(
    data_source = snakemake@config$treatmentResponse,
    annotation = "treatmentResponse",
    date = Sys.Date(),
    sessionInfo = capture.output(sessionInfo())
)


show(tre_fit)
###############################################################################
# Save OUTPUT
###############################################################################
message("Saving the treatment response experiment object...")
saveRDS(
    tre_fit,
    file = OUTPUT$tre
)
