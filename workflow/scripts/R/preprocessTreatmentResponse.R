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

zipDir <- tools::file_path_sans_ext(INPUT$tr)
dir.create(zipDir, recursive = TRUE, showWarnings = FALSE)
print(paste0("Unzipping ", INPUT$tr, " into ", zipDir))
unzip(INPUT$tr, exdir = zipDir)

files <- paste0(zipDir, "/", list.files(zipDir))


# Get the treatment metadata
# --------------------------
treatmentMetadata <- 
    files[grepl(pattern = "*per_compound*", x = files)] |>
    data.table::fread() 
treatmentMetadata[, treatmentid := cpd_name]

# Get the sample metadata
# -----------------------
sampleMetadata <- 
    files[grepl(pattern = "*per_cell_line*", x = files)] |>
    data.table::fread()

sampleMetadata <- sampleMetadata[, .(
        sampleid = ccl_name, 
        tissueid = ccle_primary_site, 
        master_ccl_id = as.character(master_ccl_id), 
        sample_availability = ccl_availability, 
        ccle_primary_hist = ccle_primary_hist)]


# Read in sensitivity data
# ------------------------
sensitivityRaw <- 
    files[grepl(pattern = "*per_cpd_post_qc*", x = files)] |>
    data.table::fread()


# code taken from https://github.com/BHKLAB-DataProcessing/PSet_CTRPv2-snakemake/blob/master/Oldscripts/downloadSensData.R
ctrp.sensitivityInfo <-  read.delim(files[grepl(pattern = "*meta.per_experiment*", x = files)])

repExps <- unique(ctrp.sensitivityInfo$experiment_id[duplicated(ctrp.sensitivityInfo$experiment_id)])

## Just collapsing the dates, everything else is identical.
for (exp in repExps) {
  myx <- ctrp.sensitivityInfo$experiment_id == exp
  duplicates <- duplicated(ctrp.sensitivityInfo$experiment_id) & myx
  first <- myx & !duplicates
  # print(ctrp.sensitivityInfo[myx,])
  # browser()
  ctrp.sensitivityInfo[first, ] <- apply(ctrp.sensitivityInfo[myx, ], 2, function(x) paste(unique(x), collapse = "//"))

  ctrp.sensitivityInfo <- ctrp.sensitivityInfo[!duplicates, ]
}


sensitivityInfo <- data.table::as.data.table(ctrp.sensitivityInfo)

sensitivityInfo <- merge(
    sensitivityInfo,
    sampleMetadata,
    by = "master_ccl_id",
    all.x = TRUE
)

# convert experiment_id to integer
sensitivityInfo[, experiment_id := as.integer(experiment_id)]

sensitivityRaw <- merge(
    sensitivityRaw,
    sensitivityInfo,
    by = "experiment_id"
)

sensitivityRaw <- merge(
    sensitivityRaw,
    treatmentMetadata,
    by = "master_cpd_id"
)

# TODO:: KEEP ALL THE OTHER COLUMNS!!! no need to remove them here 
rawdt <- sensitivityRaw[, .(
    sampleid, culture_media, treatmentid, experiment_id, 
    dose = cpd_conc_umol, viability =  100*cpd_avg_pv
)]

# Add technical replicate column 
# rawdt[, .N, by = .(sampleid, treatmentid, culture_media, dose)][order(N)]

rawdt[,
    tech_rep := seq_len(.N),
    by=.(sampleid, treatmentid, culture_media, dose)
]


data.table::fwrite(
    x = rawdt,
    file = OUTPUT$preprocessed_raw,
    sep = ",",
    col.names = TRUE,
    row.names = FALSE,
    quote = TRUE 
)