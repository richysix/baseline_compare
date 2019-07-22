#!/usr/bin/env Rscript
library('optparse')

option_list <- list(
  make_option("--species", type="character", default='Mus_musculus',
              help="Ensembl version number [default %default]" ),
  make_option("--assembly", type="character", default='GRCm38',
              help="Ensembl version number [default %default]" ),
  make_option("--ensembl_version", type="integer", default=88,
              help="Ensembl version number [default %default]" ),
  make_option("--metadata", type="character", default=NULL,
              help="YAML file of metadata [default %default]" ),
  make_option("--blacklist_file", type="character", default=NULL,
              help="Name of a file of Ensembl gene ids to remove [default %default]" ),
  make_option("--include_norm_counts", type="logical", action="store_true",
              default=FALSE,
              help="Include normalised counts [default %default]" ),
  make_option(c("-v", "--verbose"), action="store_true", default=FALSE,
              help="Print extra output [default]")
)

cmd_line_args <- parse_args(
  OptionParser(
    option_list=option_list, prog = 'create_baseline_rds.R',
    usage = "Usage: %prog [options] counts_file samples_file output_file_name\n",
    description = paste("Create SummarizedExperiment object for baseline data",
                   "from a count and sample file")
  ),
  positional_arguments = 3
)

library('SummarizedExperiment')
library('yaml')

# load count data
count_file <- cmd_line_args$args[1]
count_data <- read.delim(count_file, check.names = FALSE)
rownames(count_data) <- count_data[['Gene ID']]

# remove blacklist genes if supplied
if (!is.null(cmd_line_args$options[['blacklist_file']])) {
    blacklist_genes <- as.character(read.table(cmd_line_args$options[['blacklist_file']])[,1])
    to_remove <- intersect(rownames(count_data), blacklist_genes)
    to_keep <- setdiff(rownames(count_data), blacklist_genes)
    blacklist_only <- setdiff(blacklist_genes, rownames(count_data))
    
    # subset count data to genes to keep
    count_data <- count_data[ to_keep, ]
    # log which genes were removed to STDERR
    cat(sprintf('These genes were removed: %s',
                paste0(to_remove, collapse = ",")), sep='', file=stderr())
    cat(sprintf('These genes were in the blacklist, but not the count data: %s',
                paste0(blacklist_only, collapse = ",")), sep='', file=stderr())
}

# load sample data
samples_file <- cmd_line_args$args[2]
sample_info <- read.delim(samples_file, row.names = 1)

# set levels of stage
sample_info$stage <- factor(sample_info$stage,
                            levels = unique(sample_info$stage))

# make counts matrix
counts <-
    count_data[, grepl(' count$', colnames(count_data)) &
                !grepl('normalised count$', colnames(count_data))]
# remove ' count' from column names
colnames(counts) <-
  gsub(' count$', '', colnames(counts))
rownames(counts) <- count_data[, 'Gene ID']
# order counts by samples file
counts <- counts[ , rownames(sample_info) ]

# make normalised counts matrix
if (cmd_line_args$options[['include_norm_counts']]){
    normalised_counts <-
        count_data[, grepl('normalised count$', colnames(count_data)) ]
    # remove ' count' from column names
    colnames(normalised_counts) <-
      gsub(' normalised count$', '', colnames(normalised_counts))
    rownames(normalised_counts) <- count_data[, 'Gene ID']
    # order normalised_counts by samples file
    normalised_counts <- normalised_counts[ , rownames(sample_info) ]
}

# if sample_names column exists, rename
if ( any( names(sample_info) == "sample_name" ) ) {
    # move rownames to sample id column
    sample_info[['sample_id']] <- rownames(sample_info)
    # set rownames to sample names
    rownames(sample_info) <- sample_info[['sample_name']]
    # and also colnames of counts
    colnames(counts) <- sample_info[['sample_name']]
    if (cmd_line_args$options[['include_norm_counts']]){
        colnames(normalised_counts) <- sample_info[['sample_name']]
    }
}

# row data. chr, start etc for each gene
baseline_row_data <- count_data[, c('Gene ID', 'Chr', 'Start', 'End',
                                    'Strand', 'Biotype', 'Name',
                                    'Description')]

if (cmd_line_args$options[['include_norm_counts']]){
    baseline_data <- 
        SummarizedExperiment(assays = list(counts = as.matrix(counts),
                                           normalised_counts = as.matrix(normalised_counts)),
                            rowData = DataFrame(baseline_row_data),
                            colData = DataFrame(sample_info)
                          )
} else {
    baseline_data <- 
        SummarizedExperiment(assays = list(counts = as.matrix(counts)),
                            rowData = DataFrame(baseline_row_data),
                            colData = DataFrame(sample_info)
                          )
}

# load metadata if it exists
## want this object to be named using the assembly and ensembl_version
#sp <- sub("[a-z]+_", "", cmd_line_args$options[['species']])
#sp <- sub("([a-z])[a-z]+$", "\\1", sp)
#object_name <- paste(sp, cmd_line_args$options[['assembly']],
#                     paste0('e', cmd_line_args$options[['ensembl_version']]),
#                     'baseline', sep ="_")
#print(object_name)
if (!is.null(cmd_line_args$options[['metadata']])) {
    baseline_metadata <- yaml.load_file(cmd_line_args$options[['metadata']])
} else {
    baseline_metadata <- list(
        Species = cmd_line_args$options[['species']],
        Assembly = cmd_line_args$options[['assembly']],
        Ensembl_version = cmd_line_args$options[['ensembl_version']]
    )
}
metadata(baseline_data) <- baseline_metadata

# write out to rdata file
output_file <- cmd_line_args$args[3]
saveRDS(baseline_data, file = output_file, compress = TRUE)
