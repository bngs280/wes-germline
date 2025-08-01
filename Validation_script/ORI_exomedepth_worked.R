#!/usr/bin/env Rscript
## original_worekd
# Set up command line options
suppressMessages({
  library(ExomeDepth)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(optparse)
})

# Command line argument parsing
option_list <- list(
  make_option(c("-d", "--dir"), type="character", default="./data/BAM/",
              help="Directory containing BAM files", metavar="DIR"),
  make_option(c("-t", "--test"), type="character",
              help="Test sample name (without .bam extension)", metavar="SAMPLE"),
  make_option(c("-r", "--reference"), type="character",
              help="Comma-separated reference samples or 'default' to use pre-specified references", metavar="REFS"),
  make_option(c("-f", "--refdir"), type="character", default="./data/reference_BAMs/",
              help="Directory containing pre-specified reference BAMs", metavar="REFDIR"),
  make_option(c("-o", "--outdir"), type="character", default="./results/",
              help="Output directory for results", metavar="OUTDIR")
)

parser <- OptionParser(option_list=option_list)
args <- parse_args(parser)

# Validate arguments
if (is.null(args$test) || is.null(args$reference)) {
  print_help(parser)
  stop("Both --test and --reference must be specified", call.=FALSE)
}

# Set working directory and create output dir
setwd("/home/vgenomics/RND/CNVGerm/")
dir.create(args$outdir, showWarnings = FALSE)

# Get exon coordinates
message("Getting exon coordinates...")
hg38.txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
exons.hg38 <- exons(hg38.txdb)
exons.hg38 <- reduce(exons.hg38)
exons.hg38.df <- data.frame(
  chromosome = seqnames(exons.hg38),
  start = start(exons.hg38),
  end = end(exons.hg38),
  name = paste0("exon_", 1:length(exons.hg38))
)

# Filter to standard chromosomes
standard_chroms <- paste0("chr", c(1:22, "X", "Y", "M"))
exons.hg38.filtered <- exons.hg38.df[exons.hg38.df$chromosome %in% standard_chroms, ]

# Function to get BAM paths
get_bam_paths <- function(dir, samples) {
  bam_files <- file.path(dir, paste0(samples, ".bam"))
  missing <- !file.exists(bam_files)
  if (any(missing)) {
    stop("Missing BAM files: ", paste(bam_files[missing], collapse=", "))
  }
  return(bam_files)
}

# Handle reference samples
if (args$reference == "default") {
  ref_samples <- tools::file_path_sans_ext(list.files(args$refdir, pattern="\\.bam$"))
  message("Using default reference samples from ", args$refdir, ": ", paste(ref_samples, collapse=", "))
  ref_bams <- get_bam_paths(args$refdir, ref_samples)
} else {
  ref_samples <- unlist(strsplit(args$reference, ","))
  ref_bams <- get_bam_paths(args$dir, ref_samples)
}

# Get test sample BAM
test_bam <- get_bam_paths(args$dir, args$test)

# Combine all BAMs needed
all_bams <- c(test_bam, ref_bams)

# Get read counts
message("\nCounting reads in BAM files...")
read.counts <- getBamCounts(
  bed.frame = exons.hg38.filtered,
  bam.files = all_bams,
  min.mapq = 20,
  include.chr = FALSE
)

# Set test and reference samples
test_col <- paste0(args$test, ".bam")
ref_cols <- paste0(ref_samples, ".bam")

my.test <- read.counts[[test_col]]
my.reference.set <- as.matrix(read.counts[, ref_cols])

# Remove zero-length bins
bin.length <- (read.counts$end - read.counts$start)/1000
zero.bins <- which(bin.length == 0)

if (length(zero.bins) > 0) {
  message("Removing ", length(zero.bins), " bins with zero length")
  read.counts <- read.counts[-zero.bins, ]
  my.test <- my.test[-zero.bins]
  my.reference.set <- my.reference.set[-zero.bins, , drop = FALSE]
  bin.length <- bin.length[-zero.bins]
}

# Select optimal reference set
message("\nSelecting optimal reference set...")
my.choice <- select.reference.set(
  test.counts = my.test,
  reference.counts = my.reference.set,
  bin.length = bin.length,
  n.bins.reduced = 10000
)

message("Selected references: ", paste(unique(my.choice$reference.choice), collapse=", "))

# Create reference counts
my.matrix <- as.matrix(read.counts[, my.choice$reference.choice, drop = FALSE])
my.reference.selected <- rowSums(my.matrix)

# Create ExomeDepth object and call CNVs
message("\nRunning CNV detection...")
all.exons <- new('ExomeDepth',
                 test = my.test,
                 reference = my.reference.selected,
                 formula = 'cbind(test, reference) ~ 1')

all.exons <- CallCNVs(
  x = all.exons,
  transition.probability = 10^-4,
  chromosome = read.counts$chromosome,
  start = read.counts$start,
  end = read.counts$end,
  name = read.counts$exon
)

# Add copy number classification
classify_copy_number <- function(reads_ratio) {
  cut(reads_ratio,
      breaks = c(-Inf, 0.10, 0.75, 1.25, 1.75, 2.25, Inf),
      labels = c("0 (homozygous deletion)",
                 "1 (heterozygous deletion)",
                 "2 (normal)",
                 "3 (heterozygous duplication)",
                 "4",
                 "OTHER"),
      right = FALSE)
}

# Create output files
output_base <- file.path(args$outdir, args$test)
ranked <- all.exons@CNV.calls[order(all.exons@CNV.calls$BF, decreasing = TRUE),]
ranked$Copy_number <- classify_copy_number(ranked$reads.ratio)

# Save results
ranked_file <- paste0(output_base, "_exome_calls_ranked.csv")
#detailed_file <- paste0(output_base, "_exome_calls_annotated.csv")

write.csv(ranked, file = ranked_file, row.names = FALSE)

# # Create detailed version with CNV type
# ranked$CNV_type <- ifelse(
#   as.numeric(sub(" .*", "", ranked$Copy_number)) < 2, 
#   "Deletion", 
#   ifelse(as.numeric(sub(" .*", "", ranked$Copy_number)) > 2, 
#          "Duplication", 
#          "Normal")
# )
# write.csv(ranked, file = detailed_file, row.names = FALSE)

message("\nAnalysis complete. Results saved to:")
message("- Ranked calls: ", ranked_file)
#message("- Annotated calls: ", detailed_file)