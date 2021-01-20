### R script_Annotating VCF file(s)
rm(list = ls())
# start.time <- Sys.time()
cat("# R script for annotating VCF files", sep = "\n")

### Loading VariantAnnotation package
cat("Loading libraries ...", "\n")
suppressMessages(library(VariantAnnotation, verbose = F))

### Reading VCF file
cat("Reading VCF file ...", "\n")
vcf <- readVcf("intermediate_sh/merged_clean.vcf.gz")

### Annotation
cat("Annotation ...", sep = "\n")
var.summary <- snpSummary(vcf) # Counts and distribution statistics for SNPs.
nAA <- var.summary$g00 # Counts for genotype 00 (homozygous reference)
nAB <- var.summary$g01 # Counts for genotype 01 or 10 (heterozygous)
nBB <- var.summary$g11 # Counts for genotype 11 (homozygous alternate)
p <- round(var.summary$a0Freq, 4) # Frequency of the reference allele.
q <- round(var.summary$a1Freq, 4) # Frequency of the alternate allele

### Flags calculation for INFO field
AC <- (nAB) + (nBB*2) # Alternate allele count (AC)
AF <- q # Alternate allele frequency (AF)
MAF <- pmin(p, q) # Minor allele frequency (MAF)
AN <- (nAA + nAB + nBB)*2 # Total number of alleles (AN)
AVG_DP <- round(rowMeans(geno(vcf)$DP),4) # Average depth (AVG_DP)
N_ALT_HET <- nAB # Number of heterozygous alternate allele carriers (N_ALT_HET)
N_ALT_HOM <- nBB # Number of homozygous alternate allele carriers (N_ALT_HOM)

# AVG_DP_P
DP <- geno(vcf)$DP
DP_mean_unsorted <- rowMeans(DP)
DP_mean_sorted <- sort(DP_mean_unsorted)
N <- length(DP_mean_sorted)
AVG_DP_P <- as.numeric()

# Loop - Calculating percentile
for (j in 1:N) {  # cat("SNP_Info:", names(DP_mean_sorted)[j], "\n")
  AVG_DP_P[j] <- ((j-1)/N)*100
}

names(AVG_DP_P) <- names(DP_mean_sorted)
AVG_DP_P <- round(AVG_DP_P[order(match(names(AVG_DP_P), names(DP_mean_unsorted)))], 4) # Average depth's percentile value (AVG_DP_P)

### Adding New Flags to INFO Fields
flags <- data.frame(AC=AC, AF=AF, MAF=MAF, AN=AN, AVG_DP=AVG_DP, N_ALT_HET=N_ALT_HET, N_ALT_HOM=N_ALT_HOM, AVG_DP_P) # Making Data fram using calculated flags
suppressWarnings(info(vcf) <- cbind(info(vcf), flags))# Adding new flags data to existing INFO field data in vcf file. Please ignore the warning message.
# info(vcf)

### Adding metadata for newly added INFO flags in header
new_header <- read.csv("data/info_files/flag_description.csv", row.names = 1) # Read new flags description from a file
info(header(vcf)) <- rbind(info(header(vcf)), new_header) # Adding metadata to existing header metadata in vcf file.
# info(header(vcf))

### Removing contings info from the header that has no data in the vcf file
all_contigs <- header(vcf)@header@listData[["contig"]] # Header Info only for contigs
selected.contigs <- data.frame(length=all_contigs[19:22,], row.names = rownames(all_contigs)[19:22]) # Make a dataframe using selected contigs
header(vcf)@header@listData[["contig"]] <- selected.contigs # Adding selected contig data to VCF file contig info
header(vcf)@reference <- rownames(selected.contigs) # Adding selected contig names to VCF file reference info
#header(vcf)@header@listData[["contig"]]

### write VCF
cat("Writing output VCF file ...", "\n")
suppressWarnings(writeVcf(vcf, filename="intermediate_sh/merged_clean_annotated.vcf.gz",  index=T))
# end.time <- Sys.time()
# time.taken <- end.time - start.time
# cat("All DONE! Run Time:", round(time.taken, 2), "sec.", "\n")

### END
