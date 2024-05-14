
################################################################################
#################### Generate MiXeR-Pred Univariate Output #####################
################################################################################

rm(list=ls())

## install packages if not installed already and load them
if(!require(tidyverse)){
  install.packages("tidyverse")
  library(tidyverse)
}

if(!require(data.table)){
  install.packages("data.table")
  library(data.table)
}

if(!require(argparser)){
  install.packages("argparser")
  library(argparser)
}

## Define command line input options 
par <- arg_parser("Generate Univariate MiXeR-Pred Estimates File")

par <- add_argument(
  par, "--ref", help="path to hrc reference .bim files used by MiXeR (usually located in the reference folder of the docker/singularity image)")
par <- add_argument(
  par, "--outdir", help="path to output directory")
par <- add_argument(
  par, "--indir", help="path to intput directory (i.e., location of MiXeR .snps files)")
par <- add_argument(
  par, "--snps-file-prefix", help="file prefix for .snps files (e.g., tairt1_vs_trait2)")
par <- add_argument(
  par, "--target-bim", help="path to directory .bim file used for target sample for polygenic scoring. Used hear to filter SNPs(e.g., /full/path/target.bim)")
par <- add_argument(
  par, "--plink-dir", help="path to directory housing plink binary")
par <- add_argument(
  par, "--bfile", help="path to directory housing .bed .bim .fam files (used for clumping) along with the bfile file prefix (e.g., /full/path/1KG_EUR)")
par <- add_argument(
  par, "--snp-thresholds", nargs=Inf,
  help="number of top clumped SNPs (in thousands) to include in output file for PRS (e.g., 1 2 3 10 20 100)")

parsed <- parse_args(par)

ref <- parsed$ref
outdir <- parsed$outdir
indir <- parsed$indir
snpsfiles <- parsed$snps_file_prefix
targetbim <- parsed$target_bim
plinkdir <- parsed$plink_dir
bfile <- parsed$bfile
snpthresholds <- as.integer(parsed$snp_thresholds)

dir.create(outdir)
setwd(outdir)

## read in reference for MiXeR SNP data + order
filenames_bim2 <- c(
  paste0(ref, "/hrc_chr1_EUR_qc.bim"), paste0(ref, "/hrc_chr2_EUR_qc.bim"), 
  paste0(ref, "/hrc_chr3_EUR_qc.bim"),  paste0(ref, "/hrc_chr4_EUR_qc.bim"),
  paste0(ref, "/hrc_chr5_EUR_qc.bim"),  paste0(ref, "/hrc_chr6_EUR_qc.bim"), 
  paste0(ref, "/hrc_chr7_EUR_qc.bim"),  paste0(ref, "/hrc_chr8_EUR_qc.bim"),
  paste0(ref, "/hrc_chr9_EUR_qc.bim"),  paste0(ref, "/hrc_chr10_EUR_qc.bim"),
  paste0(ref, "/hrc_chr11_EUR_qc.bim"),  paste0(ref, "/hrc_chr12_EUR_qc.bim"),
  paste0(ref, "/hrc_chr13_EUR_qc.bim"),  paste0(ref, "/hrc_chr14_EUR_qc.bim"),
  paste0(ref, "/hrc_chr15_EUR_qc.bim"),  paste0(ref, "/hrc_chr16_EUR_qc.bim"),
  paste0(ref, "/hrc_chr17_EUR_qc.bim"),  paste0(ref, "/hrc_chr18_EUR_qc.bim"),
  paste0(ref, "/hrc_chr19_EUR_qc.bim"),  paste0(ref, "/hrc_chr20_EUR_qc.bim"),
  paste0(ref, "/hrc_chr21_EUR_qc.bim"),  paste0(ref, "/hrc_chr22_EUR_qc.bim")
)

df_snp <- rbindlist(sapply(filenames_bim2, fread, simplify = FALSE),
                    use.names=FALSE, idcol = 'filename')
df_snp <- df_snp %>%
  rename(CHR=V1,RSID=V2,PosinMorCM=V3,POS=V4,A1=V5,A2=V6)

## Read in MiXeR snp files + merge with reference + generate estimates
trait_files <- list.files(
  path = indir,
  pattern = glob2rx(paste0(snpsfiles, "*test1.snps.csv")), full.names = T)


process_files <- lapply(trait_files, function(file){
  
  run_file <- fread(file)
  run_file_snpinfo <- cbind(df_snp, run_file)
  
  # Generate posterior estimates
  run_file_snpinfo <- run_file_snpinfo %>%
    mutate(Ed=posteriordeltaC1/posteriordeltaC0) %>%  
    
    #Generate posterior significance estimates, exp(-x) where x=E(delta1^2 | z1, z2).
    mutate(Edtms2=posteriordeltaC2/posteriordeltaC0) %>% # E(delta1^2 | z1, z2) = C20/C00
    mutate(expEd=exp(-1*Edtms2))%>% #exp(-x), exp computes the exponential function.
    
    # select necessary variables -- Ask Alex which are necessary...
    select(c(RSID,CHR,POS,A1,A2,Ed,expEd))
  
})


## summarize estimates across MiXeR runs
MiXeR_output <- rbindlist(process_files, use.names=FALSE, idcol = 'run')

MiXeR_output <- MiXeR_output %>%
  group_by(RSID) %>%
  summarise(
    CHR = first(CHR), POS = first(POS),
    A1 = first(A1), A2 = first(A2),
    Ed_mean = mean(Ed, na.rm = TRUE),
    expEd_mean = mean(expEd, na.rm = TRUE)
  ) %>%
  ungroup()


## clean up MiXeR sumstats
# remove SNPs with NA for these columns
# (Ed_mean, expEd_mean)
MiXeR_output <- MiXeR_output[rowSums(is.na(MiXeR_output[,6:7]))==0,]

# remove ambiguous SNPs
MiXeR_output <- MiXeR_output%>%
  filter (! (A1=="C" & A2=="G"))%>%
  filter (! (A1=="G" & A2=="C"))%>%
  filter (! (A1=="A" & A2=="T"))%>%
  filter (! (A1=="T" & A2=="A"))

# standard header
names(MiXeR_output) <- c("RSID", "CHR", "POS", "A1", "A2", "Ed", "expEd")

# Remove duplicate SNP and CHR_POS
MiXeR_output <- MiXeR_output[duplicated(MiXeR_output$RSID) ==F,]
MiXeR_output$CHR_POS <- paste0(MiXeR_output$CHR, "_", MiXeR_output$POS)
MiXeR_output <- MiXeR_output[duplicated(MiXeR_output$CHR_POS) ==F,]
MiXeR_output <- dplyr::select(
  MiXeR_output, "RSID", "CHR", "POS", "A1", "A2", "Ed", "expEd")

## Write out MiXeR file
# outfile <- gsub("vs", "on", snpsfiles)
fwrite(MiXeR_output, file=paste0(outdir, "/", snpsfiles, "_mixer_estimates.gz"),
       sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE, na = NA)

## Filter sumstats for target sample SNPs
targetSNPS <- fread(targetbim, header=F)
MiXeR_output_filtered <- MiXeR_output[MiXeR_output$RSID %in% targetSNPS$V2,]

# Remove duplicate vid and filter for target vid
targetSNPS$vid <- paste0(targetSNPS$V1, "_", targetSNPS$V4, "_",
                         targetSNPS$V5, "_", targetSNPS$V6)
targetSNPS$vid2 <- paste0(targetSNPS$V1, "_", targetSNPS$V4, "_",
                         targetSNPS$V6, "_", targetSNPS$V5)

target_vid <- c(targetSNPS$vid, targetSNPS$vid2)

MiXeR_output_filtered$vid <- paste0(
  MiXeR_output_filtered$CHR, "_", MiXeR_output_filtered$POS, "_",
  MiXeR_output_filtered$A1, "_", MiXeR_output_filtered$A2)
MiXeR_output_filtered <- MiXeR_output_filtered[duplicated(MiXeR_output_filtered$vid) ==F,]
MiXeR_output_filtered <- MiXeR_output_filtered[MiXeR_output_filtered$vid %in% target_vid,]

MiXeR_output_filtered <- dplyr::select(
  MiXeR_output_filtered, "RSID", "CHR", "POS", "A1", "A2", "Ed", "expEd")

fwrite(MiXeR_output_filtered, 
       file=paste0(outdir, "/", snpsfiles, "_mixer_estimates_target_filtered.gz"),
       sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE, na = NA)

### PLINK clumping 
system(paste0(
  plinkdir,"/plink", " --bfile ", bfile, " --clump-p1 1", " --clump-r2 0.2",
  " --clump-kb 250", " --clump ", paste0(outdir, "/", snpsfiles, "_mixer_estimates_target_filtered.gz"),
  " --clump-snp-field RSID", " --clump-field expEd",
  " --out ", paste0(outdir, "/", snpsfiles, "_target_filtered")))

## Read in clumped file and generate p-thresholded sumstat files for PRS
dir.create(paste0(outdir, "/", snpsfiles, "_clump_topSNPs"))

lapply(c(snpthresholds), function(n_SNPs){
  
  mixer_out <- fread(paste0(outdir, "/", snpsfiles, "_mixer_estimates_target_filtered.gz"))
  indep_snp <- fread(paste0(outdir, "/", snpsfiles, "_target_filtered.clumped"))
  mixer_out <- mixer_out[mixer_out$RSID %in% indep_snp$SNP,]
  mixer_out <- mixer_out[order(mixer_out$expEd),]
  mixer_out <- mixer_out[1:(n_SNPs*1000),]
  
  write.table(
    mixer_out,file=gzfile(
      paste0(paste0(outdir, "/", snpsfiles, "_clump_topSNPs/"), 
             snpsfiles, "_top_", n_SNPs, "K.mixer.gz")),
    quote = F, sep = "\t", row.names = F)
  
})





