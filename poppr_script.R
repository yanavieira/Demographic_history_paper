library(vcfR)
library(poppr)
vcf_file <- system.file("extdata", "pinf_sc50.vcf.gz", package = "pinfsc50")

#open vcf_file of ck_dataset
vcf <- read.vcfR("/home/yana/R/x86_64-pc-linux-gnu-library/3.2/poppr/files/c94d5m5p3r5_MM80_Maf_only_ck.recode.vcf", verbose = FALSE)
vcf

#convert vcf_file to a light object
vcf_file <- vcfR2genlight(vcf)

#atribute populations to the isolates
pop(vcf_file) <- c("popAng","popAng", "popAng", "popEast", "popAng", "popEast", "popCam", "popCam", "popCam", "popCam", "popEast", "popEast", "popEast", "popEast", "popEast", "popEast", "popEast", "popEast", "popEast", "popEast", "popEast", "popEast", "popEast", "popEast", "popEast", "popEast", "popEast", "popEast", "popEast", "popEast")

#Because Ck is a haploid pathogen
ploidy(vcf_file)<- 1

vcf_file

#convert vcf_file to a different type of object
vcf_file2 <- as.snpclone(vcf_file)

#filter applied to detect the clonal lineages and perfrom the clone correction
mlg.filter(vcf_file2)<-0.02


#calculation of index association in ck_dataset with clone correction
sampia_cc <- samp.ia(vcf_file2)
sampia_cc

#convert vcf_file to a different type of object
vcf_file3 <- as.snpclone(vcf_file)

#calculation of index association in ck_dataset without clone correction
sampia <- samp.ia(vcf_file3)
sampia

#boxplot for Index association with and without clone correction
boxplot(data.frame(clone_correct = sampia_cc, no_clone_correction = sampia), ylab = "Index of Association", main = "Index of association over 100 random SNPs")

#distance calculation for both datasets
dis2 <- bitwise.dist(vcf_file2, percent = TRUE, mat = FALSE, missing_match = TRUE, differences_only = FALSE, threads = 0)

dis3 <- bitwise.dist(vcf_file3, percent = TRUE, mat = FALSE, missing_match = TRUE, differences_only = FALSE, threads = 0)

#minimum spanning networks with clone correction
poppr.msn(vcf_file2, distmat = dis2)

#minimum spanning networks without clone correction
poppr.msn(vcf_file3, distmat = dis3)

