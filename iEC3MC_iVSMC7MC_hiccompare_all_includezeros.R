#make sure all files are in the current folder and current folder is the working directory

library(HiCcompare)

##use parallel for multiple chromosome
library(BiocParallel)

register(SerialParam())

##########################
#https://bioconductor.org/packages/release/bioc/vignettes/HiCcompare/inst/doc/HiCcompare-vignette.html

#Extracting data from .cool files
iEC3MC <- cooler2bedpe(path = "iEC3MC_matrix_25kb.cool")
iVSMC7MC <- cooler2bedpe(path = "iVSMC7MC_matrix_25kb.cool")


#We can also list multiple hic.tables in order to utilize parallel computing. First we create another hic.table as was done above. Then we combine these two hic.tables into a list.

#https://bioconductor.org/packages/release/bioc/vignettes/HiCcompare/inst/doc/HiCcompare-vignette.html

##mind include pairwise interactions where one of the interaction frequencies is 0.
#run with test dataset from https://rdrr.io/bioc/HiCcompare/man/hic_loess.html, found the final test with hiccompare will test those with one of interaction frequencies is 0. But after adjust of the IF, they are not 0 anymore.

chr1<-create.hic.table(iEC3MC$cis$chr1, iVSMC7MC$cis$chr1, include.zeros = TRUE)
chr2<-create.hic.table(iEC3MC$cis$chr2, iVSMC7MC$cis$chr2, include.zeros = TRUE)
chr3<-create.hic.table(iEC3MC$cis$chr3, iVSMC7MC$cis$chr3, include.zeros = TRUE)
chr4<-create.hic.table(iEC3MC$cis$chr4, iVSMC7MC$cis$chr4, include.zeros = TRUE)
chr5<-create.hic.table(iEC3MC$cis$chr5, iVSMC7MC$cis$chr5, include.zeros = TRUE)
chr6<-create.hic.table(iEC3MC$cis$chr6, iVSMC7MC$cis$chr6, include.zeros = TRUE)
chr7<-create.hic.table(iEC3MC$cis$chr7, iVSMC7MC$cis$chr7, include.zeros = TRUE)
chr8<-create.hic.table(iEC3MC$cis$chr8, iVSMC7MC$cis$chr8, include.zeros = TRUE)
chr9<-create.hic.table(iEC3MC$cis$chr9, iVSMC7MC$cis$chr9, include.zeros = TRUE)
chr10<-create.hic.table(iEC3MC$cis$chr10, iVSMC7MC$cis$chr10, include.zeros = TRUE)
chr11<-create.hic.table(iEC3MC$cis$chr11, iVSMC7MC$cis$chr11, include.zeros = TRUE)
chr12<-create.hic.table(iEC3MC$cis$chr12, iVSMC7MC$cis$chr12, include.zeros = TRUE)
chr13<-create.hic.table(iEC3MC$cis$chr13, iVSMC7MC$cis$chr13, include.zeros = TRUE)
chr14<-create.hic.table(iEC3MC$cis$chr14, iVSMC7MC$cis$chr14, include.zeros = TRUE)
chr15<-create.hic.table(iEC3MC$cis$chr15, iVSMC7MC$cis$chr15, include.zeros = TRUE)
chr16<-create.hic.table(iEC3MC$cis$chr16, iVSMC7MC$cis$chr16, include.zeros = TRUE)
chr17<-create.hic.table(iEC3MC$cis$chr17, iVSMC7MC$cis$chr17, include.zeros = TRUE)
chr18<-create.hic.table(iEC3MC$cis$chr18, iVSMC7MC$cis$chr18, include.zeros = TRUE)
chr19<-create.hic.table(iEC3MC$cis$chr19, iVSMC7MC$cis$chr19, include.zeros = TRUE)
chr20<-create.hic.table(iEC3MC$cis$chr20, iVSMC7MC$cis$chr20, include.zeros = TRUE)
chr21<-create.hic.table(iEC3MC$cis$chr21, iVSMC7MC$cis$chr21, include.zeros = TRUE)
chr22<-create.hic.table(iEC3MC$cis$chr22, iVSMC7MC$cis$chr22, include.zeros = TRUE)
chrX<-create.hic.table(iEC3MC$cis$chrX, iVSMC7MC$cis$chrX, include.zeros = TRUE)
chrY<-create.hic.table(iEC3MC$cis$chrY, iVSMC7MC$cis$chrY, include.zeros = TRUE)

#setup experiment

hic.list <- list(chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY)

#normalize

hic.list <- hic_loess(hic.list, parallel=TRUE)

#analyze

hic.list <- hic_compare(hic.list, A.min = NA, adjust.dist = TRUE, p.method = 'fdr', parallel = TRUE)

#save the object of hic.list

save(hic.list,  file = "iEC3MC_iVSMC7MC_hiccompare_withzero_25kb_afteranalysis.RData")

#make result dataframe

hic.list <- do.call (rbind, hic.list)

#get significant results

hic.list <- hic.list[hic.list$p.adj < 0.05,]

#write results

write.csv(hic.list, file = "iEC3MC_iVSMC7MC_hiccompare_allchromosome_sigadjp_25Kresolution_withzeros.csv")



