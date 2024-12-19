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

chr1<-create.hic.table(iEC3MC$cis$chr1, iVSMC7MC$cis$chr1)
chr2<-create.hic.table(iEC3MC$cis$chr2, iVSMC7MC$cis$chr2)
chr3<-create.hic.table(iEC3MC$cis$chr3, iVSMC7MC$cis$chr3)
chr4<-create.hic.table(iEC3MC$cis$chr4, iVSMC7MC$cis$chr4)
chr5<-create.hic.table(iEC3MC$cis$chr5, iVSMC7MC$cis$chr5)
chr6<-create.hic.table(iEC3MC$cis$chr6, iVSMC7MC$cis$chr6)
chr7<-create.hic.table(iEC3MC$cis$chr7, iVSMC7MC$cis$chr7)
chr8<-create.hic.table(iEC3MC$cis$chr8, iVSMC7MC$cis$chr8)
chr9<-create.hic.table(iEC3MC$cis$chr9, iVSMC7MC$cis$chr9)
chr10<-create.hic.table(iEC3MC$cis$chr10, iVSMC7MC$cis$chr10)
chr11<-create.hic.table(iEC3MC$cis$chr11, iVSMC7MC$cis$chr11)
chr12<-create.hic.table(iEC3MC$cis$chr12, iVSMC7MC$cis$chr12)
chr13<-create.hic.table(iEC3MC$cis$chr13, iVSMC7MC$cis$chr13)
chr14<-create.hic.table(iEC3MC$cis$chr14, iVSMC7MC$cis$chr14)
chr15<-create.hic.table(iEC3MC$cis$chr15, iVSMC7MC$cis$chr15)
chr16<-create.hic.table(iEC3MC$cis$chr16, iVSMC7MC$cis$chr16)
chr17<-create.hic.table(iEC3MC$cis$chr17, iVSMC7MC$cis$chr17)
chr18<-create.hic.table(iEC3MC$cis$chr18, iVSMC7MC$cis$chr18)
chr19<-create.hic.table(iEC3MC$cis$chr19, iVSMC7MC$cis$chr19)
chr20<-create.hic.table(iEC3MC$cis$chr20, iVSMC7MC$cis$chr20)
chr21<-create.hic.table(iEC3MC$cis$chr21, iVSMC7MC$cis$chr21)
chr22<-create.hic.table(iEC3MC$cis$chr22, iVSMC7MC$cis$chr22)
chrX<-create.hic.table(iEC3MC$cis$chrX, iVSMC7MC$cis$chrX)
chrY<-create.hic.table(iEC3MC$cis$chrY, iVSMC7MC$cis$chrY)

#setup experiment

hic.list <- list(chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY)

#normalize

hic.list <- hic_loess(hic.list, parallel=TRUE)

#analyze

hic.list <- hic_compare(hic.list, A.min = NA, adjust.dist = TRUE, p.method = 'fdr', parallel = TRUE)

#make result dataframe

hic.list <- do.call (rbind, hic.list)

#get significant results

hic.list <- hic.list[hic.list$p.adj < 0.05,]

#write results

write.csv(hic.list, file = "iEC3MC_iVSMC7MC_hiccompare_allchromosome_sigadjp_25Kresolution.csv")



