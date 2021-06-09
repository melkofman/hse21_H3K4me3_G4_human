#chipSeeker 

source('lib.R')

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ChIPseeker")
BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene")
BiocManager::install("clusterProfiler")
BiocManager::install("org.Hs.eg.db")
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(clusterProfiler)

#NAME <- 'H3K4me3_GM12878.ENCFF023LTU.hg19.filtered'
NAME <- 'H3K4me3_GM12878.ENCFF432EMI.hg19.filtered'
#NAME <- 'H3K4me3_GM12878.ENCFF432EMI.hg38.filtered'
#NAME <- 'H3K4me3_GM12878.ENCFFo23LTU.hg38.filtered'


DATA_DIR <- data_dir_filtered
BED_FN <- paste0(DATA_DIR, NAME, '.bed')

###

txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

peakAnno <- annotatePeak(BED_FN, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Hs.eg.db")

#pdf(paste0(OUT_DIR, 'chip_seeker.', NAME, '.plotAnnoPie.pdf'))
png(paste0(OUT_DIR, 'chip_seeker.', NAME, '.plotAnnoPie.png'))
plotAnnoPie(peakAnno)
dev.off()


# пики на хромосомках
peak <- readPeakFile(BED_FN)
pdf(paste0(OUT_DIR, 'chip_seeker.', NAME, '.covplot.pdf'))
covplot(peak, weightCol="V5")
dev.off()
# 
