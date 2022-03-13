BiocManager::install("ChIPseeker")
BiocManager::install("ChIPpeakAnno")
BiocManager::install("diffloop")
library(ChIPseeker)
library(diffloop)
library(ChIPpeakAnno)
# BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene")
library(TxDb.Hsapiens.UCSC.hg19.knownGene)

#read file
test1 <- readPeakFile("/Users/davidseong/Box/David Seong's Files/Classwork/Stemrem 205/project files/GSE129646_BRCA_raw_counts.peaks_mda_brm_rep1_sort_dedup.counts.txt")
#annotate peaks
annotatedpeaks1 <- annotatePeak(test1, TxDb = TxDb.Hsapiens.UCSC.hg19.knownGene)
plotAnnoPie(annotatedpeaks1)
MacsGranges_Anno <- as.GRanges(annotatedpeaks1)
TSS_brm1 <- MacsGranges_Anno[abs(MacsGranges_Anno$distanceToTSS) < 
                               500]
#order dataframe of counts in descending order of peak counts and take top 20
df <- as.data.frame(TSS_brm1)
orddf <- df[order(-df$X20),]
top20 <- orddf[1:20,]

#converted gene names for top 20 genes (gene names converted using SynGO ID: https://www.syngoportal.org/convert.html)
atac20list <- c("AFTPH", "CREB3", "E2F3","EFHC1", "HMGCS1", "KHSRP", "KLHDC3", "LINC01011", "MAP3K11", "MBOAT1", "MLEC", "PIM1", "RALBP1", "RPP40", "ZNF443")



#gsea
tsvfile <- read.table("/Users/dhs37929/Downloads/gsea_report_for_na_pos_1646771646033.tsv", sep = "\t", header = TRUE)

tsvfile$NAME <- factor(tsvfile$NAME, levels = tsvfile$NAME[order(tsvfile$NES)])

ggplot(tsvfile, aes(x=NES, y=NAME, size = SIZE, color = FDR.q.val)) +
  scale_color_viridis(direction = -1) +
  geom_point(alpha=1)
