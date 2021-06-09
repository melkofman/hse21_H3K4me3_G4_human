
source('lib.R')

#NAME <- 'H3K4me3_GM12878.ENCFF023LTU.hg19'
#NAME <- 'H3K4me3_GM12878.ENCFFo23LTU.hg38'
#NAME <- 'H3K4me3_GM12878.ENCFF432EMI.hg19'
NAME <- 'H3K4me3_GM12878.ENCFF432EMI.hg38'


bed_df <- read.delim(paste0(data_dir, NAME, '.bed'), as.is = TRUE, header = FALSE)
head(bed_df)
colnames(bed_df) <- c('chrom', 'start', 'end', 'name', 'score')
bed_df$len <- bed_df$end - bed_df$start

ggplot(bed_df) +
  aes(x = len) +
  geom_histogram() +
  ggtitle(NAME, subtitle = sprintf('Number of peaks = %s', nrow(bed_df))) +
  theme_bw()
ggsave(paste0('len_hist.', NAME, 'no_filter.pdf'), path = OUT_DIR)

