
source('lib.R')


#NAME <- 'H3K4me3_GM12878.ENCFF023LTU.hg19'
#NAME <- 'H3K4me3_GM12878.ENCFFo23LTU.hg38'
#NAME <- 'H3K4me3_GM12878.ENCFF432EMI.hg19'
NAME <- 'H3K4me3_GM12878.ENCFF432EMI.hg38'
#OUT_DIR <- '/Users/melanie/Desktop/bioinf/project/gitproject/hse21_H3K4me3_G4_human/images/'

#Убрать сильно длинные участки, если есть

bed_df <- read.delim(paste0('data_dir', NAME, '.bed'), as.is = TRUE, header = FALSE)
colnames(bed_df) <- c('chrom', 'start', 'end', 'name', 'score')
bed_df$len <- bed_df$end - bed_df$start
head(bed_df)

bed_df <- bed_df %>%
  arrange(-len) %>%
  filter(len < 50000)

ggplot(bed_df) +
  aes(x = len) +
  geom_histogram() +
  ggtitle(NAME, subtitle = sprintf('Number of peaks = %s', nrow(bed_df))) +
  theme_bw()
ggsave(paste0('filtered.', NAME, '.filtered.pdf'), path = OUT_DIR)

bed_df %>%
  select(-len) %>%
  write.table(file=data_dir_filtered,
              col.names = FALSE, row.names = FALSE, sep = '\t', quote = FALSE)
