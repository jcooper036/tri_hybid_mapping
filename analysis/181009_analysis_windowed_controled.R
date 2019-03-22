library(readr)
library(ggplot2)
library(plotrix)
library(gridExtra)

## for reading TSV
rdTSV <- function(fileName) read_delim(fileName, "\t", escape_double = FALSE, col_types = cols(mel = col_number(), sec = col_number(), sim = col_number()),trim_ws = TRUE)

## read all the TSVs, bind them by sex, raw vcfs
rep1 <- rdTSV("~/Desktop/genome_analysis/controled_windows/181009_adjusted_rep1_.tsv")
rep2 <- rdTSV("~/Desktop/genome_analysis/controled_windows/181009_adjusted_rep2_.tsv")
rep3 <- rdTSV("~/Desktop/genome_analysis/controled_windows/181009_adjusted_rep3_.tsv")

## merge the data together
repMerge <- merge(rep1, rep2, by=c("CHROM", "POS"))
repMerge <- merge(repMerge, rep3, by=c("CHROM", "POS"))
repMerge$mel.av <- (repMerge$mel.x + repMerge$mel.y + repMerge$mel) / 3
repMerge$sim.av <- (repMerge$sim.x + repMerge$sim.y + repMerge$sim) / 3
repMerge$sec.av <- (repMerge$sec.x + repMerge$sec.y + repMerge$sec) / 3

## subset just the major chromosomes
repMerge_mini <- repMerge[repMerge$CHROM %in% c("2L", "2R", "3L", "3R"),]
# rep1_mini <- rep1[rep1$CHROM %in% c("2L", "2R", "3L", "3R"),]

## strip the outliers
repMerge_mini<-repMerge_mini[!(abs(repMerge_mini$mel.av)>0.3),]
repMerge_mini<-repMerge_mini[!(abs(repMerge_mini$sim.av)>0.3),]
repMerge_mini<-repMerge_mini[!(abs(repMerge_mini$sec.av)>0.3),]

## make a masking frame that has the position for lines for Lhr and gfzf
dummy <- data.frame(CHROM = c('2R', '3R'), POS = c(17429032, 7145880))

## for plotting the graph
ggplot(data = repMerge_mini) +
  ## line plot data
  geom_line(aes(x=POS, 
                 y=mel.av),
             color='blue') +
  geom_line(aes(x=POS, 
                 y=sim.av),
             color='orange') +
  geom_line(aes(x=POS, 
                 y=sec.av),
             color='red') +  
  
  ## axis and title labels
  labs(title = 'Recomb Allele Frequencies',
       x="Position",
       y="Allele Frequency") +
  
  coord_cartesian(ylim=c(-0.28, 0.28)) +
  scale_y_continuous(breaks=seq(-0.5, 0.5, 0.05)) +
  facet_wrap(~CHROM, scales = 'free', nrow = 1, ncol = 4) +
  geom_vline(data=dummy, aes(xintercept=POS), color='black', linetype='dotted') +
  theme_void()

## for plotting the graph
ggplot(data = repMerge_mini) +
  ## line plot data
  geom_point(aes(x=POS, 
                 y=mel.av),
             color='grey') +
  geom_point(aes(x=POS, 
                 y=sim.av),
             color='orange') +
  
  ## axis and title labels
  labs(title = 'Recomb Allele Frequencies',
       x="Position",
       y="Allele Frequency") +
  
  coord_cartesian(ylim=c(-0.4, 0.4)) +
  scale_y_continuous(breaks=seq(0, 10, 0.25)) +
  facet_wrap(~CHROM, scales = 'free', nrow = 2, ncol = 2) +
  geom_vline(data=dummy, aes(xintercept=POS), color='black', linetype='dotted')


# , xlim=c(7500000, 11500000)
## for plotting the graph
ggplot(data = rep3) +
  ## line plot data
  geom_point(aes(x=POS, 
                 y=mel),
             color='blue') +
  geom_point(aes(x=POS, 
                 y=sim),
             color='orange') +
  geom_point(aes(x=POS, 
                 y=sec),
             color='red') +  
  
  ## axis and title labels
  labs(title = 'Recomb Allele Frequencies',
       x="Position",
       y="Allele Frequency") +
  
  coord_cartesian(ylim=c(-0.4, 0.4)) +
  scale_y_continuous(breaks=seq(0, 10, 0.25)) +
  facet_wrap(~CHROM, scales = 'free') +
  geom_vline(data=dummy, aes(xintercept=POS), color='black', linetype='dotted')


# geom_vline(xintercept=17429032, color='black', linetype="dotted") +
#, xlim=c(7500000, 11500000)
# plot2f = plot2f + geom_vline(xintercept=17429032, color='black')
# plot4f = plot4f + geom_vline(xintercept=7145880, color='black')

## for plotting the graph
ggplot(data = rep1) +
  ## line plot data
  geom_point(data = rep1, aes(x=POS, 
                 y=sim),
             color='blue') +
  geom_point(data = rep2, aes(x=POS, 
                 y=sim),
             color='orange') +
  geom_point(data = rep3, aes(x=POS, 
                 y=sim),
             color='red') +  
  
  ## axis and title labels
  labs(title = 'Recomb Allele Frequencies',
       x="Position",
       y="Allele Frequency") +
  
  coord_cartesian(ylim=c(-0.4, 0.4)) +
  scale_y_continuous(breaks=seq(0, 10, 0.25)) +
  facet_wrap(~CHROM, scales = 'free') +
  geom_vline(data=dummy, aes(xintercept=POS), color='black', linetype='dotted')
