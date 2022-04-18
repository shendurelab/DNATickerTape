library(tidyverse)

P1E16_bigram12 <- read.delim2('TAPE-Prog1-Epoch16-repAll.csv', sep = ',', strip.white = TRUE)
P1E16_bigram12$ReadCounts <- as.numeric(as.character(P1E16_bigram12$ReadCounts))

P1E16_bigram12$Position1 <- factor(P1E16_bigram12$Position1, levels = rev(c('CA','GC','TA','GT',
                                                                            'CG','AT','TG','AC',
                                                                            'GG','CT','AA','TT',
                                                                            'AG','CC','GA','TC')))

P1E16_bigram12$Position2 <- factor(P1E16_bigram12$Position2, levels = c('CA','GC','TA','GT',
                                                                        'CG','AT','TG','AC',
                                                                        'GG','CT','AA','TT',
                                                                        'AG','CC','GA','TC'))


P1E16_bigram12_RowColnorm <- P1E16_bigram12 %>%
  group_by(Position2) %>%
  mutate(ReadCountsColNorm = ReadCounts / sum(ReadCounts)) %>%
  group_by(Position1) %>%
  mutate(ReadCountsRowColNorm = ReadCountsColNorm / sum(ReadCountsColNorm))

P1E16_bigram12_RowColnorm <- filter(P1E16_bigram12_RowColnorm, Position1 != Position2)
P1E16_bigram12_RowColnorm <- mutate(P1E16_bigram12_RowColnorm, ReadCountsRowColNorm = replace(ReadCountsRowColNorm, ReadCountsRowColNorm > 0.25, 0.25))


ggplot(P1E16_bigram12_RowColnorm, aes(y = Position1, x = Position2, fill= ReadCountsRowColNorm)) + 
  geom_tile()+
  scale_fill_distiller(palette = "YlGnBu") +
  theme(axis.text.x = element_text(angle = -90, vjust = 1, hjust=0.3, size = 12), 
        axis.text.y = element_text(size = 12),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "black",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        axis.line = element_line(size = 0.5, color = rgb(0,0,0,max=255)))
ggsave(file = 'Figure2C_Program1_2Dhistogram_plot.eps', height = 3, width = 5)



