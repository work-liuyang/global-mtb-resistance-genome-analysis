library(ggplot2)
library(ggsci)

prevdir<-'./table'

inhprev <- read.table(paste0(prevdir, '/isoniazid.prevalence.stats.xls'), sep = '\t', 
                      header = T, check.names = F,quote = "")
rifprev <- read.table(paste0(prevdir, '/rifampicin.prevalence.stats.xls'), sep = '\t', 
                      header = T, check.names = F,quote = "")
flqprev <- read.table(paste0(prevdir, '/fluoroquinolones.prevalence.stats.xls'), sep = '\t', 
                      header = T, check.names = F,quote = "")
bedprev <- read.table(paste0(prevdir, '/bedaquiline.prevalence.stats.xls'), sep = '\t', 
                      header = T, check.names = F,quote = "")
linprev <- read.table(paste0(prevdir, '/linezolid.prevalence.stats.xls'), sep = '\t', 
                      header = T, check.names = F,quote = "")

allprev <- rbind(
  data.frame(drug=rep('Isoniazid',5),inhprev[1:5,1:2]),
  data.frame(drug=rep('Rifampicin',5),rifprev[1:5,1:2]),
  data.frame(drug=rep('Fluoroquinolones',5),flqprev[1:5,1:2]),
  data.frame(drug=rep('Bedaquiline',5),bedprev[1:5,1:2]),
  data.frame(drug=rep('Linezolid',5),linprev[1:5,1:2])
)

allprev <- na.omit(allprev)
allprev$Variants <- factor(allprev$Variants, levels = unique(allprev$Variants))
allprev$drug <- factor(allprev$drug, levels = unique(allprev$drug))

pprev <- ggplot(allprev,aes(Variants,count,fill=drug)) + 
  geom_bar(stat = 'identity',position = 'dodge', width=0.8, alpha=0.6) + 
  scale_fill_npg() + 
  scale_x_discrete(limits=rev(allprev$Variants)) +
  scale_y_continuous(expand = c(0,0),limits = c(0,4000)) +
  coord_flip() +
  theme_bw() +
  theme(panel.border = element_blank()) +
  theme(legend.title=element_blank(), 
        legend.text=element_text(angle=270, size=10)) +
  theme(legend.key.width=unit(2, 'mm'), 
        legend.key.height = unit(5,'cm'),
        legend.text.align=0.5) +
  theme(axis.title = element_text(size=12),
        axis.text = element_text(size=10)) +
  geom_text(aes(label = count),hjust=-0.5,size=3.5)

ggsave(filename = "variants.top5.pdf",plot = pprev,width = 6, height = 8)

