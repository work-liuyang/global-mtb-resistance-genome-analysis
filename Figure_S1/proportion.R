library(ggplot2)
library(patchwork)
library(stringr)

#propdir<-'C:/Users/刘杨/Desktop/04.DR_lineage_stats/proportion/'
propdir <- './'

center_dot <- function(x) {
  str_replace(sprintf("%.1f", x), "\\.", "·")
}

y_limits <- c(0, 1)
y_breaks <- seq(0, 1, 0.2)

#isoniazid
inhprop <- read.table(paste0(propdir, 'isoniazid.prop.xls'), sep = '\t', 
                      header = T, check.names = F,quote = "")
inhprop$country <- factor(inhprop$country, levels = unique(inhprop$country))
inhprop$proportion[inhprop$proportion==0] <- NA

pinh <- ggplot(inhprop,aes(country,proportion,fill=inhprop)) + 
  geom_bar(stat = 'identity',position = 'stack', fill='#3A89C1',alpha=0.8,width = 0.8)  + 
  theme_minimal() + labs(y='Isoniazid',x='Country') + 
  theme(panel.grid = element_blank(), axis.text=element_text(size=9), axis.title=element_text(size=9, face="bold")) +
  scale_y_continuous(expand = c(0,0), labels=center_dot, limits=y_limits, breaks=y_breaks)  + coord_flip()

#rifampicin
rifprop <- read.table(paste0(propdir, 'rifampicin.prop.xls'), sep = '\t', 
                      header = T, check.names = F,quote = "")
rifprop$country <- factor(rifprop$country, levels = unique(rifprop$country))
rifprop$proportion[rifprop$proportion==0] <- NA

prif <- ggplot(rifprop,aes(country,proportion,fill=rifprop)) + 
  geom_bar(stat = 'identity',position = 'stack', fill='#3A89C1',alpha=0.8,width = 0.8)  + 
  theme_minimal() + labs(y='Rifampicin',x='') + 
  theme(axis.text.y = element_blank()) +
  theme(panel.grid = element_blank(), axis.text=element_text(size=9), axis.title=element_text(size=9, face="bold")) +
  scale_y_continuous(expand = c(0,0), labels=center_dot, limits=y_limits, breaks=y_breaks) + coord_flip()


#fluoroquinolones
flqprop <- read.table(paste0(propdir, 'fluoroquinolones.prop.xls'), sep = '\t', 
                      header = T, check.names = F,quote = "")
flqprop$country <- factor(flqprop$country, levels = unique(flqprop$country))
flqprop$proportion[flqprop$proportion==0] <- NA

pflq <- ggplot(flqprop,aes(country,proportion,fill=flqprop)) + 
  geom_bar(stat = 'identity',position = 'stack', fill='#3A89C1',alpha=0.8,width = 0.8)  + 
  theme_minimal() + labs(y='Fluoroquinolones',x='') + 
  theme(axis.text.y = element_blank()) +
  theme(panel.grid = element_blank(), axis.text=element_text(size=9), axis.title=element_text(size=9, face="bold")) +
  scale_y_continuous(expand = c(0,0), labels=center_dot, limits=y_limits, breaks=y_breaks) + coord_flip()

#bedaquiline
bedprop <- read.table(paste0(propdir, 'bedaquiline.prop.xls'), sep = '\t', 
                      header = T, check.names = F,quote = "")
bedprop$country <- factor(bedprop$country, levels = unique(bedprop$country))
bedprop$proportion[bedprop$proportion==0] <- NA

pbed <- ggplot(bedprop,aes(country,proportion,fill=bedprop)) + 
  geom_bar(stat = 'identity',position = 'stack', fill='#3A89C1',alpha=0.8,width = 0.8)  + 
  theme_minimal() + labs(y='Bedaquiline',x='') + 
  theme(axis.text.y = element_blank()) +
  theme(panel.grid = element_blank(), axis.text=element_text(size=9), axis.title=element_text(size=9, face="bold")) +
  scale_y_continuous(expand = c(0,0), labels=center_dot, limits=y_limits, breaks=y_breaks) + coord_flip()

#linezolid
linprop <- read.table(paste0(propdir, 'linezolid.prop.xls'), sep = '\t', 
                      header = T, check.names = F,quote = "")
linprop$country <- factor(linprop$country, levels = unique(linprop$country))
linprop$proportion[linprop$proportion==0] <- NA

plin <- ggplot(linprop,aes(country,proportion,fill=linprop)) + 
  geom_bar(stat = 'identity',position = 'stack', fill='#3A89C1',alpha=0.8,width = 0.8)  + 
  theme_minimal() + labs(y='Linezolid',x='') + 
  theme(axis.text.y = element_blank()) +
  theme(panel.grid = element_blank(), axis.text=element_text(size=9), axis.title=element_text(size=9, face="bold")) +
  scale_y_continuous(expand = c(0,0), labels=center_dot, limits=y_limits, breaks=y_breaks) + coord_flip()

pprop <- pinh + prif + pflq + pbed + plin + plot_layout(nrow = 1)

ggsave(filename = paste0(propdir, "resistance.mutation.proportion.png"),plot = pprop,width = 12, height = 5.7,dpi = 300)
ggsave(filename = paste0(propdir, "resistance.mutation.proportion.pdf"),plot = pprop,width = 12, height = 5.7)
