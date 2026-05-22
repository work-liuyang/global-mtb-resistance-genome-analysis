library(ggplot2)

snparisedir <- "./table"
rifsnp <- read.table(paste0(snparisedir, "/rifampicin_snp_change.lineplot.xls"),
		                          sep = "\t", header = T,quote = "",check.names = F)

inhsnp <- read.table(paste0(snparisedir, "/isoniazid_snp_change.lineplot.xls"),
		                          sep = "\t", header = T,quote = "",check.names = F)

flqsnp <- read.table(paste0(snparisedir, "/fluoroquinolones_snp_change.lineplot.xls"),
		                          sep = "\t", header = T,quote = "",check.names = F)

bedsnp <- read.table(paste0(snparisedir, "/bedaquiline_snp_change.lineplot.xls"),
		                          sep = "\t", header = T,quote = "",check.names = F)

linsnp <- read.table(paste0(snparisedir, "/linezolid_snp_change.lineplot.xls"),
		                          sep = "\t", header = T,quote = "",check.names = F)


snplineplot <- function(x, prefix, width, height){
  
  psnp <- ggplot(x, aes(x = year, y = snpnum, group = country)) + 
    geom_point(show.legend = F,color='#65B8E1',size=0.5) +
    geom_line(show.legend = F,color='#65B8E1') +
    theme(axis.text.x = element_text(angle = 45,hjust = 1)) +
    theme(axis.text = element_text(size = 9)) +
    theme(axis.title = element_text(size=9,face = "bold")) + 
    theme(strip.text = element_text(size=9)) +
    theme(strip.background = element_rect(fill = '#65B8E1')) +
    labs(x='Years',y='The cumulative number of distinct drug resistance-associated SNPs',title = prefix) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size=9)) +
    scale_x_continuous(breaks=c(1994,2000,2005,2010,2015,2020,2025),
                       labels = c(1994,2000,2005,2010,2015,2020,2025))
  if(max(na.omit(x$snpnum))<=3){
    psnp <- psnp + scale_y_continuous(breaks = seq(0, max(na.omit(x$snpnum)), by = 1))
  }

  psnp <- psnp + facet_wrap(~country,scales = 'fixed')

  ggsave(filename = paste0(prefix,".snp.rise.pdf"),plot = psnp,width = width, height = height)
}

snplineplot(rifsnp,'Rifampicin',9,9)
snplineplot(inhsnp,'Isoniazid',9,9)
snplineplot(flqsnp,'Fluoroquinolones',9,9)
snplineplot(bedsnp,'Bedaquiline',9,9)
snplineplot(linsnp,'Linezolid',9,9)
