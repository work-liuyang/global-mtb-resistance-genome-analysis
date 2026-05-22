#!/usr/bin/Rscript
library(ggplot2)
library(dplyr)
library(maps)
library(grid)

colors <- c('#01539E', '#F4805D', '#CB6667', '#9BD5E7', '#CDB2AD', '#1A71B3', '#DEE89A', '#23A4D8', '#FF94C6', '#6AFFB7', '#D27FFF')
index <- c('Lineage1', 'Lineage2', 'Lineage3', 'Lineage4', 'Lineage5', 'Lineage6', 'Lineage7', 'La1', 'La2', 'La3', 'Ambiguous')


lag <- read.table("lineage_stats.xls",sep = "\t",header = T)
names(lag) <- c("Area", "lineage", "count")
lag <- arrange(lag,Area,lineage)

colorDF <-  data.frame(lineage=index, color=colors, row.names=1)

world_map<- map_data("world") |> subset(lat > -60)

gg <- ggplot() +
  geom_map(
    data = world_map, map = world_map,
    aes(long, lat, map_id = region),
    color = "white", fill = "lightgray", linewidth = 0.5) +
  theme(panel.background = element_blank(),
    panel.border = element_blank(),
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    plot.margin = unit(c(1, 5, 1, 1), "cm")
  )


size <- arrange(aggregate(lag$count,by=list(lag$Area),sum),desc(x))
names(size) <- c("area","size")

plinshi <- ggplot(size,aes(area,size)) + 
  geom_point(aes(size=size)) + 
  scale_size_continuous(range = c(2,24),
    breaks = c(10,500,1000,1500,2000))

size[,'plotsize'] <- (ggplot_build(plinshi)$data)[[1]]['size']
rm(plinshi)

for (i in 1:dim(size)[1]) {
  submap<-world_map[world_map$region==size$area[i],]
  long<-sum(submap$long)/length(submap$long)
  lat<-sum(submap$lat)/length(submap$lat)
  size$long[i] <- long
  size$lat[i] <- lat
}


for (x in 1:nrow(size)) {
  sublineages <- lag[lag$Area == size$area[x], ]$lineage
  sub_plot <- ggplot(lag[lag$Area == size$area[x], ], aes(x = 1, y = count, fill = lineage)) +
    geom_bar(width = 1, stat = "identity") +
    coord_polar("y", start = 0) +
    theme_void() +
    theme(legend.position = "None") +
    xlim(0, 1.5) +
    scale_fill_manual(values=colorDF[sublineages, "color"])
  
  sub_grob <- ggplotGrob(sub_plot)

  xmin <- max(size$long[x] - (size$plotsize[x]/2), -Inf)
  xmax <- min(size$long[x] + (size$plotsize[x]/2), Inf)
  ymin <- max(size$lat[x] - (size$plotsize[x]/2), -Inf)
  ymax <- min(size$lat[x] + (size$plotsize[x]/2), Inf)

  gg <- gg + annotation_custom(sub_grob,
                               xmin = xmin,
                               xmax = xmax,
                               ymin = ymin,
                               ymax = ymax)
  
  gg <- gg + annotate("text",
                      x = size$long[x]-5,
                      y = size$lat[x]-6,
                      label = paste(size$size[x], "\n", size$area[x]),
                      hjust = 0.03, size=3, alpha=0.7,
                      vjust = 0)
}


gg <- gg + 
  geom_bar(data = lag, 
           aes(x = 0, fill = factor(lineage, levels=index), y = 0), stat = "identity") +
  scale_fill_manual(name = "Lineages", values = colorDF$color)+
  theme(legend.position = c(1.02, 0.25), 
       legend.text = element_text(size=15),
       legend.title = element_text(size = 16))


gg <- gg + 
  geom_point(data=size, aes(x=0,y=0,size=size), pch=21) +
  scale_size_continuous(range = c(2,24),breaks = c(10,500,1000,1500,2000),labels = c('10', '500','1000','1500','>=2000')) + 
  labs(size='Number of isolates') +
  guides(size = guide_legend(order = 1), fill = guide_legend(order = 2)) +
  theme(legend.position = c(1.02, 0.55),
       legend.text = element_text(size=15),
       legend.title = element_text(size = 16))


ggsave('lineage_piemap.pdf', gg, height=10, width=18)

if(file_test("-f", "Rplots.pdf")) file.remove("Rplots.pdf")
