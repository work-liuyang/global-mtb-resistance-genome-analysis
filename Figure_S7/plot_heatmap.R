heatmap_data <- read.delim(
  "plot_heatmap.xls",
  sep = "\t",
  header = TRUE,row.names = 1,
  check.names = FALSE
)

library(pheatmap)

#assocwr <- heatmap_data[heatmap_data$Confidence=="AssocwR",]
# breaks<-c(0.000, 0.0001, 0.001, 0.005, 0.010, 0.020, 0.050, 0.100, 0.200, 0.500)
# color <- colorRampPalette(c("#f7fbff", "#6baed6", "#08519c"))(9)
# pheatmap(
#   assocwr[,-9],
#   cluster_rows = F,
#   cluster_cols = F, 
#   breaks = breaks,
#   color = color, 
#   border_color = "grey100", 
#   na_col = "#FFFFF7", 
#   angle_col = 0,
#   fontsize = 9,
#   filename = "AssocwR_features_heatmap.pdf", 
#   width = 10, height = 12
# )

#uncertain  <- heatmap_data[heatmap_data$Confidence=="Uncertain",]
breaks<-c(0.000, 0.00001, 0.0001, 0.001, 0.005, 0.010, 0.020, 0.050, 0.100, 0.200, 0.500)
color <- colorRampPalette(c("#f7fbff", "#6baed6", "#08519c"))(10)
# pheatmap(
#   uncertain[,-9],
#   cluster_rows = F,
#   cluster_cols = F, 
#   breaks = breaks,
#   color = color, 
#   border_color = "grey100", 
#   na_col = "#FFFFF7", 
#   angle_col = 0,
#   fontsize = 9,
#   filename = "Uncertain_features_heatmap.pdf", 
#   width = 10, height = 12
# )


pheatmap(
  heatmap_data[,-9],
  cluster_rows = F,
  cluster_cols = F, 
  breaks = breaks,
  color = color, 
  border_color = "grey100", 
  na_col = "#FFFFF7", 
  angle_col = 0,
  fontsize = 9,
  annotation_row = heatmap_data[,"Confidence",drop=F],
  annotation_colors = list(Confidence=c(AssocwR = "#ED0000", Uncertain="#00468B", Mutation_from_literature="#42B540")),
  filename = "All_features_heatmap.pdf", 
  width = 12, height = 16
)
