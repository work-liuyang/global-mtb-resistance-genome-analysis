
#library(tidyverse)
library(dplyr)
library(readxl)
library(ggpubr)
library(patchwork)
library(stringr)

file_plot <- "fold_metrics.xlsx"   
file_test <- "t_test_results.xlsx"


sheets_plot <- excel_sheets(file_plot)
data_list <- lapply(sheets_plot, function(s){
  read_excel(file_plot, sheet = s)
})
names(data_list) <- sheets_plot

sheets_test <- excel_sheets(file_test)
test_list <- lapply(sheets_test, function(s){
  read_excel(file_test, sheet = s)
})
names(test_list) <- sheets_test

center_dot <- function(x) {
  str_replace(sprintf("%.2f", x), "\\.", "·")
}

plot_one_sheet <- function(data, test_df, sheet_name){
  
  data$strategy <- factor(data$strategy, levels=unique(data$strategy))
  metrics <- c("AUPRC", "F1", "Precision", "Recall")
  plot_list <- list()
  
  for(m in metrics){
    
    sub_test <- test_df %>% filter(metric == !!m)
    print(sub_test)
 
    p <- ggboxplot(
      data,
      x = "strategy",
      y = m,
      color = "strategy",
      apha = 0.6,
      palette = "lancet",
      bxp.errorbar = TRUE,
      bxp.errorbar.width = 0.2,
      x.text.angle = 45,
      add = "jitter",
      add.params = list(shape=1)
    ) +
      ggtitle(toupper(m)) +
      theme_bw() +
      theme(
        axis.text.x = element_text(face = "bold",size = 9),
        axis.text.y = element_text(size = 9),
        
        panel.grid = element_blank(),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5, size=9, face = "bold")
      ) +
      labs(x=NULL,y=NULL) +
      scale_y_continuous(labels = center_dot)#, limits = c(NA,1))

   

    if(nrow(sub_test) > 0){
      yposition = max(data[[m]], na.rm=TRUE) * seq(1.01, 1.1, length.out=as.numeric((table(sub_test$metric)[m])))
      if(s == "Rifampicin" && m == "Precision"){
        yposition = max(data[[m]], na.rm=TRUE) * seq(1.01, 1.015, length.out=as.numeric((table(sub_test$metric)[m])))
      }
      if(s == "Ethambutol" && m == "Precision"){
        yposition = max(data[[m]], na.rm=TRUE) * seq(1.01, 1.02, length.out=as.numeric((table(sub_test$metric)[m])))
      }
      if(s == "Isoniazid" && m == "Precision"){
        yposition = max(data[[m]], na.rm=TRUE) * seq(1.01, 1.05, length.out=as.numeric((table(sub_test$metric)[m])))
      }


      p <- p + stat_pvalue_manual(
        hide.ns = F,
        sub_test,
        label = "significance",
#        y.position = max(data[[m]], na.rm=TRUE) * seq(1.01, 1.1, length.out=as.numeric((table(sub_test$metric)[m])))
        y.position = yposition
      )
    }
    
    plot_list[[m]] <- p
  }
  
  combined_plot <- (plot_list$AUPRC | plot_list$F1) /
                   (plot_list$Precision | plot_list$Recall) +
    plot_annotation(
      title = paste0("Drug: ", sheet_name),
      theme = theme(plot.title = element_text(hjust = 0.5, size = 9, face = "bold"))
    )
  
  return(combined_plot) 
}

for(s in sheets_plot){
  cat("Plot：", s, "\n")
  p <- plot_one_sheet(data_list[[s]], test_list[[s]], s)
  ggsave(paste0("Model_Compare_", s, ".pdf"), p, width=14, height=10)
  ggsave(paste0("Model_Compare_", s, ".png"), p, width=14, height=10, dpi=300)
}

