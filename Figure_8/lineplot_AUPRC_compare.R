library(ggplot2)
library(ggrepel)
library(dplyr)
library(tidyr)


drug_order <- c("Rifampicin", "Isoniazid", "Ethambutol", "Ethionamide", "Amikacin", "Levofloxacin", "Moxifloxacin", "Kanamycin")
color1 <- c("#EA8379","#7DAEE0","#B395BD")
color2 <- c("#F16C23", "#2B6A99", "#1B7C3D")


all_data <- read.delim("cv.xls")
all_data$drug <- factor(all_data$drug, levels = drug_order)

p1 <- ggplot(all_data, aes(x = drug, y = val, color = group, group = group)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2.5) +
  geom_text_repel(
    aes(label = round(val, 3)),
    size = 3,
    seed = 123,
    force = 2,
    max.overlaps = 100
  ) +
  scale_color_manual(
    breaks = c("AssowR_and_Uncertain", "AssowR", "AssowR_predominant"),
    labels = c("AssowR_and_Uncertain", "AssowR", "AssowR_predominant"),
    values = color2
    ) +
  ylim(-0.05, 1.05) +
  labs(x = "Drug", y = "Area Under the Precision-Recall Curve", color="Model") +
  theme_bw() +
  theme(
    legend.position = "right",
    panel.grid.major.y = element_line(linetype = "dashed"), 
    panel.grid.minor.y = element_line(linetype = "dashed"),
    panel.grid.major.x = element_blank()
  )

p2 <- ggplot(all_data, aes(x = drug, y = mean, color = group, group = group)) +
  geom_errorbar(
    aes(ymin = mean - std, ymax = mean + std),
    width = 0.1,
    linewidth = 0.5
  ) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.5) +
  geom_text_repel(
    aes(label = round(mean, 3)),
    size = 3,
    seed = 123,
    force = 2,
    max.overlaps = 100
  ) +
  scale_color_manual(
    breaks = c("AssowR_and_Uncertain", "AssowR", "AssowR_predominant"),
    labels = c("AssowR_and_Uncertain", "AssowR", "AssowR_predominant"),
    values = color2
    ) +
  ylim(-0.05, 1.05) +
  labs(x = "Drug", y = "Area Under the Precision-Recall Curve", color="Model") +
  theme_bw() +
  theme(
    legend.position = "right",
    panel.grid.major.y = element_line(linetype = "dashed"),
    panel.grid.minor.y = element_line(linetype = "dashed"),
    panel.grid.major.x = element_blank()
  )

ggsave("AUPRC_val_plot.pdf", p1, width = 10, height = 6, device = "pdf")
ggsave("AUPRC_mean_error_plot.pdf", p2, width = 10, height = 6, device = "pdf")

