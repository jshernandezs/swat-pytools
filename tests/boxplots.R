# Script for making grouped Boxplots
rm(list = ls())

library(ggplot2)
library(dplyr)
library(ggpubr)
library(gridExtra)

theme_set(theme_pubclean())

input_dir <- '../figures'
opt_val <- T

if (!opt_val){
  files <- c('IHA1', 'IHA2', 'IHA3_4_5','MAG7')
  ylim <- c(-80, 80)
} else {
  files <- c('IHA1', 'IHA2', 'IHA3_4_5')
  ylim <- c(-80, 200)
}

plots <- plots_val <- list()
i <- 1
for (f in files){
  file <- paste0(input_dir,'/',f,'.csv')
  data <- read.csv(file)
  data$Strategy <- as.factor(data$Strategy)
  data$variable <- toupper(data$variable)
  
  file <- paste0(input_dir,'/',f,'_val.csv')
  data_val <- read.csv(file)
  data_val$Strategy <- as.factor(data_val$Strategy)
  data_val$variable <- toupper(data_val$variable)
  
  e <- ggplot(data, aes(x = variable, y = value)) + 
    geom_hline(yintercept=30, linetype="dashed", color = "black") +
    geom_hline(yintercept=-30, linetype="dashed", color = "black") +
    stat_boxplot(geom ='errorbar', aes(col = Strategy), 
                 position = position_dodge(0.8),
                 width = 0.4) + 
    geom_boxplot(aes(fill = Strategy),
                 position = position_dodge(0.8),
                 width = 0.5, size = 0.4) +
    coord_cartesian(ylim = ylim) +
    scale_fill_manual(values=c('lightblue', 'salmon')) + 
    scale_color_manual(values = c("black", "black")) +
    labs(y = '', x = '', subtitle = 'Calibration') +
    theme(legend.position = 'none',
          plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 1, face = 'italic'),
          axis.text.x = element_blank())
  
  e2 <- ggplot(data_val, aes(x = variable, y = value)) + 
    geom_hline(yintercept=30, linetype="dashed", color = "black") +
    geom_hline(yintercept=-30, linetype="dashed", color = "black") +
    stat_boxplot(geom ='errorbar', aes(col = Strategy), 
                 position = position_dodge(0.8),
                 width = 0.4) + 
    geom_boxplot(aes(fill = Strategy),
                 position = position_dodge(0.8),
                 width = 0.5, size = 0.4) +
    coord_cartesian(ylim = ylim) +
    scale_fill_manual(name = expression(italic("Strategy")), values=c('lightblue', 'salmon')) + 
    scale_color_manual(name = expression(italic("Strategy")), values = c("black", "black")) +
    labs( y = '', x = '', subtitle = 'Validation') +
    theme(legend.position = 'none',
          plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 1, face = 'italic'),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
  
  if (i == 2 || i == 4){
    e <- e + theme(axis.text.y=element_blank(),
                   axis.ticks.y=element_blank())
    
    e2 <- e2 + theme(axis.text.y=element_blank(),
                     axis.ticks.y=element_blank())
  }
  
  plots[[f]] <- e
  plots_val[[f]] <- e2
  i <- i + 1
}

if (!opt_val){
  f3 <- ggarrange(plots[[3]], plots[[4]],
                  labels = c('c)', 'd)'),
                  font.label = list(face='bold'),
                  ncol = 2, nrow = 1)
  
  f4 <- ggarrange(plots_val[[3]], plots_val[[4]],
                  labels = c('', ''),
                  font.label = list(face='plain'),
                  ncol = 2, nrow = 1,
                  common.legend = TRUE, legend = "bottom")
} else {
  f3 <- ggarrange(plots[[3]],
                  labels = c('c)', 'd)'),
                  font.label = list(face='bold'),
                  ncol = 1, nrow = 1)
  
  f4 <- ggarrange(plots_val[[3]],
                  labels = c('', ''),
                  font.label = list(face='plain'),
                  ncol = 1, nrow = 1,
                  common.legend = TRUE, legend = "bottom")
}

f1 <- ggarrange(plots[[1]], plots[[2]],
                labels = c('a)', 'b)'),
                font.label = list(face='bold'),
                ncol = 2, nrow = 1)

f2 <- ggarrange(plots_val[[1]], plots_val[[2]],
                labels = c('', ''),
                font.label = list(face='plain'),
                ncol = 2, nrow = 1)

fa <- ggarrange(f1, f2,
                ncol = 1, nrow = 2,
                heights = c(1, 1.13))
fa <- annotate_figure(fa, left = text_grob("Relative error (%)", vjust = 1, color = "black", rot = 90))

fb <- ggarrange(f3, f4,
                ncol = 1, nrow = 2,
                heights = c(1, 1.32))
fb <- annotate_figure(fb, left = text_grob("Relative error (%)", vjust = 1, color = "black", rot = 90))

if (!opt_val){
  figure <- ggarrange(fa, fb,
                      ncol = 1, nrow = 2,
                      heights = c(1, 1.1))
} else {
  layout_matrix <- matrix(c(1, 1, 1, 1, 3, 2, 2, 3), nrow = 2, byrow = TRUE)
  figure <- grid.arrange(fa, fb, layout_matrix = layout_matrix, heights = c(1, 1.1))
}

if (!opt_val){
  ggsave(paste0(input_dir,'/ERHI_perf.tiff'), plot = figure, dpi = 600, width = 12, height = 14, units = "cm", scale = 2, 
         device = 'tiff', compression = "lzw")
} else {
  ggsave(paste0(input_dir,'/ERHI_var.tiff'), plot = figure, dpi = 600, width = 12, height = 14, units = "cm", scale = 2, 
         device = 'tiff', compression = "lzw")
}