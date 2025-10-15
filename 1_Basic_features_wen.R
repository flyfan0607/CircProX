#!/usr/bin/env Rscript
# =============================================================================
# circRNA基础特征可视化脚本
# 文件: 1_Basic_features.R
# 作者: 可视化生物信息工程师
# 日期: 2024年
# 描述: 对circRNA数据中的基础特征进行可视化分析
# 包括染色体分布、链方向、circHORF类型等
# =============================================================================

# 提高运行速度的设置
options(stringsAsFactors = FALSE)
options(future.globals.maxSize = 8000 * 1024^2)

# 加载必要的R包
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(scales)
  library(patchwork) # 用于组合多个图表
})

# 设置工作目录和创建输出目录
#setwd("f:/★文件/文章/1.circRNA-算法/2.circProX/Data/example and code/afterCICADA")
dir.create("Basic_features", showWarnings = FALSE, recursive = TRUE)

cat("=== Beginning! ===\n")
cat("Begin time:", as.character(Sys.time()), "\n\n")

# 读取数据
cat("Reading data...\n")
data <- read.csv("merged_results.csv", stringsAsFactors = FALSE)

# 1. 染色体分布可视化
# =============================================================================
if("Chromosome" %in% names(data)) {
  cat("1. Chromosome distribution\n")
  
  # 数据预处理
  chr_data <- data %>%
    count(Chromosome, sort = TRUE) %>% arrange(desc(n)) %>%
    mutate(
      Percentage = round(n / sum(n) * 100, 2),
      Chromosome = factor(Chromosome, levels = Chromosome)
    )
  
  # 创建水平条形图
  p1_horizontal <- ggplot(chr_data, aes(x = reorder(Chromosome, n), y = n, fill = n)) +
    geom_col(alpha = 0.8) +
    scale_fill_viridis_c(option = "turbo", name = "Number") +  # 使用viridis配色
    coord_flip() +
    labs(
      title = "Distribution of translatable circRNA' chromosome",
      #subtitle = "水平条形图展示",
      x = "Chromosome", 
      y = "circRNA number"
    ) +
    scale_y_discrete(expand = expansion(mult = 0.01))+
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 12),
      panel.grid = element_blank(),
      legend.position = "right",
      legend.margin = margin(r = 10, l = -50),
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 12)
    ) +
    geom_text(aes(label = paste0(n, " (", Percentage, "%)"), 
              hjust = -0.1), size = 3.5, show.legend = FALSE)+
    expand_limits(y = max(chr_data$n) * 1.15)

  # 保存高质量图片
  ggsave(
    "Basic_features/chromosome_horizontal_bar.pdf",  # 保存为PDF以保证质量
    p1_horizontal, 
    width = 10,
    height = 8,
    dpi = 300)
  
  cat("Chromosome distribution has been saved (PDF and PNG)\n")
}

# 2. 链方向分布可视化
# =============================================================================
if("Strand" %in% names(data)) {
  cat("\n2. Strand distribution\n")
  
  # 数据预处理
  strand_data <- data %>%
    count(Strand) %>%
    mutate(
      Percentage = round(n / sum(n) * 100, 2),
      Label = paste0(Strand, "\n", n, " (", Percentage, "%)"))
  
  p2_pie <- ggplot(strand_data, aes(x = "", y = n, fill = Strand)) +
    geom_col(width = 1, color = "white", linewidth = 1.5) +  
    coord_polar("y", start = 0) +
    scale_fill_viridis_d(option = "turbo", name = "链方向",alpha = 0.7) +
    labs(
      title = "Distribution of translatable circRNA' strand"
    ) +
    theme_void() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold", 
                                margin = margin(b = -30)),  
      legend.position = "none",
      plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")  
    ) +
    geom_text(aes(label = Label), 
              position = position_stack(vjust = 0.5), 
              size = 4,
              color = "white",  # 设置文字为白色，增强与背景的对比度
              fontface = "bold") +  # 文字加粗
    guides(fill = guide_legend(title.position = "top",  # 图例标题居上
                               title.hjust = 0.5,       # 图例标题居中
                               nrow = 1))  
  
  # 保存高质量图片
  ggsave(
    "Basic_features/strand_pie_chart.pdf",
    p2_pie,
    width = 8,
    height = 8,
    dpi = 300
  )

  
  cat("Strand distribution has been saved（PDF and PNG）\n")
}
# 3. 编码概率分布分析
# =============================================================================
if("Coding.probability.score" %in% names(data)) {
  cat("\n3. Coding probability score distribution\n")
  
  # 数据预处理和统计
  coding_prob <- as.numeric(data$Coding.probability.score)
  coding_stats <- summary(coding_prob)
  print(coding_stats)
  
  # 创建编码概率分布图
  p3_dist <- ggplot(data, aes(x = as.numeric(Coding.probability.score))) +
    geom_histogram(aes(y = after_stat(density)), bins = 30,
                   fill = "#4ECDC4", alpha = 0.7, color = "#45B7D1") +
    geom_density(color = "#FF6B6B", size = 1) +
    labs(
      title = "Distribution of translatble circRNAs' coding probability score",
      subtitle = paste("Mean:", round(mean(coding_prob, na.rm = TRUE), 2),
                      "Median:", round(median(coding_prob, na.rm = TRUE), 2)),
      x = "Coding probability score",
      y = "Density"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 12),
      panel.grid = element_blank(),
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 12)
    ) +
    scale_x_continuous(labels = scales::number_format(accuracy = 0.01)) +
    scale_y_continuous(labels = scales::number_format(accuracy = 0.01))
  
  # 保存高质量图片
  ggsave(
    "Basic_features/coding_probability_distribution.pdf",
    p3_dist,
    width = 10,
    height = 6,
    dpi = 300
  )
  
  cat("Coding probability distribution has been saved (PDF and PNG)\n")
}
# 4. 编码产物数量分布分析
# =============================================================================
if("Product.number" %in% names(data)) {
  cat("\n4. Coding product number distribution\n")
  
  # 数据预处理和统计
  product_num <- as.numeric(data$Product.number)
  product_stats <- summary(product_num)
  print(product_stats)
  
  # 创建编码产物数量分布图
  p4_dist <- ggplot(data, aes(x = as.numeric(Product.number))) +
    geom_histogram(aes(y = after_stat(density)), bins = 30,
                   fill = "#4ECDC4", alpha = 0.7, color = "#45B7D1") +
    geom_density(color = "#FF6B6B", size = 1) +
    labs(
      title = "Distrubition of translatble circRNAs' coding product number",
      subtitle = paste("Mean:", round(mean(product_num, na.rm = TRUE), 2),
                      "Median:", round(median(product_num, na.rm = TRUE), 2)),
      x = "Number of coding products",
      y = "Density"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 12),
      panel.grid = element_blank(),
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 12)
    ) +
    scale_x_continuous(labels = scales::number_format(accuracy = 1)) +
    scale_y_continuous(labels = scales::number_format(accuracy = 0.01))
  
  # 保存高质量图片
  ggsave(
    "Basic_features/product_number_distribution.pdf",
    p4_dist,
    width = 10,
    height = 6,
    dpi = 300
  )
  
  # 同时保存PNG格式用于快速预览
  ggsave(
    "Basic_features/product_number_distribution.png",
    p4_dist,
    width = 10,
    height = 6,
    dpi = 300
  )
  
  cat("Coding product number distribution has been saved (PDF and PNG)\n")
}
# 5. circHORF类型分布可视化
# =============================================================================
if("circHORF.type" %in% names(data)) {
  cat("\n5. circHORF type distribution\n")
  
  # 数据预处理
  type_data <- data %>%
    count(circHORF.type) %>%
    mutate(
      Percentage = round(n / sum(n) * 100, 2),
      Label = paste0(circHORF.type, "\n", n, " (", Percentage, "%)")
    )

  p5_rose <- ggplot(type_data, aes(x = circHORF.type, y = n, fill = circHORF.type)) +
    geom_col(width = 1, color = "white", size = 1) +
    coord_polar("x", start = 0, direction = -1) +
    # 使用指定的绿色、黄色、紫色配色，并设置饱和度为0.7
    scale_fill_manual(
      values = alpha(c("#4CAF50", "#FFC107", "#9C27B0"), alpha=0.7),  # 绿、黄、紫
      name = "circHORF类型"
    ) +
    scale_y_continuous(expand = c(0, 0)) +
    labs(
      title = "Distribution of circHORF type",
      x = NULL,
      y = NULL
    ) +
    theme_void() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold",
                                margin = margin(b = -30)),
      # 移除图例
      legend.position = "none",
      plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")
    ) +
    geom_text(aes(label = Label), 
              position = position_stack(vjust = 0.5), 
              size = 3.5,
              fontface = "bold")
  
  # 保存高质量图片
  ggsave(
    "Basic_features/circHORF_type_pie_chart.pdf",
    p5_rose,
    width = 8,
    height = 8,
    dpi = 300
  )
  
  # 同时保存PNG格式用于快速预览
  ggsave(
    "Basic_features/circHORF_type_pie_chart.png",
    p5_rose,
    width = 8,
    height = 8,
    dpi = 300
  )
  
  cat("circHORF type distribution has been saved (PDF and PNG)\n")
}
# 6. 合并所有图表为一张大图
# =============================================================================
cat("\n6. Combining all plots into one figure\n")

# 创建一个列表存储所有图表
plot_list <- list()

# 检查各个图表是否存在并添加到列表中
if(exists("p1_horizontal")) plot_list$p1 <- p1_horizontal
if(exists("p2_pie")) plot_list$p2 <- p2_pie
if(exists("p3_dist")) plot_list$p3 <- p5_rose
if(exists("p4_dist")) plot_list$p4 <- p4_dist
if(exists("p5_rose")) plot_list$p5 <- p3_dist

# 检查是否有图表可以合并
if(length(plot_list) > 0) {
  # 计算每行应该有多少图表
  n_plots <- length(plot_list)
  plots_per_row <- ceiling(n_plots / 2)
  
  # 使用patchwork包合并图表
  # 第一行图表
  row1 <- NULL
  for(i in 1:min(plots_per_row, n_plots)) {
    if(i == 1) {
      row1 <- plot_list[[i]]
    } else {
      row1 <- row1 + plot_list[[i]]
    }
  }
  
  # 第二行图表（如果有）
  row2 <- NULL
  if(n_plots > plots_per_row) {
    for(i in (plots_per_row+1):n_plots) {
      if(i == plots_per_row+1) {
        row2 <- plot_list[[i]]
      } else {
        row2 <- row2 + plot_list[[i]]
      }
    }
  }
  
  # 合并两行
  if(!is.null(row2)) {
    combined_plot <- row1 / row2
  } else {
    combined_plot <- row1
  }

  # 保存合并后的图表
  ggsave(
    "Basic_features/combined_basic_features.pdf",
    combined_plot,
    width = 18,
    height = 12,
    dpi = 300
  )
  cat("Combined figure has been saved (PDF and PNG)\n")
}

# 清理内存
rm(list = ls())
gc()
