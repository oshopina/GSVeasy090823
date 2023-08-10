library(gsveasyr)
library(vegan)
library(dplyr)

# Read the environmental data
env.all <- read.table('Data/mag_env_for_shotgun_samples.csv', header=TRUE, sep = ",")
rownames(env.all) <- env.all$SampleID

N_folder <- 'Data/N_cycle/'
urea_folder = 'Data/urea/'

#outliers "15_SD7466_H059_pH7" "35_SD7486_H082_pH6" "92_SD7543_S145"

# Load GSV data and process
load_graftm_gene_count(urea_folder, N_folder,
                    min_num_reads = 1000000,
                    changeSampleName = TRUE, env.df = env.all,
                    refColumn = env.all$gsveasy_sample)

names(GSV_gene_count_list)[1] = "urea_raw"
names(GSV_gene_count_list)[5] = "urea_per_million"

# Prepare pH group label as a factor
GSV_gene_count_list$env$pH.group.label <- as.factor(GSV_gene_count_list$env$pH.group.label)

rm(N_folder, urea_folder)

######################### Alpha diversity #####################################
alphas <- list()
for (i in c('N_per_million', 'urea_per_million')) {
  abundance_table <- GSV_gene_count_list[[i]]
  shannon <- diversity(abundance_table, index = "shannon")
  sobs <- rowSums(decostand(abundance_table, method = 'pa'))
  chao <- estimateR(round(abundance_table))[2,]
  alpha <- data.frame(shannon, sobs, chao, GSV_gene_count_list$env$pH, samples = names(shannon))
  colnames(alpha) <- c('shannon', 'sobs', 'chao', 'pH', 'samples')
  alphas[[i]] <- alpha
}

rm(abundance_table, alpha, chao, i, shannon, sobs)

library(ggplot2)
library(ggpmisc)

alpha_plots <- list()
for (i in c('N_per_million', 'urea_per_million')) {
  df <- alphas[[i]]
  df_filtered <- na.omit(df)
  plots <- list()
  for (j in c('shannon', 'sobs', 'chao')) {
    model <- lm(paste(j, "~ pH"), data = df)
    p_value <- summary(model)$coefficients["pH", "Pr(>|t|)"]
    color <- ifelse(coef(model)["pH"] > 0, "red", "blue")
    if (p_value > 0.05) {
      color <- "gray"
    }
    # Create the scatter plot with trend line using the lm method
    plot <- ggplot(df_filtered, aes(x = pH, y = !!sym(j))) +
      geom_point() +
      geom_smooth(method = 'loess', color = color) +
      stat_poly_eq(use_label(c("R2", "p"))) +
      geom_point() +
      ggtitle(i)
    plots[[j]] <- plot
  }
  alpha_plots[[i]] <- plots
}
rm(df, df_filtered, model, plot, plots, color, i, j, p_value)

# pdf("Figures/urea_alpha_plots.pdf")
# for (plot in alpha_plots$urea_per_million) {
#   print(plot)
# }
# dev.off()

## Gene abundance 
gene_plots <- list()
anova_results = list()
for (i in c('N_per_million', 'urea_per_million')) {
  anova <- auto_aov_fixed(GSV_gene_count_list[[i]], ~ pH, GSV_gene_count_list$env)
  anova_results[[i]] = anova$Results
  genes_list <- unique(anova$Results$Data)
  plots <- list()
  for (j in genes_list) {
    df <- cbind(GSV_gene_count_list[[i]][[j]], GSV_gene_count_list$env$pH) %>% as.data.frame()
    model <- lm(V2 ~ V1, data = df)
    p_value <- summary(model)$coefficients["V1", "Pr(>|t|)"]
    color <- ifelse(coef(model)["V1"] > 0, "red", "blue")
    if (p_value > 0.05) {
      color <- "gray"
    }
    # Create the scatter plot with trend line using the lm method
    plot <- ggplot(df, aes(x = V2, y = V1)) +
      geom_point() +
      geom_smooth(method = "loess", color = color) +
      stat_poly_eq(use_label(c("R2", "p"))) +
      xlab('pH') +
      ylab(j) +
      ggtitle(i)
    plots[[j]] <- plot
  }
  gene_plots[[i]] <- plots
}

rm(anova, df, model, plot, plots, color, genes_list, i, j, p_value)

# pdf("Figures/urea_gene_plots.pdf")
# for (plot in gene_plots$urea_per_million) {
#   print(plot)
# }
# dev.off()
