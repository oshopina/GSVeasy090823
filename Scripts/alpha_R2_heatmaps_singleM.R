library(gsveasyr)

env = read.csv('Data/mag_env_for_shotgun_samples.csv')
load_all_gsv_output('Data/urea/', min_num_reads = 1000000,
                    changeSampleName = TRUE, env_df = env,
                    refColumn = env$gsveasy_sample,
                    clustlvls = c(55, 60, 65, 70, 75, 80, 85, 90, 95, 100))

metrics = c('C_sobs', 'C_shannon', 'C_chao1')
alpha_results = list()
alpha_clusters = list(GSV_alpha_div_cluster55,
                      GSV_alpha_div_cluster60,
                      GSV_alpha_div_cluster65,
                      GSV_alpha_div_cluster70,
                      GSV_alpha_div_cluster75,
                      GSV_alpha_div_cluster80,
                      GSV_alpha_div_cluster85,
                      GSV_alpha_div_cluster90,
                      GSV_alpha_div_cluster95,
                      GSV_alpha_div_cluster100)
names(alpha_clusters) = c('c55', 'c60', 'c65', 'c70', 'c75', 'c80', 'c85', 'c90', 'c95', 'c100')

for (v in metrics) {
  alpha_r2 = data.frame()
  for (j in names(alpha_clusters)) {
    anova = auto_aov_fixed(alpha_clusters[[j]][[v]], ~ pH, alpha_clusters[[j]]$C.env)
    signif = subset(anova$Results, Signif != " ", select = Data)
    gene_names = as.character(signif$Data)
    
    r2s = c()
    for (i in gene_names) {
      r2 = summary(lm(
        alpha_clusters[[j]][[v]][[i]] ~ pH,
        alpha_clusters[[j]]$C.env
      ))$r.squared
      r2s[i] = r2
    }
    alpha_r2 = rbind(alpha_r2, cbind(
      Cluster = j,
      R2 = r2s,
      Gene = names(r2s)
    ))
  }
  alpha_results[[v]] = alpha_r2
}

heatmap_tables = list()

for (i in metrics) {
  heatmap_table = pivot_wider(alpha_results[[i]], names_from = Cluster, values_from = R2) %>% as.data.frame()
  rownames(heatmap_table) = heatmap_table$Gene
  heatmap_table = heatmap_table[, -1]
  heatmap_table[] = lapply(heatmap_table, as.numeric)
  heatmap_tables[[i]] = heatmap_table
}

library(ComplexHeatmap)

col_fun = colorRampPalette(c('white', 'yellow', 'brown'))
order_genes = rownames(heatmap_tables$C_sobs)

pdf('Figures/alpha_R2_heatmaps_singleM.pdf')
Heatmap(heatmap_tables$C_sobs, cluster_rows = F, cluster_columns = F, row_order = order_genes,
        name = 'R2', column_title = 'Sobs', col = col_fun(100))
Heatmap(heatmap_tables$C_shannon, cluster_rows = F, cluster_columns = F, row_order = order_genes,
        name = 'R2', column_title = 'Shannon', col = col_fun(100))
Heatmap(heatmap_tables$C_chao1, cluster_rows = F, cluster_columns = F, row_order = order_genes,
        name = 'R2', column_title = 'Chao1', col = col_fun(100))
dev.off()
        