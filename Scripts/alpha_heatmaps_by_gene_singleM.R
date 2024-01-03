library(gsveasyr)

env = read.csv('Data/mag_env_for_shotgun_samples.csv')
load_all_gsv_output('Data/urea/', min_num_reads = 1000000,
                    changeSampleName = TRUE, env_df = env,
                    refColumn = env$gsveasy_sample,
                    clustlvls = c(55, 60, 65, 70, 75, 80, 85, 90, 95, 100))

genes = unique(GSV_rarefaction_list$C[-nrow(GSV_rarefaction_list$C), ]$Gene)
metrics = c('C_sobs', 'C_shannon', 'C_chao1')
genes_tables = list()
metrics_tables = list()

for (j in metrics) {
  for (i in genes) {
    gene_table = cbind(
      GSV_alpha_div_cluster55[[j]][i],
      GSV_alpha_div_cluster60[[j]][i],
      GSV_alpha_div_cluster65[[j]][i],
      GSV_alpha_div_cluster70[[j]][i],
      GSV_alpha_div_cluster75[[j]][i],
      GSV_alpha_div_cluster80[[j]][i],
      GSV_alpha_div_cluster85[[j]][i],
      GSV_alpha_div_cluster90[[j]][i],
      GSV_alpha_div_cluster95[[j]][i],
      GSV_alpha_div_cluster100[[j]][i]
    )
    colnames(gene_table) = c('c55',
                             'c60',
                             'c65',
                             'c70',
                             'c75',
                             'c80',
                             'c85',
                             'c90',
                             'c95',
                             'c100')
    genes_tables[[i]] = gene_table
  }
  metrics_tables[[j]] = genes_tables
}

heatmaps = list()
heatmaps_genes = list()

for (j in metrics) {
for (i in genes) {
  heatmap_table = metrics_tables[[j]][[i]]
  env = env[order(env$pH), ]
  heatmap_table = heatmap_table[rownames(env),]
  
  col_fun = colorRampPalette(c('white', 'yellow', 'brown'))
  hm = Heatmap(t(heatmap_table), cluster_rows = F, cluster_columns = F, 
          name = j, column_title = i, col = col_fun(100))
  heatmaps_genes[[i]] = hm
}
  heatmaps[[j]] = heatmaps_genes
}

pdf('Figures/alpha_plots_by_gene_singleM.pdf')
heatmaps
dev.off()
