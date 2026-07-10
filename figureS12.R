# Portable project path. Set TLS_PROJECT_ROOT to the directory containing the input data.
PROJECT_ROOT <- normalizePath(Sys.getenv("TLS_PROJECT_ROOT", unset = "."), winslash = "/", mustWork = FALSE)


####胃#####


library(igraph)

library(Matrix)
library(Biobase)
library(BiocGenerics)
library(VGAM)
library(stats4)
library(splines)
library(DDRTree)
library(irlba)
library(monocle)
library(Seurat)
library(Matrix)

library(sp)
library(SeuratObject)
library(SeuratDisk)
library(Seurat) 

library(harmony)
library(Seurat)
library(SeuratDisk)
library(rhdf5)
library(readr)

pbmc = readRDS(file.path(PROJECT_ROOT, "泛癌S", "h5ad", "探究癌前_癌", "胃", "wei_with_Normalwei_subset.rds"))


expr_matrix <- as(pbmc@assays$RNA@counts, "CsparseMatrix")
f_data <- data.frame(gene_short_name = rownames(pbmc), row.names = rownames(pbmc))
p_data <- pbmc@meta.data
p_data$celltype <- pbmc@active.ident

pd <- new("AnnotatedDataFrame", data = p_data)
fd <- new("AnnotatedDataFrame", data = f_data)

cds <- newCellDataSet(expr_matrix,
                      phenoData = pd,
                      featureData = fd,
                      lowerDetectionLimit = 0.5,
                      expressionFamily = negbinomial.size())

# 预处理
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
cds <- detectGenes(cds, min_expr = 0.1)

# 差异分析
disp_table = dispersionTable(cds)

disp.gene = subset(disp_table, mean_expression >= 0.1 & dispersion_empirical >= 1 * dispersion_fit)$gene_id
expressed_genes = disp.gene
cds = setOrderingFilter(cds,disp.gene)
plot_ordering_genes(cds)

cds <- reduceDimension(cds, 
                       max_components = 2,
                       method = 'DDRTree')
cds <- orderCells(cds)




saveRDS(cds, file.path(PROJECT_ROOT, "泛癌S", "h5ad", "上皮细胞", "172", "monocel2", "重新拆分的数据", "上皮细胞Monocle2", "上皮细胞_胃癌_monocle.rds"))


HSMM_myo <- readRDS(file = file.path(PROJECT_ROOT, "泛癌S", "h5ad", "上皮细胞", "172", "monocel2", "重新拆分的数据", "上皮细胞Monocle2", "上皮细胞_胃癌_monocle.rds"))



plot_cell_trajectory(HSMM_myo, color_by = "celltype_3") +
  facet_wrap(~celltype_3, nrow = 2)

plot_cell_trajectory(HSMM_myo, color_by = "group") +
  facet_wrap(~group, nrow = 2)
plot_cell_trajectory(HSMM_myo, color_by = "cancer") +
  facet_wrap(~cancer, nrow = 2)

cds = HSMM_myo
Time_diff <- differentialGeneTest(cds, cores = 32, 
                                  fullModelFormulaStr = "~sm.ns(Pseudotime)")
Time_diff <- Time_diff[,c(5,2,3,4,1,6,7)] #把gene放前面，也可以不改
Time_diff <- subset(Time_diff, qval < 0.05)
Time_diff<- Time_diff[Time_diff$use_for_ordering == TRUE, ]
write.csv(Time_diff, "~/daima/文章整理代码/补充图S12/胃相关疾病.csv", row.names = T)

Time_genes <- Time_diff %>% pull(gene_short_name) %>% as.character()
any(is.na(Time_genes))
p=plot_pseudotime_heatmap(cds[Time_genes,],   return_heatmap=T)#num_clusters=4,show_rownames=T,
p
# ggsave("~/daima/文章整理代码/补充图S12/胃相关疾病Time_heatmapAll.pdf", p, width = 5, height = 6)

ggsave(
  filename = "~/daima/文章整理代码/补充图S12/胃相关疾病Time_heatmapAll.png",
  plot = p,
  width = 5,
  height = 5,
  dpi = 300,
  bg = "white"
)

dev.off()


p$tree_row

clusters <- cutree(p$tree_row, k = 6)
clustering <- data.frame(clusters)
clustering[,1] <- as.character(clustering[,1])
colnames(clustering) <- "Gene_Clusters"
table(clustering)

head(clustering)
genes_cluster1 <- subset(clustering, Gene_Clusters %in%c(1,6))
write.csv(genes_cluster1, "~/daima/文章整理代码/补充图S12/胃相关疾病_所欲沿着肿瘤分化的路径的升高基因.csv", row.names = T)

file_path <- "~/daima/文章整理代码/补充图S12/胃相关疾病_所欲沿着肿瘤分化的路径的升高基因.csv"
genes_cluster1 <- read.csv(file_path, row.names = 1, check.names = FALSE)
head(genes_cluster1)
table(genes_cluster1$Gene_Clusters)
genes_cluster1 = row.names(genes_cluster1)


library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)
library(ggplot2)
library(stringr)

gtp_gene_list <- genes_cluster1

gtp_gene_ids <- bitr(
  gtp_gene_list,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

gtp_gene_ids <- gtp_gene_ids[!duplicated(gtp_gene_ids$ENTREZID), ]

kegg_result <- enrichKEGG(
  gene = gtp_gene_ids$ENTREZID,
  organism = "hsa",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2
)

kegg_result <- setReadable(
  kegg_result,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID"
)

kegg_df <- as.data.frame(kegg_result)

print(kegg_df$Description)
head(kegg_df[, c("Description", "geneID")])


library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)
library(ggplot2)
library(stringr)

kegg_result <- setReadable(
  kegg_result,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID"
)

kegg_df <- as.data.frame(kegg_result)

select_terms <- c(
  "Ribosome",
  "Lipid and atherosclerosis",
  "Fluid shear stress and atherosclerosis",
  "NF-kappa B signaling pathway",
  "TNF signaling pathway",
  "IL-17 signaling pathway",
  "Epithelial cell signaling in Helicobacter pylori infection",
  "NOD-like receptor signaling pathway",
  "Viral protein interaction with cytokine and cytokine receptor",
  "Antigen processing and presentation",
  "Hematopoietic cell lineage",
  "Intestinal immune network for IgA production",
  "Th17 cell differentiation",
  "Th1 and Th2 cell differentiation",
  "Inflammatory bowel disease"
)

pathway_group <- c(
  "Fluid shear stress and atherosclerosis" = "Hypoxia",
  
  "Ribosome" = "Metabolism",
  "Lipid and atherosclerosis" = "Metabolism",
  
  "NF-kappa B signaling pathway" = "Tumor progression",
  "TNF signaling pathway" = "Tumor progression",
  "IL-17 signaling pathway" = "Tumor progression",
  "Epithelial cell signaling in Helicobacter pylori infection" = "Tumor progression",
  
  "NOD-like receptor signaling pathway" = "Tumor microenvironment",
  "Viral protein interaction with cytokine and cytokine receptor" = "Tumor microenvironment",
  
  "Antigen processing and presentation" = "Tumor immune microenvironment",
  "Hematopoietic cell lineage" = "Tumor immune microenvironment",
  "Intestinal immune network for IgA production" = "Tumor immune microenvironment",
  "Th17 cell differentiation" = "Tumor immune microenvironment",
  "Th1 and Th2 cell differentiation" = "Tumor immune microenvironment",
  
  "Inflammatory bowel disease" = "Inflammation-associated progression"
)

group_cols <- c(
  "Hypoxia" = "#B65E2E",
  "Metabolism" = "#4C9DCB",
  "Tumor progression" = "#D95F9F",
  "Tumor microenvironment" = "#4DAF4A",
  "Tumor immune microenvironment" = "#E6B85C",
  "Inflammation-associated progression" = "#9E78B8"
)

group_order <- names(group_cols)

kegg_plot_df <- kegg_df %>%
  dplyr::filter(Description %in% select_terms) %>%
  dplyr::mutate(
    pathway_class = pathway_group[Description],
    pathway_class = factor(pathway_class, levels = group_order),
    neg_log10_padj = -log10(pmax(p.adjust, .Machine$double.xmin))
  ) %>%
  dplyr::arrange(
    pathway_class,
    dplyr::desc(p.adjust)
  )

# 这里非常关键：
# 当前 data.frame 排列顺序 = 你想要的“从上到下”顺序
# ggplot 的 y 轴因子默认是从下到上显示，所以要反转 levels
kegg_plot_df$Description <- factor(
  kegg_plot_df$Description,
  levels = rev(kegg_plot_df$Description)
)

p <- ggplot(
  kegg_plot_df,
  aes(x = neg_log10_padj, y = Description, fill = pathway_class)
) +
  geom_col(
    width = 0.75,
    color = "white",
    linewidth = 0.4
  ) +
  geom_text(
    aes(label = Description),
    x = 0.15,
    hjust = 0,
    size = 4.3,
    color = "black",
    family = "Arial"
  ) +
  scale_fill_manual(
    values = group_cols,
    breaks = group_order,
    name = NULL,
    drop = FALSE
  ) +
  scale_x_continuous(
    expand = expansion(mult = c(0, 0.08))
  ) +
  labs(
    x = expression(-log[10]("adjusted p value")),
    y = NULL,
    title = "Enrichment of selected KEGG pathways"
  ) +
  theme_classic(base_size = 15, base_family = "Arial") +
  theme(
    plot.title = element_text(
      size = 18,
      face = "bold",
      hjust = 0,
      color = "black"
    ),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),
    axis.text.x = element_text(
      size = 13,
      color = "black"
    ),
    axis.title.x = element_text(
      size = 15,
      color = "black",
      margin = margin(t = 8)
    ),
    axis.line.x = element_line(
      linewidth = 0.45,
      color = "black"
    ),
    axis.ticks.x = element_line(
      linewidth = 0.45,
      color = "black"
    ),
    legend.position = "right",
    legend.text = element_text(
      size = 12,
      color = "black"
    ),
    legend.key.size = unit(0.45, "cm"),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(8, 14, 8, 8)
  )

print(p)


ggsave(
  filename = "~/daima/文章整理代码/补充图S12/胃KEGG.pdf",
  plot = p,
  width = 8,
  height = 8,
  device = cairo_pdf,
  bg = "white"
)


library(nichenetr) # Please update to v2.0.4
library(Seurat)
library(SeuratObject)
library(tidyverse)


pbmc = readRDS(file.path(PROJECT_ROOT, "泛癌S", "h5ad", "nichenet3", "胃癌_筛选.rds"))
table(pbmc$cancer)
Idents(pbmc) = pbmc$celltype_3
table(pbmc$celltype_1)
table(pbmc$celltype_3)
# table(pbmc$combined_batch)

# table(pbmc$celltype_3)
# selected_cell_types <- c("AT2","AT2 proliferating","LUAD_Malignant cells","LUAD_Malignant cells_low_cnv")
# pbmc <- subset(pbmc, subset = celltype_3 %in% selected_cell_types)
# table(pbmc$celltype_3)
# table(pbmc$cancer)


pbmc@meta.data %>% head()
pbmc@meta.data$celltype_3 %>% table() 



lr_network <- readRDS("~/bao/nichenetr/lr_network_human_21122021.rds")
ligand_target_matrix <- readRDS("~/bao/nichenetr/ligand_target_matrix_nsga2r_final.rds")
weighted_networks <- readRDS("~/bao/nichenetr/weighted_networks_nsga2r_final.rds")


lr_network <- lr_network %>% distinct(from, to)
head(lr_network)

ligand_target_matrix[1:5,1:5]
head(weighted_networks$lr_sig) 
head(weighted_networks$gr)


receiver = c("STAD_Malignant cells")
expressed_genes_receiver <- get_expressed_genes(receiver, pbmc, pct = 0.05)

all_receptors <- unique(lr_network$to)  
expressed_receptors <- intersect(all_receptors, expressed_genes_receiver)
potential_ligands <- lr_network %>% filter(to %in% expressed_receptors) %>% pull(from) %>% unique()



sender_celltypes <- unique(pbmc$celltype_3)
list_expressed_genes_sender <- sender_celltypes %>% unique() %>% lapply(get_expressed_genes, pbmc, 0.05)
expressed_genes_sender <- list_expressed_genes_sender %>% unlist() %>% unique()

potential_ligands_focused <- intersect(potential_ligands, expressed_genes_sender) 
length(expressed_genes_sender)
length(potential_ligands)
length(potential_ligands_focused)




file <- "~/daima/文章整理代码/补充图S12/胃相关疾病_所欲沿着肿瘤分化的路径的升高基因.csv"
geneset_oi <- read.csv(file, row.names = 1, check.names = FALSE)
geneset_oi = row.names(geneset_oi)

background_expressed_genes <- expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]
length(background_expressed_genes)
length(geneset_oi)

ligand_activities <- predict_ligand_activities(geneset = geneset_oi,
                                               background_expressed_genes = background_expressed_genes,
                                               ligand_target_matrix = ligand_target_matrix,
                                               potential_ligands = potential_ligands)

# ligand_activities <- ligand_activities %>% arrange(-aupr_corrected) %>% mutate(rank = rank(desc(aupr_corrected)))
# ligand_activities

ligand_activities <- ligand_activities %>%
  arrange(desc(aupr_corrected)) %>%      # 用 dplyr::desc 只在 arrange 里
  mutate(rank = row_number()) 
ligand_activities



p_hist_lig_activity <- ggplot(ligand_activities, aes(x=aupr_corrected)) + 
  geom_histogram(color="black", fill="darkorange")  + 
  geom_vline(aes(xintercept=min(ligand_activities %>% top_n(30, aupr_corrected) %>% pull(aupr_corrected))),
             color="red", linetype="dashed", size=1) + 
  labs(x="ligand activity (PCC)", y = "# ligands") +
  theme_classic()
p_hist_lig_activity

best_upstream_ligands <- ligand_activities %>% top_n(30, aupr_corrected) %>% arrange(-aupr_corrected) %>% pull(test_ligand)
best_upstream_ligands

# vis_ligand_aupr <- ligand_activities %>% filter(test_ligand %in% best_upstream_ligands) %>%
#   column_to_rownames("test_ligand") %>% select(aupr_corrected) %>% arrange(aupr_corrected) %>% as.matrix(ncol = 1)
# (make_heatmap_ggplot(vis_ligand_aupr,
#                      "Prioritized ligands", "Ligand activity", 
#                      legend_title = "AUPR", color = "darkorange") + 
#     theme(axis.text.x.top = element_blank()))  


vis_ligand_aupr <- ligand_activities %>% dplyr::filter(test_ligand %in% best_upstream_ligands) %>%
  column_to_rownames("test_ligand") %>% dplyr::select(aupr_corrected) %>% dplyr::arrange(aupr_corrected) %>% as.matrix(ncol = 1)
(make_heatmap_ggplot(vis_ligand_aupr,
                     "Prioritized ligands", "Ligand activity", 
                     legend_title = "AUPR", color = "darkorange") + 
    theme(axis.text.x.top = element_blank()))   


active_ligand_target_links_df <- best_upstream_ligands %>%
  lapply(get_weighted_ligand_target_links,
         geneset = geneset_oi,
         ligand_target_matrix = ligand_target_matrix,
         n = 100) %>%
  bind_rows() %>% drop_na()
nrow(active_ligand_target_links_df)
head(active_ligand_target_links_df)


active_ligand_target_links <- prepare_ligand_target_visualization(
  ligand_target_df = active_ligand_target_links_df,
  ligand_target_matrix = ligand_target_matrix,
  cutoff = 0.33) 
nrow(active_ligand_target_links)
head(active_ligand_target_links)
order_ligands <- intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev()
order_targets <- active_ligand_target_links_df$target %>% unique() %>% intersect(rownames(active_ligand_target_links))
vis_ligand_target <- t(active_ligand_target_links[order_targets,order_ligands])
p = make_heatmap_ggplot(vis_ligand_target, "Prioritized ligands", "Predicted target genes",
                        color = "purple", legend_title = "Regulatory potential") +
  scale_fill_gradient2(low = "whitesmoke",  high = "purple")
p

ggsave("~/daima/文章整理代码/补充图S12/胃癌.pdf", plot = p, width = 14, height = 6)

ligand_receptor_links_df <- get_weighted_ligand_receptor_links(
  best_upstream_ligands, expressed_receptors,
  lr_network, weighted_networks$lr_sig) 

vis_ligand_receptor_network <- prepare_ligand_receptor_visualization(
  ligand_receptor_links_df,
  best_upstream_ligands,
  order_hclust = "both") 
(make_heatmap_ggplot(t(vis_ligand_receptor_network), 
                     y_name = "Ligands", x_name = "Receptors",  
                     color = "mediumvioletred", legend_title = "Prior interaction potential"))




n_ligands <- 25
n_targets <- 50

active_ligand_target_links <- prepare_ligand_target_visualization(
  ligand_target_df = active_ligand_target_links_df,
  ligand_target_matrix = ligand_target_matrix,
  cutoff = 0.33
)

cat("Number of active ligand-target links:", nrow(active_ligand_target_links), "\n")

order_ligands_all <- intersect(
  best_upstream_ligands,
  colnames(active_ligand_target_links)
)

order_ligands <- head(order_ligands_all, n_ligands) %>% rev()

target_score <- apply(
  active_ligand_target_links[, order_ligands, drop = FALSE],
  1,
  max,
  na.rm = TRUE
)

target_score <- target_score[target_score > 0]

order_targets <- names(sort(target_score, decreasing = TRUE))
order_targets <- intersect(
  order_targets,
  active_ligand_target_links_df$target %>% unique()
)

order_targets <- head(order_targets, n_targets)

vis_ligand_target <- t(
  active_ligand_target_links[order_targets, order_ligands, drop = FALSE]
)

p <- make_heatmap_ggplot(
  vis_ligand_target,
  "Prioritized ligands",
  "Predicted target genes",
  color = "purple",
  legend_title = "Regulatory potential"
) +
  scale_fill_gradient2(
    low = "whitesmoke",
    high = "purple"
  ) +
  theme(
    axis.text.x = element_text(size = 8, angle = 45, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 10),
    axis.title.x = element_text(size = 13),
    axis.title.y = element_text(size = 13),
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10)
  )

p

ggsave("~/daima/文章整理代码/补充图S12/胃癌_缩减.pdf", plot = p, width = 10, height = 6)

