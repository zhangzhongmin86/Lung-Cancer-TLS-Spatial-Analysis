# Portable project path. Set TLS_PROJECT_ROOT to the directory containing the input data.
PROJECT_ROOT <- normalizePath(Sys.getenv("TLS_PROJECT_ROOT", unset = "."), winslash = "/", mustWork = FALSE)



####H5AD转rds#####

library(devtools)
library(usethis)

library(CellChat)
library(patchwork)
library(SeuratDisk)
library(Seurat) 

library(harmony)
library(Seurat)
library(SeuratDisk)
library(rhdf5)
library(readr)
library(Seurat)
library(reticulate)




pbmc = readRDS(file.path(PROJECT_ROOT, "转h5ad", "取处组1的所有样本取验证", "group1_adata.rds"))

table(pbmc$celltype_3)
pbmc <- NormalizeData(pbmc)
pbmc.matrix = pbmc@assays$RNA$data
meta <- pbmc@meta.data


cellchat=createCellChat(object = pbmc.matrix,meta = meta,group.by = "celltype_3")
cellchat <- setIdent(cellchat, ident.use = "celltype_3") 
levels(cellchat@idents) 
groupSize <- as.numeric(table(cellchat@idents)) 

CellChatDB <- CellChatDB.human
showDatabaseCategory(CellChatDB)
CellChatDB.use <- CellChatDB
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat) 

cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

cellchat <- computeCommunProb(cellchat,trim = 0.1)
cellchat <- filterCommunication(cellchat, min.cells = 1)
names(cellchat)
# 提取关注得细胞间通讯关系
df.net <- subsetCommunication(cellchat)
head(df.net)
class(df.net)
cellchat <- computeCommunProbPathway(cellchat)
cellchat@netP$pathways

# 整合通讯网络结果
cellchat <- aggregateNet(cellchat)
groupSize <- as.numeric(table(cellchat@idents))
saveRDS(cellchat,file.path(PROJECT_ROOT, "转h5ad", "取处组1的所有样本取验证", "group1_adata_cellchat.rds"))




cellchat = readRDS(file.path(PROJECT_ROOT, "superscc", "结果位置", "cellcahat", "免疫内皮细胞.rds"))
cellchat@netP$pathways


ptm = Sys.time()
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")


pathways.show <- c("CCL") 
# Hierarchy plot
# Here we define `vertex.receive` so that the left portion of the hierarchy plot shows signaling to fibroblast and the right portion shows signaling to immune cells 
vertex.receiver = seq(1,4) # a numeric vector. 
netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver)
# Circle plot
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")


netAnalysis_contribution(cellchat, signaling = pathways.show)


pairLR.CXCL <- extractEnrichedLR(cellchat, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR.CXCL[2,] # show one ligand-receptor pair
LR.show
# Hierarchy plot
vertex.receiver = seq(1,4) # a numeric vector
netVisual_individual(cellchat, signaling = pathways.show,  pairLR.use = LR.show, vertex.receiver = vertex.receiver)
#> [[1]]
# Circle plot
netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")


table(cellchat@idents)
netVisual_bubble(cellchat, sources.use = 2, targets.use = c(1:14), remove.isolate = FALSE)



netVisual_bubble(cellchat, sources.use = 4, targets.use = c(5:11), signaling = c("CCL","CXCL"), remove.isolate = FALSE)

cellchat@netP$pathways
pathways.show <- c("CXCL") 
ptm = Sys.time()
# Compute the network centrality scores
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)



netAnalysis_contribution(cellchat, signaling = pathways.show)


pairLR.CXCL <- extractEnrichedLR(cellchat, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR.CXCL[3,] # show one ligand-receptor pair
LR.show
# Hierarchy plot
vertex.receiver = seq(1,4) # a numeric vector
netVisual_individual(cellchat, signaling = pathways.show,  pairLR.use = LR.show, vertex.receiver = vertex.receiver)


# 1. 提取 CXCL 通路的所有显著通讯
df.cxcl <- subsetCommunication(cellchat, signaling = "CXCL")

# 2. 只保留 ligand = CXCL13
df.cxcl13 <- df.cxcl[df.cxcl$ligand == "CXCL13", ]

# 3. 查看是否存在
print(df.cxcl13[, c("source", "target", "ligand", "receptor", "interaction_name", "prob", "pval")])

if (nrow(df.cxcl13) == 0) {
  stop("当前 CellChat 结果中没有显著的 CXCL13 通讯。")
}

# 4. 构造 pairLR.use
pairLR.use <- unique(df.cxcl13["interaction_name"])

# 5A. 画 bubble plot
netVisual_bubble(
  cellchat,
  sources.use = levels(cellchat@idents),
  targets.use = levels(cellchat@idents),
  pairLR.use = pairLR.use,
  remove.isolate = FALSE,
  title.name = "CXCL13-related interactions"
)

# 5B. 画单个 LR 对的 circle plot
netVisual_individual(
  cellchat,
  signaling = "CXCL",
  pairLR.use = pairLR.use,
  layout = "circle"
)





par(mfrow=c(1,1))
netVisual_heatmap(cellchat, signaling = "CXCL", color.heatmap = "Reds")
#> Do heatmap based on a single object
#> 
par(mfrow=c(1,1))
netVisual_heatmap(cellchat, signaling = "CCL", color.heatmap = "Reds")
#> Do heatmap based on a single object

plotGeneExpression(cellchat, signaling = "CXCL", enriched.only = TRUE, type = "violin")


