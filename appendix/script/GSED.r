library(clusterProfiler)
library(org.EcK12.eg.db)
library(AnnotationDbi)

library(clusterProfiler)
library(AnnotationDbi)
library(org.EcK12.eg.db)
library(enrichplot)
library(ggplot2)

# ==== (A) Load your DE table ====
# Example path: de_results.csv with columns: SYMBOL, log2FC, pvalue
# tbl <- read.csv("de_results.csv")

# --- Demo data (remove when using your real table) ---
tbl <- data.frame(
  SYMBOL = c("lacZ","rpoS","oxyR","recA","dnaK","groL","aceA","aceB","fumC","sodA","katG"),
  log2FC = c( 2.3, -1.4,  1.1, 0.8, 1.6, 1.9, -1.2, -1.0, 0.9, 1.3, 0.7),
  pvalue = c(1e-6, 5e-4, 2e-3, 0.04, 1e-5, 2e-6, 1e-3, 6e-3, 0.02, 3e-4, 0.05)
)

# ==== (B) Map to ENTREZ IDs ====
map_to_entrez <- function(keys, from_keytype="SYMBOL"){
  na.omit(mapIds(org.EcK12.eg.db, keys=keys, keytype=from_keytype,
                 column="ENTREZID", multiVals="first"))
}
tbl$ENTREZID <- map_to_entrez(tbl$SYMBOL)

# ==== (C) Build ranking metric ====
# Common choice: signed significance = log2FC * -log10(pvalue)
tbl$rank_score <- with(tbl, log2FC * -log10(pvalue))

# Named numeric vector: names = ENTREZID; values = rank metric
ranks <- with(na.omit(tbl[,c("ENTREZID","rank_score")]),
              setNames(rank_score, ENTREZID))
ranks <- sort(ranks, decreasing = TRUE)

length(ranks); head(ranks)
ranks
#=================================================================
#+++++++++++++++++++++++Step Run GSEA (GO & KEGG) ++++++++++++++++
#=================================================================
# ---- GO (Biological Process) ----
gsea_go <- gseGO(
  geneList      = ranks,
  OrgDb         = org.EcK12.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",
  minGSSize     = 10,
  maxGSSize     = 500,
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  eps           = 1e-10,      # numerical stability
  verbose       = FALSE
)

# ---- KEGG (organism = eco) ----
gsea_kegg <- gseKEGG(
  geneList      = ranks,
  organism      = "eco",
  minGSSize     = 10,
  maxGSSize     = 500,
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  eps           = 1e-10,
  verbose       = FALSE
)


library(org.EcK12.eg.db)
GENE <- c("lacZ","rpoS","oxyR","recA","dnaK","groL","aceA","aceB","fumC","sodA","katG")
mapped <- AnnotationDbi::select(org.EcK12.eg.db,
                                keys = GENE,
                                columns = "GO",
                                keytype = "GID")
