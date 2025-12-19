library(DESeq2)
library(ggplot2)
library(dplyr)
library(ggrepel)
# Parse command line arguments
library(optparse)

option_list = list(
    make_option(c("-c", "--count_table"), type="character", default=NULL, 
                            help="path to count table file", metavar="character"),
    make_option(c("-m", "--metadata"), type="character", default=NULL,
                            help="path to metadata file", metavar="character"),
    make_option(c("-r", "--ref_name"), type="character", default=NULL,
                            help="reference name", metavar="character"),
    make_option(c("-g", "--gff"), type="character", default=NULL,
                            help="path to GFF file", metavar="character"),
    make_option(c("-o", "--output"), type="character", default=NULL,
                            help="output directory path", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

# Check if all required arguments are provided
if (is.null(opt$count_table) || is.null(opt$metadata) || is.null(opt$ref_name) || 
        is.null(opt$gff) || is.null(opt$output)) {
    print_help(opt_parser)
    stop("All arguments must be supplied.", call.=FALSE)
}

count_table_file <- opt$count_table
metadata_file <- opt$metadata
ref_name <- opt$ref_name
gff_file <- opt$gff
output_path <- opt$output

# count_table_file <- "/research/groups/ma1grp/home/zyu/work_2025/RNA_direct_10_Oct/bai_project/5_DESeq2/count_table/overlap_CDS/sample_1_gene.count"
# metadata_file <- "/research/groups/ma1grp/home/zyu/work_2025/RNA_direct_10_Oct/bai_project/appendix/meta_data.csv"
# ref_name <- "sample_1"
# gff_file <- "/research/groups/ma1grp/home/zyu/work_2025/RNA_direct_10_Oct/bai_project/appendix/ref/sample_1.gff3"
# output_path <- "/research/groups/ma1grp/home/zyu/work_2025/RNA_direct_10_Oct/bai_project/5_DESeq2"

print(("================ Starting DESeq2 Analysis ================"))
print(paste("Count table file:", count_table_file))
print(paste("Metadata file:", metadata_file))
print(paste("Reference name:", ref_name))
print(paste("GFF file:", gff_file))
print(paste("Output path:", output_path))
print(("========================================================="))

#=================================================================
#+++++++++++++++++++++++Step 0-1 read gff files ++++++++++++++++++
#=================================================================
gff <- rtracklayer::import(gff_file)
gff <- as.data.frame(gff)
gff_frame <- subset(gff, select = c("ID","Name","gene","start","end","strand","seqnames"))
gff_frame <- mutate(gff_frame, gene = if_else(is.na(gene), Name, gene))
# rownames(gff_frame) <- gff_frame$ID
# gff_frame$ID <- NULL
# head(gff_frame)


if(!dir.exists(paste0(output_path,"/",ref_name))){
  dir.create(paste0(output_path,"/",ref_name), recursive = TRUE)
}
setwd(paste0(output_path,"/",ref_name))
write.csv(gff_frame, "0_gff_annotation.csv", row.names = FALSE)
#=================================================================
#+++++++++++++++++++++++Step 0-2 user-defined functions ++++++++++
#=================================================================
# Function to create volcano plot
create_volcano_plot <- function(res_gff, group1, group2, cutoff_fc = 1.8) {
    df <- res_gff
    volcano_name <- paste0(group1, "_vs_", group2)
    
    # Add significance column
    df$significance <- "No"
    df$significance[df$padj < 0.05 & df$log2FoldChange > cutoff_fc] <- "Up"
    df$significance[df$padj < 0.05 & df$log2FoldChange < -cutoff_fc] <- "Down"
    df$label <- ifelse(df$padj < 0.01 & abs(df$log2FoldChange) > cutoff_fc, df$gene, NA)
    diff_gene_count <- nrow(subset(df, df$padj < 0.05 & abs(df$log2FoldChange) > cutoff_fc))
    print(paste0("Number of differentially expressed genes between ", group1, " and ", group2, ": ", diff_gene_count))
    highlight_points <- subset(df, df$padj < 0.01 & abs(df$log2FoldChange) > cutoff_fc)
    max_padj <- max(-log10(df$padj[is.finite(-log10(df$padj))]))
    
    p <- ggplot(df, aes(x = log2FoldChange, y = -log10(padj), color = significance)) +
        geom_point(alpha = 0.8, size = 2.5) +
        geom_point(data = highlight_points, aes(x = log2FoldChange, y = -log10(padj)),
                             size = 4, alpha = 0.5, shape = 21, stroke = 2) +
        geom_vline(xintercept = c(-cutoff_fc, cutoff_fc), linetype = "dashed", color = "black") +
        geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray50") +
        geom_text_repel(aes(label = label),
                                        size = 4,
                                        max.overlaps = 10,
                                        box.padding = 0.3,
                                        point.padding = 0.2) +
        scale_color_manual(values = c("Up" = "#d20815", "Down" = "#0072B2", "No" = "gray80")) +
        theme_bw() +
        theme(
            plot.title = element_text(size = 16, face = "bold"),
            axis.title = element_text(size = 15, face = "bold"),
            axis.text = element_text(size = 16),
            legend.position = "right",
            legend.title = element_text(size = 16),
            legend.text = element_text(size = 14)
        ) +
        labs(
            title = paste0(volcano_name, ": labeling Padj < 0.05 and |log2FC| > ", cutoff_fc),
            subtitle = paste0("Number of differentially expressed genes: ", diff_gene_count),
            x = "log2(Fold Change)",
            y = "-log10(Adjusted P-value)",
            color = "Regulation"
        ) +
        xlim(c(-10, 10)) +
        ylim(c(0, max_padj)) +
        scale_x_continuous(breaks = seq(-9, 9, by = 1))
    
    ggsave(paste0(volcano_name, ".pdf"), plot = p, width = 10, height = 6)
    return(p)
}

#=================================================================
#+++++++++++++++++++++++Step 1 Count matrix input ++++++++++++++++
#=================================================================
count_table <- read.table(count_table_file, row.names = 1, header = TRUE, comment.char = "#")
colnames(count_table) <- gsub("_.*$", "", colnames(count_table)) 

metadata <- read.csv(metadata_file, stringsAsFactors = FALSE)
rownames(metadata) <- metadata$SampleName  # Ensure row names match the count table


#=================================================================
#+++++++++++++++++++++++Step 2 Check the row and col name ++++++++
#=================================================================
id <- rownames(metadata)
cts <- count_table[, id]

coldata <- metadata[id, ]

check_row_col <- all(rownames(coldata) == colnames(cts))
print(paste("Row names in coldata match column names in cts:",as.character(check_row_col)))

#=================================================================
#+++++++++++++++++++++++Step 3 Read count table ++++++++++++++++++
#=================================================================
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ Group)

#=================================================================
#+++++++++++++++++++++++Step 3 Filter low count genes ++++++++++++
#=================================================================        
smallestGroupSize <- 3
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]
dds <- DESeq(dds)


#=================================================================
#+++++++++++++++++++++++Step 4 PCA plot +++++++++++++++++++++++++++
#=================================================================
vsd <- vst(dds, blind=FALSE)
rld <- rlog(dds, blind=FALSE)

plotPCA(vsd, intgroup=c("Group"))

pcaData <- plotPCA(vsd, intgroup=c("Group"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))




p <- ggplot(pcaData, aes(PC1, PC2, color=Group)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + theme_bw() +theme_bw() +
  theme(
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 14),
    legend.position = "right",
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    axis.line  = element_line(color = "black", linewidth = 0.6),
    axis.ticks = element_line(color = "black", linewidth = 0.6),
    panel.border = element_rect(color = "black",linewidth = 1.2, fill = NA)
  ) +
  ggtitle(paste("The PCA plot for", ref_name," as reference")) +
  coord_fixed()
ggsave(filename = paste0("1_PCA_", ref_name, ".pdf"), plot = p, width = 6, height = 5)


#=================================================================
#+++++++++++++++++++++++Step 4 Differential expression +++++++++++
#=================================================================

if (dir.exists("2_volcano_results")) {
  } else {
    dir.create("2_volcano_results")
  }
setwd(paste0(output_path,"/",ref_name,"/2_volcano_results") )
group1 <- unique(coldata$Group)[1]
group2 <- unique(coldata$Group)[2]
group3 <- unique(coldata$Group)[3]
rld_assay <- assay(rld)
rld_assay <- as.data.frame(rld_assay)
rld_assay$ID <- rownames(rld_assay)
norm_counts <- counts(dds, normalized = TRUE)
write.csv(as.data.frame(norm_counts), file = "Normalized_counts.csv")

all_res_gff <- data.frame()
for (i in 1:(length(unique(coldata$Group))-1)) {
  for (j in (i+1):length(unique(coldata$Group))) {
    group1 <- unique(coldata$Group)[i]
    group2 <- unique(coldata$Group)[j]
    res <- results(dds,contrast = c("Group", group1, group2))
    sample_name <- paste0(group1, "_vs_", group2)
    res_df <- as.data.frame(res)
    res_df$ID <- rownames(res_df)
    res_gff <- left_join(res_df, gff_frame, by = "ID")
    all_res_gff <- rbind(all_res_gff, res_gff)

    write.csv(res_gff, file = paste0("DESeq2_results_", sample_name, ".csv"), row.names = FALSE)
    
    diff_gene_tab <- subset(res_gff, padj < 0.05 & abs(log2FoldChange) > 1.8)
    diff_gene_tab$regulation <- ifelse(diff_gene_tab$log2FoldChange > 0, "Up", "Down")
    diff_gene_tab <- left_join(diff_gene_tab, rld_assay, by = "ID")
    # diff_gene_tab$ID <- NULL
    write.csv(diff_gene_tab, file = paste0("DESeq2_diff_genes_", sample_name, ".csv"), row.names = FALSE)
    # Save DESeq2 result
    # Create volcano plot
    create_volcano_plot(res_gff, group1, group2, cutoff_fc = 1.8)
  }
}



# #=================================================================
# #+++++++++++++++++++++++Step 5 Heatmap ++++++++++++++++++++++++++
# #=================================================================
# library(circlize)
# all_diff_genes <- all_res_gff %>% filter(padj < 0.05 & abs(log2FoldChange) > 1.8) 
# unique_diff_genes <- unique(all_diff_genes$ID)
# length(unique_diff_genes)
# # vsd_assay <- assay(vsd)[rownames(vsd) %in% unique_diff_genes, ]
# rld_assay <- assay(rld)[rownames(rld) %in% unique_diff_genes, ]
# rld_assay <- as.data.frame(rld_assay)
# rld_assay$ID <- rownames(rld_assay)
# rld_assay_gff <- left_join(rld_assay, gff_frame, by = "ID")
# rownames(rld_assay_gff) <- paste0(rld_assay_gff$ID)


# for (i in 1:(length(unique(coldata$Group))-1)) {
#   for (j in (i+1):length(unique(coldata$Group))) {
#     group1 <- unique(coldata$Group)[i]
#     group2 <- unique(coldata$Group)[j]
#     res <- results(dds,contrast = c("Group", group1, group2))
#     sample_name <- paste0(group1, "_vs_", group2)
#     res_df <- as.data.frame(res)
#     res_df$ID <- rownames(res_df)
#     ref_select <- subset(res_df, select = c("ID", "log2FoldChange", "padj"))
#     colnames(ref_select)[2:3] <- paste0(colnames(ref_select)[2:3], "-", sample_name)
#     rld_assay_gff <- left_join(rld_assay_gff, ref_select, by = "ID")


#   }
# }
# head(rld_assay_gff)
# rld_assay_gff$ID <- NULL

# setwd(paste0(output_path,"/",ref_name) )
# write.csv(rld_assay_gff, "3_heatmap_expression_values_with_annotation.csv", row.names = TRUE)
# # vsd_scaled <- t(scale(t(vsd_assay)))
# # annote_col_table <- subset(coldata, select = c("Group"))
# # annote_col_color <- list(Group = c("JE1" = "#ef7a70", "JE2" = "#6f9afd", "JE65" = "#f2d16b"))

# # annote_col <- HeatmapAnnotation(
# #   df = annote_col_table,
# #   col = annote_col_color,  # list of color mappings for annotation factors
# #   show_annotation_name = TRUE,
# #   show_legend = FALSE)

# # #++++++++++++++++++++++++Cell Setting ++++++++++++++++++++++++++++++++
# # col_fun = colorRamp2(c(8, 12, 20), c("#5474b1", "white", "#c1402d"))
# # lgd = Legend(col_fun = col_fun, title = "foo")
# # cell_fun_set = function(j, i, x, y, width, height, fill) {
# #     grid.rect(x, y, width, height, 
# #               gp = gpar(col = "black", fill = NA, lwd = 0.5))  # Black border
# # }
# # write.csv(vsd_assay, "heatmap_expression_values.csv", row.names = TRUE)
# # ht <- Heatmap(
# #   vsd_scaled,
# #   col = col_fun, # Color mapping for the heatmap
# #   name = "Expression", # Name of the heatmap
# #   cluster_rows = TRUE, # Disable row clustering
# #   cluster_columns = FALSE, # Disable column clustering
# #   show_column_names = FALSE, # Hide column names
# #   show_row_names = TRUE, # Show row names
# #   column_split = annote_col_table$Group, # Split columns by the 'Group' annotation
# #   row_title = NULL, # No title for the rows
# #   top_annotation = annote_col, # Top annotation for columns
# #   column_title_gp = gpar(fontsize = 16, fontface = "bold"), # Title appearance for columns
# #   heatmap_width = unit(15, "cm"),  # match pheatmap's cellwidth
# #   heatmap_height = unit(100, "cm"), # match pheatmap's cellheight
# #   row_names_gp = gpar(fontface = "italic", fontsize = 5)  # Adjust row names appearance

# # )
# # pdf(file = paste0("heatmap_v1.pdf"),width = 10,height = 80.5)
# # draw(ht)
# # dev.off()
