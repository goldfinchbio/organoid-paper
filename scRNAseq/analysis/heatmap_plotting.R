#if want to use other font - need to use this package - but didn't fully work
#install.packages('showtext')

#load some libraries
library(dplyr)
library(pheatmap)
library(Seurat)
library(RColorBrewer)
library(showtext)


# this was a good resource
# https://davetang.org/muse/2018/05/15/making-a-heatmap-in-r-with-the-pheatmap-package/
# for choosing colors good reference
# https://sahirbhatnagar.com/heatmap
# http://eyetracking.upol.cz/color/

# read .csvs when working locally 
avg_reduced_clusters <- read.csv('/mnt/cellranger/scRNAseq/output_eva/write/avg_reduced_clusters_oct20.csv', header = TRUE, row.names = 'gene')
avg_reduced_clusters$X <- NULL

# change the columnnames (there must be a better way of doing this)
avg_reduced_clusters <- avg_reduced_clusters %>% 
                            rename(
                            Pluripotent_D0 = Pluripotent,
                            Stroma_D7 = Stroma.x,
                            Glial_D7 = Glial.x,
                            Nephron.Progenitor_D10 = Nephron.Progenitor.x,
                            Stroma_D10 = Stroma.y,
                            Glial_D10 = Glial.y,
                            Nephron.Progenitor_D12 = Nephron.Progenitor.y,
                            Stroma_D12 = Stroma.x.x,
                            Glial_D12 = Glial.x.x,
                            Nephron.Progenitor_D14 = Nephron.Progenitor.x.x,
                            Stroma_D14 = Stroma.y.y,
                            Glial_D14 = Glial.y.y,
                            Nephron.Progenitor_D24 = Nephron.Progenitor.y.y,
                            Podocyte_D24 = Podocyte.x,
                            Proximal.Tubule_D24 = Proximal.Tubule.x,
                            Thick.Ascending.Limb_D24 = Thick.Ascending.Limb.x,
                            Distal.Tubule_D24 = Distal.Tubule.x,
                            Stroma_D24 = Stroma.x.x.x,
                            Glial_D24 = Glial.x.x.x,
                            Nephron.Progenitor_D26 = Nephron.Progenitor.x.x.x,
                            Podocyte_D26 = Podocyte.y,
                            Proximal.Tubule_D26 = Proximal.Tubule.y,
                            Thick.Ascending.Limb_D26 = Thick.Ascending.Limb.y,
                            Distal.Tubule_D26 = Distal.Tubule.y,
                            Stroma_D26 = Stroma.y.y.y,
                            Glial_D26 = Glial.y.y.y,
                            Melanocyte_D26 = Melanocyte,
                            Nephron.Progenitor_D28 = Nephron.Progenitor.y.y.y,
                            Podocyte_D28 = Podocyte,
                            Proximal.Tubule_D28 = Proximal.Tubule,
                            Thick.Ascending.Limb_D28 = Thick.Ascending.Limb,
                            Distal.Tubule_D28 = Distal.Tubule,
                            Stroma_D28 = Stroma,
                            Glial_D28 = Glial
                            )

# replace the NAs with zero
avg_reduced_clusters <- avg_reduced_clusters %>% replace(is.na(.), 0)

# function for calculating z-score
cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}

# apply function (accross rows) and make into dataframe
avg_reduced_clusters_z <- as.data.frame(t(apply(avg_reduced_clusters, 1, cal_z_score)))

#drop the NaN rows
avg_reduced_clusters_z <- avg_reduced_clusters_z[complete.cases(avg_reduced_clusters_z), ]

# make annotations for the columns
#sample_col_red <- data.frame(sample = rep(c("Day_0", "Day_07","Day_10", "Day_12", "Day_14","Day_24", "Day_26", "Day_28"), c(1,2,3,3,3,8,9,9))) old
sample_col_red <- data.frame(sample = rep(c("D0", "D7","D10", "D12", "D14","D24", "D26", "D28"), c(1,2,3,3,3,7,8,7)))
row.names(sample_col_red) <- colnames(avg_reduced_clusters_z)


# define some colors
my_colour = list(
  sample = c(D0 = "#FAFAFA", D7 = "#EEEEEE",D10 = "#C8C8C8", 
             D12 = "#ABABAB", D14 = "#8C8C8C",D24 ="#6E6E6E", 
             D26 ="#565656", D28 = '#000000')
)

# define function to save heatmaps as .pdf
save_pheatmap_pdf <- function(x, filename, width=7, height=7) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}


# set color palette
paletteLength <- 50
color=colorRampPalette(c("blue", "white", "red"))(paletteLength)

# make zero white
myBreaks <- c(seq(min(avg_reduced_clusters_z), 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(max(avg_reduced_clusters_z)/paletteLength, max(avg_reduced_clusters_z), length.out=floor(paletteLength/2)))


#replace all labels after the '_'
labels_col = colnames(avg_reduced_clusters_z)
labels_col = gsub("\\_.+$", "", labels_col)
labels_col = gsub(".", " ", labels_col, fixed=TRUE)

# plot all genes with dendrogram and gene names
figure <- pheatmap(avg_reduced_clusters_z, 
          border_color = "grey70",
          color=color,
          annotation_colors = my_colour,
          annotation_col = sample_col_red,
          cluster_cols = FALSE,
          gaps_col = c(1,3,6,9,12,19,27),
          fontsize = 4,
          labels_col = labels_col,
          breaks = myBreaks,
          angle_col = 90,
          )

save_pheatmap_pdf(figure, "/mnt/cellranger/scRNAseq/output_eva/figures/heatmap_scRNA_all_genes_genenames_with_annot.pdf", width= 10, height = 12)


# plot all genes with dendrogram and gene names -  no day annotation
figure <- pheatmap(avg_reduced_clusters_z, 
                   border_color = "grey70",
                   color=color,
                   cluster_cols = FALSE,
                   gaps_col = c(1,3,6,9,12,19,27),
                   fontsize = 4,
                   labels_col = labels_col,
                   breaks = myBreaks,
                   angle_col = 90,
                    )

save_pheatmap_pdf(figure, "/mnt/cellranger/scRNAseq/output_eva/figures/heatmap_scRNA_all_genes_genenames_without_annot.pdf", width= 10, height = 12)


# reduced - all genes - remove labels and clustertree
figure <- pheatmap(avg_reduced_clusters_z, 
           border_color = "grey70",
           color=color,
           annotation_colors = my_colour,
           annotation_col = sample_col_red,
           show_rownames = FALSE,
           treeheight_row = 0,
           cluster_cols = FALSE,
           gaps_col = c(1,3,6,9,12,19,27),
           fontsize = 12,
           labels_col = labels_col,
           breaks = myBreaks,
           angle_col = 90,
            )

save_pheatmap_pdf(figure, "/mnt/cellranger/scRNAseq/output_eva/figures/heatmap_scRNA_all_genes_nov20_with_annot.pdf", width= 8, height = 10)

# reduced - all genes - remove labels and clustertree
figure <- pheatmap(avg_reduced_clusters_z, 
                   border_color = "grey70",
                   color=color,
                   show_rownames = FALSE,
                   treeheight_row = 0,
                   cluster_cols = FALSE,
                   gaps_col = c(1,3,6,9,12,19,27),
                   fontsize = 12,
                   labels_col = labels_col,
                   breaks = myBreaks,
                   angle_col = 90
)

save_pheatmap_pdf(figure, "/mnt/cellranger/scRNAseq/output_eva/figures/heatmap_scRNA_all_genes_nov20_without_annot.pdf", width= 8, height = 10)



################## plot reduced dataset with selected genelists ##################

genes_subset = read.csv('/mnt/cellranger/scRNAseq/output_eva/data/selected_genes_oct20.csv', header = TRUE)

avg_reduced_clusters_z_column <- tibble::rownames_to_column(avg_reduced_clusters_z, "gene")
avg_reduced_clusters_z_subset <- left_join(genes_subset, avg_reduced_clusters_z_column, by = "gene")
rownames(avg_reduced_clusters_z_subset) <- avg_reduced_clusters_z_subset$gene
avg_reduced_clusters_z_subset$gene <- NULL

# plot heatmap with selected genes and italisized gene names
newnames <- lapply(
  rownames(avg_reduced_clusters_z_subset),
  function(x) bquote(italic(.(x))))


myBreaks <- c(seq(min(avg_reduced_clusters_z_subset), 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(max(avg_reduced_clusters_z_subset)/paletteLength, max(avg_reduced_clusters_z_subset), length.out=floor(paletteLength/2)))



figure <- pheatmap(avg_reduced_clusters_z_subset, 
           #border_color = "grey60",
           border_color = "grey70",
           color=color,
           annotation_colors = my_colour,
           annotation_col = sample_col_red,
           cluster_cols = FALSE,
           cluster_rows = FALSE,
           gaps_col = c(1,3,6,9,12,19,27),
           fontsize = 14,
           angle_col = 90,
           labels_col = labels_col,
           breaks = myBreaks,
           labels_row = as.expression(newnames),
          )

save_pheatmap_pdf(figure, "/mnt/cellranger/scRNAseq/output_eva/figures/heatmap_scRNA_selected_genes_nov20_with_annot.pdf", width= 10, height = 7)

figure <- pheatmap(avg_reduced_clusters_z_subset, 
                   #border_color = "grey60",
                   border_color = "grey70",
                   color=color,
                   cluster_cols = FALSE,
                   cluster_rows = FALSE,
                   gaps_col = c(1,3,6,9,12,19,27),
                   fontsize = 14,
                   angle_col = 90,
                   labels_col = labels_col,
                   breaks = myBreaks,
                   labels_row = as.expression(newnames),
)

save_pheatmap_pdf(figure, "/mnt/cellranger/scRNAseq/output_eva/figures/heatmap_scRNA_selected_genes_nov20_without_annot.pdf", width= 10, height = 7)



########################################   replot everything with Roboto font     #########################################################



font_add_google(name = "Roboto", family = "Roboto")
showtext_auto()

# reduced - all genes - remove labels and clustertree
figure <- pheatmap(avg_reduced_clusters_z, 
                   border_color = "grey70",
                   color=color,
                   annotation_colors = my_colour,
                   annotation_col = sample_col_red,
                   show_rownames = FALSE,
                   treeheight_row = 0,
                   cluster_cols = FALSE,
                   gaps_col = c(1,3,6,9,12,19,27),
                   fontsize = 12,
                   labels_col = labels_col,
                   breaks = myBreaks,
                   angle_col = 90,
                   family = 'Roboto'
)
save_pheatmap_pdf(figure, "/mnt/cellranger/scRNAseq/output_eva/figures/heatmap_scRNA_all_genes_nov20_with_annot_roboto.pdf", width= 8, height = 10)


# reduced - all genes - remove labels and clustertree - no day annotation 
figure <- pheatmap(avg_reduced_clusters_z, 
                   border_color = "grey70",
                   color=color,
                   show_rownames = FALSE,
                   treeheight_row = 0,
                   cluster_cols = FALSE,
                   gaps_col = c(1,3,6,9,12,19,27),
                   fontsize = 12,
                   labels_col = labels_col,
                   breaks = myBreaks,
                   angle_col = 90,
                   family = 'Roboto'
)
save_pheatmap_pdf(figure, "/mnt/cellranger/scRNAseq/output_eva/figures/heatmap_scRNA_all_genes_nov20_without_annot_roboto.pdf", width= 8, height = 10)



################## plot reduced dataset with selected genelists ##################


# plot heatmap with selected genes and italisized gene names
newnames <- lapply(
  rownames(avg_reduced_clusters_z_subset),
  function(x) bquote(italic(.(x))))

figure <- pheatmap(avg_reduced_clusters_z_subset, 
                   #border_color = "grey60",
                   border_color = "grey70",
                   color=color,
                   annotation_colors = my_colour,
                   annotation_col = sample_col_red,
                   cluster_cols = FALSE,
                   cluster_rows = FALSE,
                   gaps_col = c(1,3,6,9,12,19,27),
                   fontsize = 12,
                   angle_col = 90,
                   labels_col = labels_col,
                   breaks = myBreaks,
                   labels_row = as.expression(newnames),
                   family = 'Roboto'
)

save_pheatmap_pdf(figure, "/mnt/cellranger/scRNAseq/output_eva/figures/heatmap_scRNA_selected_genes_nov20_with_annot_roboto.pdf", width= 10, height = 7)


figure <- pheatmap(avg_reduced_clusters_z_subset, 
                   #border_color = "grey60",
                   border_color = "grey70",
                   color=color,
                   cluster_cols = FALSE,
                   cluster_rows = FALSE,
                   gaps_col = c(1,3,6,9,12,19,27),
                   fontsize = 12,
                   angle_col = 90,
                   labels_col = labels_col,
                   labels_row = as.expression(newnames),
                   breaks = myBreaks,
                   family = 'Roboto'
)

save_pheatmap_pdf(figure, "/mnt/cellranger/scRNAseq/output_eva/figures/heatmap_scRNA_selected_genes_nov20_without_annot_roboto.pdf", width= 10, height = 7)

