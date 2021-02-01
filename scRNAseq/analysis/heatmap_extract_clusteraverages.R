install.packages("textshape")

#load some libraries
library(dplyr)
library(pheatmap)
library(Seurat)
library(RColorBrewer)
library(textshape)


# this was a good resource
# https://davetang.org/muse/2018/05/15/making-a-heatmap-in-r-with-the-pheatmap-package/
# for choosing colors good reference
# https://sahirbhatnagar.com/heatmap
# http://eyetracking.upol.cz/color/

### this is to extract only Nanostring genes ####

# read in the Seurat objects and make the average expression DF
list_timpoints = list("Day_0", "Day_07","Day_10", "Day_12", "Day_14","Day_24", "Day_26", "Day_28")
basefile = "/mnt/cellranger/scRNAseq/output/"
midfile = "_in_vitro_integration/"
endfile = "_in_vitro_integration.RDS"

# read in files with genenames for Nanostring
nanost_genes_all <- read.csv('/mnt/cellranger/scRNAseq/output_eva/data/nanostring_genes.csv', 
                             header = FALSE, col.names = 'gene')

# create dataframe names for averages
avg_reduced_clusters <- nanost_genes_all

# how to collapse cell types
collapse <- c("S-Shaped Body", "Comma-Shaped Body", "Renal Vessicle", "Cap Mesenchyme", "Podocyte Precursor")
cell_levels_collapse <- c("Podocyte", "Proximal Tubule", "Distal Tubule", "Thick Ascending Limb", "Podocyte Precursor", "Nephron Progenitor","Epithelial",
                          "Endothelial", "Stroma", "Glial", "Neural", "Melanocyte", "Cartilage", "Muscle",  "Pluripotent", "Cell Cycle")


#loop to select nanostring genes for heatmap - reduced clusters (= Nephron progenitors only)
for (i in list_timpoints){
  file_name = paste(basefile,i,midfile,i,endfile, sep = "") #merge the string
  seur <- readRDS(file_name) # read the file
  temp_meta <- tibble(cell_name = rownames(seur@meta.data),
                      cell_type_manual = seur$cell_type_manual)
  temp_meta <- temp_meta %>% mutate(cell_type_manual = ifelse(cell_type_manual == "Neural", "Neuronal", cell_type_manual),
                                    cell_type_simple = ifelse(cell_type_manual %in% collapse, "Nephron Progenitor", cell_type_manual),
                                    cell_type_simple = factor(cell_type_simple, levels = cell_levels_collapse)) %>% 
    select(cell_name, cell_type_simple) %>% column_to_rownames("cell_name")
  
  seur <- AddMetaData(seur, temp_meta)
  Idents(object = seur) <- "cell_type_simple" # if want to read in all the clusters have to delete this line
  day_select <- AverageExpression(seur)
  avg_expression <- day_select$RNA
  avg_expression <- tibble::rownames_to_column(avg_expression, "gene")
  ## add in column names for the day
  avg_reduced_clusters = left_join(avg_reduced_clusters, avg_expression, by = "gene")
}


# save the raw file as .csvs
write.csv(avg_reduced_clusters,'/mnt/cellranger/scRNAseq/output_eva/write/avg_reduced_clusters_oct20.csv')

### this is to extract all genes ####

list_timpoints = list("Day_07","Day_10", "Day_12", "Day_14","Day_24", "Day_26", "Day_28")

# create the first D_0 file "manually"

file_name = paste(basefile,"Day_0",midfile,"Day_0",endfile, sep = "") #merge the string
seur <- readRDS(file_name) # read the file
temp_meta <- tibble(cell_name = rownames(seur@meta.data),
                    cell_type_manual = seur$cell_type_manual)
temp_meta <- temp_meta %>% mutate(cell_type_manual = ifelse(cell_type_manual == "Neural", "Neuronal", cell_type_manual),
                                  cell_type_simple = ifelse(cell_type_manual %in% collapse, "Nephron Progenitor", cell_type_manual),
                                  cell_type_simple = factor(cell_type_simple, levels = cell_levels_collapse)) %>% 
  select(cell_name, cell_type_simple) %>% column_to_rownames("cell_name")

seur <- AddMetaData(seur, temp_meta)
Idents(object = seur) <- "cell_type_simple" # if want to read in all the clusters have to delete this line
day_select <- AverageExpression(seur)
avg_expression <- day_select$RNA
avg_expression <- tibble::rownames_to_column(avg_expression, "gene")



#loop to select nanostring genes for heatmap - reduced clusters (= Nephron progenitors only)
for (i in list_timpoints){
  file_name = paste(basefile,i,midfile,i,endfile, sep = "") #merge the string
  seur <- readRDS(file_name) # read the file
  temp_meta <- tibble(cell_name = rownames(seur@meta.data),
                      cell_type_manual = seur$cell_type_manual)
  temp_meta <- temp_meta %>% mutate(cell_type_manual = ifelse(cell_type_manual == "Neural", "Neuronal", cell_type_manual),
                                    cell_type_simple = ifelse(cell_type_manual %in% collapse, "Nephron Progenitor", cell_type_manual),
                                    cell_type_simple = factor(cell_type_simple, levels = cell_levels_collapse)) %>% 
    select(cell_name, cell_type_simple) %>% column_to_rownames("cell_name")
  
  seur <- AddMetaData(seur, temp_meta)
  Idents(object = seur) <- "cell_type_simple" # if want to read in all the clusters have to delete this line
  day_select <- AverageExpression(seur)
  avg_expression_t <- day_select$RNA
  avg_expression_t <- tibble::rownames_to_column(avg_expression_t, "gene")
  ## add in column names for the day
  avg_expression = full_join(avg_expression, avg_expression_t, by = "gene")
}

# save the raw file as .csvs
write.csv(avg_expression,'/mnt/cellranger/scRNAseq/output_eva/write/avg_expression_all_genes_oct20.csv')
