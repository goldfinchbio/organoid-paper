
# Functions:
## boxplot of seurat
seur_box_plot <- function(object, features, group){
  melt <- reshape2::melt(object@meta.data %>% select(features, group), id.vars=group)
  
  p1 <- ggplot(data = melt, aes(x=melt[,group], y=value, fill= melt[,group]))+
    geom_boxplot() +
    facet_wrap(~variable, scales = "free") +
    theme(legend.position = "none") +
    xlab(group) +
    ylab("")
  p1 + theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
  
}


## quick function to save rds files and csvs
saveDat <- function(dat, filename=NULL, csv=FALSE){
  if(save){
    if(is.null(filename)){filename = deparse(substitute(dat)) }
    readr::write_rds(dat, path = paste0(saveDir, "/", filename, ".RDS"))
    if(csv){
      write.csv(dat, file = file.path(paste0(saveDir, "/", filename, ".csv")))}}
}



## Customized DataTables
downloadableDT <- function(atable, rownames=NULL, ...) {
  require(DT)
  datatable(data = atable,
            rownames = rownames, ...,
            extensions = 'Buttons', 
            options = list( 
              dom = "Blfrtip", 
              buttons = list("copy", list(
                extend = "collection", 
                buttons = c("csv", "excel", "pdf"), 
                text = "Download"
              ) ), # end of buttons customization
              
              # customize the length menu
              lengthMenu = list( c(10, 20, -1) # declare values
                                 , c(10, 20, "All") # declare titles
              ), # end of lengthMenu customization
              pageLength = 10
              
              
            ) # end of options
            
  )
}


## replicate ggplot color palette
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}


# custom themes for jsd plots
theme_jsd1 <- function(base_size = 16, base_family = "DejaVu Sans",
                       base_line_size = base_size/22, base_rect_size = base_size/22){
  theme_grey(base_size = base_size, base_family = base_family) %+replace% 
    theme(axis.text = element_text(size = rel(0.8)), 
          axis.ticks = element_line(colour = "black"), 
          legend.key = element_rect(colour = "grey80"), 
          panel.background = element_rect(fill = "white", colour = NA), 
          panel.border = element_rect(fill = NA, colour = "grey50"), 
          panel.grid.major = element_line(colour = "grey90", size = 0.2), 
          panel.grid.minor = element_line(colour = "grey98", size = 0.5), 
          strip.background = element_rect(fill = "grey80", colour = "grey50", size = 0.2),
          axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank()
    )
}

theme_jsd2 <- function(base_size = 16, base_family = "DejaVu Sans",
                       base_line_size = base_size/22, base_rect_size = base_size/22){
  theme_grey(base_size = base_size, base_family = base_family) %+replace% 
    theme(axis.text = element_text(size = rel(0.8)), 
          axis.ticks = element_line(colour = "black"), 
          legend.key = element_rect(colour = "grey80"), 
          panel.background = element_rect(fill = "white", colour = NA), 
          panel.border = element_rect(fill = NA, colour = "grey50"), 
          panel.grid.major = element_line(colour = "grey90", size = 0.2), 
          panel.grid.minor = element_line(colour = "grey98", size = 0.5), 
          strip.background = element_rect(fill = "grey80", colour = "grey50", size = 0.2),
          axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          legend.key.size = unit(4,"mm"),
          axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank()
    )
}

## modified theme_bw
theme_default <-  function (base_size = 16, base_family = "", base_line_size = base_size/22, 
          base_rect_size = base_size/22) 
{
  theme_grey(base_size = base_size, base_family = base_family, 
             base_line_size = base_line_size, base_rect_size = base_rect_size) %+replace% 
    theme(panel.background = element_rect(fill = "white",colour = NA), 
          panel.border = element_rect(fill = NA, colour = "grey20"), 
          panel.grid = element_line(colour = "grey92"), 
          panel.grid.minor = element_line(size = rel(0.5)), 
          strip.background = element_rect(fill = "grey85",colour = "grey20"), 
          legend.key = element_rect(fill = "white", colour = NA), complete = TRUE)
}




# centralized cell levels for plot order
cell_levels <- c("Podocyte", "Proximal Tubule", "Distal Tubule", "Thick Ascending Limb", "Collecting Duct", "Epithelial", "Podocyte Precursor", "Early Proximal Tubule", "Nephron Progenitor",
                 "S-Shaped Body", "Comma-Shaped Body", "Renal Vessicle", "Cap Mesenchyme",
                 "Endothelial", "Stroma", "Glial", "Cartilage", "Melanocyte", "Muscle", "Neuronal",
                 "Pluripotent", "Cell Cycle", "High Ribosomal", "Apoptotic")

cell_levels_collapse <- c("Podocyte", "Proximal Tubule", "Distal Tubule", "Thick Ascending Limb", "Nephron Progenitor","Epithelial",
                          "Endothelial",
                          "Stroma", "Glial", "Neuronal", "Melanocyte", "Cartilage", "Muscle", 
                          "Pluripotent")

colors <- c(gg_color_hue(8), "grey90", "grey80", "grey70", "grey60", "grey50", "grey40", "grey30", "grey20")

names(colors) <- cell_levels_collapse
