# Load needed packages
library(readxl)
library(dada2)
library(ShortRead)
library(tidyverse)
library(phyloseq)
library(ggrepel)
library(ggsci)
library(microbiome)
library(stringr)
library(vegan)
library(gplots)
library(viridis)
library(readxl)
library(cowplot)
library(ComplexHeatmap)
library(ggVennDiagram)
library(gridExtra)
library(ggpubr)
library(ape)
library(psadd)
library(scales)
library(ampvis2)
library(broom)
# Build Community matrix ----
# This function takes the output of the dada2 and creates a community
# matrix A(m n) where (m) are location sampled and (n) are the 
# chosen taxonomy lavel
# dat = dada2 taxonomy and otu table
# OTU_level = taxonomy lavel at whch level the data will be reduced
communi_matrix <- function(dat = csv_dt, OTU_level = 'Order', 
                           tlrn = 1e-5, rt = 1, Field = 'Per') {
  community_dt1 <-  dat%>% group_by_('Sample', OTU_level)
  switch (Field,
          Per = {community_dt1<- community_dt1 %>%
            summarise(Aggregated = sum(X._hits, na.rm = T))},
          Count = {community_dt1<- community_dt1 %>%
            summarise(Aggregated = sum(num_hits, na.rm = T))}
  )
  
  community_dt <- community_dt1 %>%
    group_by_('Sample', OTU_level) %>%
    summarise(Aggregated = sum(Aggregated, na.rm = T)) %>% 
    filter(Aggregated > tlrn) %>%
    spread_(OTU_level, 'Aggregated') %>%
    mutate_all(funs(ifelse(is.na(.),
                           yes = 0,.))) %>%
    as.data.frame()
  rownames(community_dt) <- as.character(community_dt$Sample)
  community_mt <- community_dt
  community_mt <- community_mt[,-1]
  ifelse(rt == 1, 
         return(community_dt),
         return(as.matrix(community_mt[, -c(length(community_dt),
                                            length(community_dt)-1)]))) 
}

# Box-plot ----
## This functions takes the phyloseq object and pltos a boxplot
# phylo = phyloseq object
# site = the column name of differences (default = "Localidad")
# indexs = biodiversity indexes
# n_col = number of column for each bozplot
plot_boxplot <- function(phylo, site, indexes, n_col =1) {
  richness_data <- plot_richness(phylo, measures = indexes)
  ggplot(richness_data$data, aes(Type, value, fill = Type)) +
    geom_boxplot(width = 0.7, alpha = 0.5 ) +
    facet_wrap("variable", scales = "free", ncol = n_col) +
    scale_fill_aaas() + guides(fill = F) +
    labs(x = "", y = "") +
    theme_test() +
    theme(axis.text.y = element_text(size = 12),
          axis.text.x = element_text(size = 10))
    
}

# Relatuve abundance plot ---
## This function takes the phyloseq object and creates a filled barchart
## of the top families
# phylo = phyloseq object
# n_top = top n families 
relative_abun <- function(phylo, n_top){
  top20 <- names(sort(taxa_sums(phylo), decreasing = TRUE))[1:n_top]
  ps.top20 <- transform_sample_counts(phylo, function(OTU) OTU/sum(OTU))
  ps.top20 <- prune_taxa(top20, ps.top20)
  
  mdf <- psmelt(ps.top20) %>% 
    mutate_if(is.factor, as.character) %>% 
    select(Names, Type, Family, Abundance) %>% 
    group_by(Names, Type, Family) %>% 
    summarise(Abundance = sum(Abundance)) %>% ungroup %>% 
    arrange(Type) %>%
    mutate(Names = factor(Names, levels = unique(.$Names) )) %>% 
    mutate(Family = case_when(
      grepl("midas", .$Family) ~ NA_character_,
      TRUE ~ Family )
    ) %>% 
    mutate(Family = ifelse(is.na(Family), "Unclussified", Family) ) |> 
    mutate(Names = factor(Names, 
                          levels = c("M3", "M2", "M1", "M11", "M13", "M15",
                                     "M10", "M12", "M14", "M8", "M9", "M7",
                                     "M6", "M5") ))
  # Bar plot
  ggplot(mdf, aes(x = Names, y = Abundance, fill = Family)) +
    geom_bar(stat = "identity", position = "fill", color = "black", width = 0.9)  +
   labs(y = "Relative abundance") + 
    scale_fill_d3(alpha = 0.8) + 
    theme_classic() +
    theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5)) 
  
}


# Venn diagram ----
## This function is used by the phylo_venn function
## This function aggregates an OTU level and filter per type and tolerance
# phylo = phyloseq object
# OTU_level = taxonomy lavel to bu melt
# Type = sample type to be filtered
# tol = relative percentage to filter low families abundance
melt_phylo <- function(phylo, OTU_level, Type, tol) {
  top <- names(sort(taxa_sums(phylo), decreasing = TRUE))
  ps.top <- transform_sample_counts(phylo, function(OTU) OTU/sum(OTU))
  ps.top <- prune_taxa(top, ps.top)
  psmelt(ps.top)  %>% 
    distinct_(OTU_level, .keep_all = T) %>%  na.omit() %>% 
    filter(Abundance > tol) %>% filter(Type %in% Type)
}

## THis functions aggregates the phyloseq in uniqes names to plot the venn diagram
# phylo1 = phyloseq annotated with the silva database
# phylo2 = phyloseq annotated with the midas database
# OTU_level = taxonomy lavel to bu melt
# Type = sample type to be filtered
# tol = relative percentage to filter low families abundance
# ... = other arguments to be passed to ggVennDiagram function
phylo_venn <- function(phylo1, phylo2, OTU_level, Tipo, tol, ...) {
  phylo1 <- phyloseq::subset_samples(phylo1, Type %in% Tipo)
  phylo2 <- phyloseq::subset_samples(phylo2, Type %in% Tipo)
  phylo_silva <- melt_phylo(phylo1, OTU_level, Type = Tipo , tol)[, OTU_level]
  phylo_midas <- melt_phylo(phylo2, OTU_level, Type = Tipo, tol)[, OTU_level]
  phylo_midas <- phylo_midas[!grepl('midas', phylo_midas)]
  venn_list <- list(`S` = tolower(phylo_silva),
                    `M` = tolower(phylo_midas) )
  ggVennDiagram(venn_list,..., set_color = "black", 
                edge_size = 0.1, ) + 
    guides(fill = F) + 
    scale_fill_gradient(low = "#FFFFFF", high = "#FF0000")
}


# Theme Publication ----
# Source: https://rpubs.com/Koundy/71792
theme_Publication <- function(base_size=14, base_family="helvetica") {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(face = "bold",
                                      size = rel(1.2), hjust = 0.5),
            text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(face = "bold",size = rel(1)),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text = element_text(), 
            axis.line = element_line(colour="black"),
            axis.ticks = element_line(),
            panel.grid.major = element_line(colour="#f0f0f0"),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.position = "bottom",
            legend.direction = "horizontal",
            legend.key.size= unit(0.2, "cm"),
            legend.margin = unit(0, "cm"),
            legend.title = element_text(face="italic"),
            plot.margin=unit(c(10,5,5,5),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold")
    ))
  
}

scale_fill_Publication <- function(...){
  library(scales)
  discrete_scale("fill","Publication",manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)
  
}

scale_colour_Publication <- function(...){
  library(scales)
  discrete_scale("colour","Publication",manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)
  
}
