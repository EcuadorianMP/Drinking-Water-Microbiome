#Load Libraries and Fuctions
source("Functions.R")
#### Functions for Phyloseq ####
load("RData/Taxonomy.RData") # OTU and Taxonomy for silva database
load("RData/taxonomy_midas.RData") # OTU and Taxonomy for midas database

# Metadata creatino ----
samples_out <- rownames(seqtab_nochim) 
subject <- sapply(strsplit(samples_out, "_"), "[", 1)
Location <- c("Cuenca", "Uyumbicho", "Uyumbicho", 
               "Uyumbicho", "Uyumbicho", "Uyumbicho",
               "Uyumbicho", "Cuenca", "Cuenca",
               "Guayllabamba", "Guayllabamba", "Guayllabamba", 
               "Uyumbicho", "Uyumbicho")
Type <- factor(c("School", "DW", "RW",
                 "Reservoir", "RW", "Reservoir",
                 "RW", "Reservoir", "RW",
                 "School", "Reservoir", "DW",
                 "School", "School"))

extra_dt <- data.frame(Names = subject, Location = Location, Type = Type) %>% 
  mutate(Type = factor(Type, levels = c("RW", "DW", "Reservoir", "School"))) %>% 
  dplyr::arrange(Type)
rownames(extra_dt) <- subject 
rownames(seqtab_nochim) <-  rownames(seqtab_midas) <- subject

# Phyloseq creation ----
ps <- phyloseq(otu_table(seqtab_nochim, taxa_are_rows = F),
               sample_data(extra_dt), tax_table(Final_taxonomy) )
tree_silva <- rtree(ntaxa(ps), rooted = T, tip.label = taxa_names(ps))
ps <- merge_phyloseq(ps, tree_silva)
ps_midas <- phyloseq(otu_table(seqtab_midas, taxa_are_rows = F),
                     sample_data(extra_dt), tax_table(taxonomy_midas) )
tree_midas <- rtree(ntaxa(ps_midas), rooted = T, tip.label = taxa_names(ps_midas))
ps_midas <- merge_phyloseq(ps_midas, tree_midas)
# Write Kronas ----
# Kronas per Sample type
plot_krona(subset_taxa(ps, (Class!="Chloroplast") ),
           output = "./plots/Supplmental/Kronas/Type", "Type")
# Kronas per location
plot_krona(subset_taxa(ps, (Class!="Chloroplast") ),
           output = "./plots/Supplmental/Kronas/Location", "Location")
#Kronas per sample
plot_krona(subset_taxa(ps, (Class!="Chloroplast") ),
           output = "./plots/Supplmental/Kronas/Samples", "Names")


ps <- tax_glom(ps, taxrank = "Family")
(ace_in <- plot_boxplot(phylo = ps, site = "Location", n_col = 3,
                              indexes = c("ACE")) +
    theme(strip.text = element_text(size = 9)) + rotate_x_text())

(shannon_in <- plot_boxplot(phylo = ps, site = "Localidad", n_col = 3,
                               indexes = c("Shannon")) +
    stat_compare_means(comparisons = list(c("RW", "DW")),
                       test = "wilcox.test", 
                       method.args = list(alternative = "greater"),
                       label.y = c(2.8) ) +  ylim(c(0.8, 3)) +
    theme(strip.text = element_text(size = 9)) + rotate_x_text())

(simpson_in <- plot_boxplot(phylo = ps_midas, site = "Localidad", n_col = 3,
                            indexes = c("Simpson")) +
    stat_compare_means(comparisons = list(c("RW", "DW")),
                       method.args = list(alternative = "greater"),
                       test = "wilcox.test") + ylim(c(0.5, 1.05)) +
    theme(strip.text = element_text(size = 9)) + rotate_x_text()) 

## Barplot
abun_silva <- relative_abun(phylo = ps, n_top = 11)  +
  labs(x = "Sample")
abun_midas <- relative_abun(phylo = ps_midas, n_top = 25) 

# Ven diagramm
rw_ven <- phylo_venn(phylo1 = ps, phylo2 = ps_midas, 
                     OTU_level = "Family", Tipo = "RW",
                     tol = 0.001, label = "count") 

# ggtree ----
## Subset desired taxa
# "Ruminococcaceae"
ps_interst <- subset_taxa(ps, Family %in% c("Lachnospiraceae",
                                            "Ruminococcaceae"))
ps_interst <- tax_glom(ps_interst, "Genus", NArm = TRUE)

(phy_tree <- plot_tree(ps_interst,label.tips = "Genus", 
                       shape = "Type",
                      color = "Location", base.spacing = 0.05,
                      size = "Abundance",
                      plot.margin = 0.4, text.size = 3) + 
  ggsci::scale_color_locuszoom(alpha = 0.8)    )

#### Outside Phyloseq domine ####
### Exporting phyloseq objects
# Uncomment the following two lines
#write_phyloseq(ps, "OTU")
#write_phyloseq(ps, "TAXONOMY")

# Importing raw data and formmating
otu <- read_csv("data/otu_table.csv")
taxonomy <- read_csv("data/taxonomy_table.csv")

all_data <- cbind(taxonomy, otu)
data_gt <- pivot_longer(all_data, cols = M1:M9,
                        values_to = "num_hits", names_to = "Sample") %>% 
  group_by(Sample) %>% mutate(X._hits = num_hits/sum(num_hits) )

comm_mt <- communi_matrix(dat = data_gt, OTU_level = "Family",
                          rt = 2, Field = "Count")

## Import Metadata
metadata <- read_csv("data/metadata.csv") %>% 
  separate(Muestra, sep = " ", into = c("Sample", "Trash")) %>% 
  select(-Trash, -Localidad, -`ID Inicial`, -Type,
         -Characteristics, -Estado, -Contenido, - NO2) %>% 
  # Removing LOD values
  mutate(Cl = ifelse(Cl <= 5, 0, Cl),
         PO4 = ifelse(PO4 <= 0.1, 0, PO4), 
         Fe = ifelse(Fe <= 0.07, 0, Fe) )
            
metadata_mt <- as.matrix(metadata[, 2:11])
ambien_dt <- metadata %>% select(-Sample)
rownames(metadata_mt) <- metadata$Sample


cca_cedia <- cca(comm_mt~., data = ambien_dt)
anova.cca(cca_cedia, by="term", permutations=1E4) 
# p-values
# Mg: 0.001
# Na: 0.012
# pH: 0.019
# Tm: 0.007
# K: 0.038
plot(cca_cedia, type = "p")
ambien_cca <- cca_cedia$CCA$biplot %>% 
  as.data.frame() %>%
  mutate(Variable = rownames(.),
         CCA1 = ifelse(CCA1 > 0.5, CCA1*0.5, CCA1),
         CCA2 = ifelse(Variable%in% "Fosfatos", 0.20, CCA2)) %>% 
  mutate(Var_label = case_when( 
    Variable == "PO4" ~ as.character(expression( "PO[4]^-3" )),
    Variable == "NO2" ~ as.character(expression("NO[2]^-{}")),
    Variable == "NO3" ~ as.character(expression("NO[3]^-{}")),
    Variable == "pH" ~ as.character(expression("pH^{p < 0.05}")),
    Variable == "Mg" ~ as.character(expression("Mg^{p < 0.05}")),
    Variable == "Na" ~ as.character(expression("Na^{p < 0.05}")),
    Variable == "K" ~ as.character(expression("K^{p < 0.05}")),
    Variable == "Tm" ~ as.character(expression("Tm^{p < 0.05}")),
    T ~ Variable) ) %>% 
  mutate(CCA2 = case_when(
    Variable == "Fe" ~ CCA2+0.03,
    Variable == "Cl" ~ CCA2-0.08,
    Variable == "Na" ~ CCA2-0.045,
    T ~ CCA2  ))

spec_cca <- cca_cedia$CCA$wa %>% 
  as.data.frame() %>%
  mutate(Names = rownames(.),
         CCA1 = ifelse(CCA1 > 0, CCA1*0.1, CCA1), 
         CCA2 = ifelse(CCA2 < -2, CCA2*0.1, CCA2))%>% 
  left_join(extra_dt)

( cca_ambien <- ggplot(spec_cca, aes(CCA1, CCA2)) +
  geom_point(aes( shape= Type,color  =Location), size = 3.5) +
  theme_test() +
  geom_vline(xintercept = 0, color = "gray60") +
  geom_hline(yintercept = 0, color = "gray60") +
  geom_segment(data = ambien_cca , aes(x =0, y =0, 
                                    xend= CCA1, 
                                    yend = CCA2),
               arrow = arrow(length = unit(0.05, 'cm')), alpha = 0.08,
               color = 'blue', lineend = 'square') +
  geom_text(data = ambien_cca, aes(CCA1, CCA2, label = Var_label),
            alpha = 0.7, size = 4, parse = T, color = "blue") +
  ggsci::scale_color_locuszoom(alpha = 0.8) +
  #scale_shape_manual(values=c(0, 12, 9, 10))+
  labs(x = "CCA1 (24.86%)", y = "CCA2 (21.77%)",
       shape = "Source", color = "Location")  +
  xlim(c(-0.6, 0.5)) )


# Radar plots
radar_type <- sample_data(ps) %>% data.frame %>%
  select(Sample = Names, Type, Location) %>% 
  right_join(metadata) %>% select(-Sample, -Location) %>% 
  mutate_each(funs(rescale), -Type) %>% group_by(Type) %>% 
  summarise_all(~mean(.)) %>%
  ggradar::ggradar(group.line.width = 0.5, grid.max = 0.8,
                   group.point.size = 1, background.circle.colour = "white",
                   gridline.min.colour = "white", label.gridline.max = F,
                   label.gridline.mid = F, label.gridline.min = F,
                   gridline.mid.colour = "white", legend.position = "right") +
  scale_color_aaas() +
  ggsci::scale_color_aaas(alpha = 0.6) + ggtitle("Radar plot per sample type")

radar_location <- sample_data(ps) %>% data.frame %>%
  select(Sample = Names, Type, Location) %>% 
  right_join(metadata) %>% select(-Sample, -Type) %>% 
  mutate_each(funs(rescale), -Location) %>% group_by(Location) %>% 
  summarise_all(~mean(.)) %>%
  ggradar::ggradar(group.line.width = 0.5, 
                   group.point.size = 1, background.circle.colour = "white",
                   gridline.min.colour = "white", label.gridline.max = F,
                   label.gridline.mid = F, label.gridline.min = F,
                   gridline.mid.colour = "white", legend.position = "right") +
  scale_color_aaas() +
  ggsci::scale_color_lancet(alpha = 0.6) +
  ggtitle("Radar plot per sampling location")

supplemental_radar <- arrangeGrob(radar_type, radar_location)
suplemental_fig2 <- ggpubr::as_ggplot(supplemental_radar)
ggsave(filename = "radar_plot.jpeg", plot = suplemental_fig2,
       path = "./plots/Supplmental", dpi = 300, width = 150, height = 200,
       units = "mm")

# Heatap with clusters ----
com_mt_hm <- communi_matrix(dat = data_gt, Field = 'Count', rt = 2,
                            OTU_level = 'Family')

reduced_hm <- log(com_mt_hm + 1)[, colSums(log(com_mt_hm + 1)) > 20] %>%  t 
col <- colorRampPalette(c('green3','black','red', 'red'))(35)
col <- colorRampPalette(c("gray90", "gray50", "gray10"))(35)

# Heatmap 2 description of family importance
info <- read_xlsx("data/familyht.xlsx", sheet = 1)
info_col <- c("IR" = "gray80", "AR" =  "red", "DWSS" = "orange",
              "WW" = "blue", "FC" = "yellow", "WTP" = "violet")
info2 <- info[, 2]
info2 <- as.matrix(info2)
rownames(info2) <- info$A

ha = HeatmapAnnotation(Source = extra_dt$Type,
                       Location = extra_dt$Location,
                       `Mg (mg/L)` =anno_barplot( metadata_mt[, 5]),
                       col = list(Source = c("RW" = "#1F77B4FF",
                                           "DW" =  "#FF7F0EFF",
                                           "Reservoir" = "#2CA02CFF",
                                           "School" = "#9467BDFF"),
                                  Location = c("Cuenca" = "turquoise",
                                                "Guayllabamba" = "tan",
                                                "Uyumbicho" = "turquoise4")))

# Heatmap 3, families found in midas database
tax_midas <- taxonomy_midas[, 5]
tax_midas <- tax_midas[!grepl("midas", tax_midas)]
names(tax_midas) <- NULL
in_midas <- rownames(reduced_hm) %in% tax_midas %>% as.character() %>% 
  str_replace("TRUE", "Yes") %>% str_replace("FALSE", "No")
found_midas <- matrix(in_midas, ncol = 1)
rownames(found_midas) <- rownames(reduced_hm)
colnames(found_midas) <- "In MIDAS"
col_midas <- c("Yes" = ggplot2::alpha("tomato", 0.8), "No" = "gray60")


# Drowing heatmaps
h1 <- Heatmap(reduced_hm, name = "Abundance\nlog(x +1)", col = viridis(n = 20),
              column_km = 3, row_km = 3, 
               use_raster = T,
                heatmap_legend_param = list(ncol = 1), top_annotation = ha)

h2 <- Heatmap(info2, na_col = "white", col = info_col, name = "Annotation",
              use_raster = T)
h3 <- Heatmap(found_midas, col = col_midas, name = "In MIDAS", 
              use_raster = T)
ht_list <- h1 + h2 + h3
pdf(file = "./plots/Main/Figure 3.pdf", width = 8.5, height = 7)
draw(ht_list)
dev.off()

# OTU differential analysis ----
# Subset RW and DW
ps_treatment <- subset_samples(ps, Type %in% c("RW", "DW"))
ps_treatment <- tax_glom(ps_treatment, taxrank = "Genus")
# Keep otus that were detected at least five times and at least in the
# half of the samples
criteria_5counts = genefilter_sample(ps_treatment,
                                        filterfun_sample(function(x) x > 5),
                                        A=0.5*nsamples(ps_treatment))
# Filter samples
treatment_clean = prune_taxa(criteria_5counts, ps_treatment)
treatment_clean # 17 OTUS preserved

# Transform to relative abundance
treatment_clean = transform_sample_counts(treatment_clean,
                                             function(x) 100 * x/sum(x))

otu_table(treatment_clean) <- t(otu_table(treatment_clean))

significant_OTU <- abundant_otu$data %>% group_by(Display, Group) %>% nest %>% 
  spread(key = Group, value = data) %>% 
  mutate(
    t_test = map2(DW, RW, ~{wilcox.test(.x$Abundance, .y$Abundance,
                                        alternative = "less") %>% tidy()}),
    DW = map(DW, nrow),
    RW = map(RW, nrow)
  ) %>% 
  unnest()  

rizhobium <- subset_taxa(treatment_clean, Genus %in% c("Rhizobium"))
rizhobium_box <- amp_rabund(rizhobium, group = "Type",
                            plot.flip = T,tax.add = "Phylum") +
  stat_compare_means(method = "wilcox.test", paired = T, 
                     method.args = list(alternative = "great")) +
  theme_test() + labs(color = "Type") + rotate_x_text()
ggsave(filename = "rizhobium.jpeg", path = "./plots/Supplmental", 
       plot = rizhobium_box, width = 3, height = 4, units = "cm",dpi = 300, 
       scale = 3.5)

abundant_otu <- amp_rabund(treatment_clean, group = "Type", tax.show = 5,
                           tax.add = "Phylum", plot.flip = T) +
  theme_test() + rotate_x_text(angle = 45) + 
  stat_compare_means(method = "wilcox.test", 
                     method.args = list(alternative = "less"))
ggsave(filename = "5_abundant_OTU.jpeg", path = "./plots/Supplmental", 
       plot = abundant_otu, width = 6, height = 3, units = "cm",dpi = 300, 
       scale = 4)
