source("Functions.R")
load("taxonomy_midas.RData")
# Look for extra_dt in the clean_DADA.R file
ps_midas <- phyloseq(otu_table(seqtab_nochim, taxa_are_rows = F),
               sample_data(extra_dt), tax_table(taxonomy) )

plot_boxplot(ps_midas, "Localidad", c("Shannon", "Simpson"))
relative_abun(phylo = ps_midas, n_top = 25)
