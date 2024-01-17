# Figure 1 ----

# Create multiple plots
all_plt <- arrangeGrob(ace_in, shannon_in, simpson_in,
                       rw_ven, abun_silva,  
                       layout_matrix = rbind(c(1, 2, 5, 5, 5),
                                             c(3, 4, 5, 5, 5)) )

figure1 <- ggpubr::as_ggplot(all_plt) + 
  draw_plot_label(label = LETTERS[1:5], 
                  x = c(0, 0.2, 0, 0.2, 0.42),
                  y = c(1, 1, .5, .5, 1))

ggsave(filename = "Figure 1.pdf", plot = figure1, scale = 2.5,
       path = "./plots/Main", width = 5, height = 2.5, units = "in", dpi = 450)

# Figure 2 ----

fig2 <- arrangeGrob(cca_ambien, phy_tree, 
                    layout_matrix = rbind(c(1, 1, 1, 1, 2, 2,  2) ))
figure2 <- ggpubr::as_ggplot(fig2) +
  draw_plot_label(label = LETTERS[1:2], 
                  x = c(0, 0.55),
                  y = c(1, 1))

ggsave(filename = "figure2.pdf", plot = figure2, width = 5,
       height = 2.5, units = "in", dpi = 450, path = "./plots/Main", 
       scale = 2.6)
