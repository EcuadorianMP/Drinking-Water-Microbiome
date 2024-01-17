library(readxl)
library(reshape2)
library(ggplot2)
library (RColorBrewer)
library (wesanderson)

Matriz = read_excel("data//Matriz1.xlsx")
View(Matriz)

#Transponer matriz
Matriz1 = melt (Matriz, id.vars="Localidad")

#Transponer matriz
rmMatriz1 = t(Matriz)

#Delete rows
Matriz2 = Matriz1[-(1:90), ]
colnames(Matriz2)[2]="Substance"
colnames(Matriz2)[3]="Lecture"

#As numeric
Matriz2$Lecture = as.numeric(Matriz2$Lecture)

#BAR
Boxplot1 = ggplot(Matriz2, aes(x=Localidad, y=Lecture)) +
  facet_wrap (Matriz2$Substance, scales ="free_y") +
  geom_bar(stat = "identity") +
  theme_Publication() +
  theme(strip.text = element_text(size=7), axis.text.x = element_text(angle = 45, hjust = 1), axis.title.x=element_blank())
Boxplot1

#BOXPLOT
Boxplot1 = ggplot(Matriz2, aes(x=Localidad, y=Lecture)) +
  facet_wrap (Matriz2$Substance, scales ="free_y") +
  geom_boxplot () +
  theme_Publication() +
  theme(strip.text = element_text(size=7), axis.text.x = element_text(angle = 45, hjust = 1), axis.title.x=element_blank())
Boxplot1

#DOTPLOT
Boxplot1 = ggplot(Matriz2, aes(x=Localidad, y=Lecture)) +
  facet_wrap (Matriz2$Substance, scales ="free_y") +
  geom_dotplot (binwidth = 20) +
  theme_Publication() +
  theme(strip.text = element_text(size=7), axis.text.x = element_text(angle = 45, hjust = 1), axis.title.x=element_blank())
Boxplot1

#VIOLIN
Boxplot1 = ggplot(Matriz2, aes(x=Localidad, y=Lecture)) +
  facet_wrap (Matriz2$Substance, scales ="free_y") +
  geom_violin () +
  geom_jitter(width=0.001) +
  theme_Publication() +
  theme(strip.text = element_text(size=7), axis.text.x = element_text(angle = 45, hjust = 1), axis.title.x=element_blank())
Boxplot1

#Dotplot
Boxplot1 = ggplot(Matriz2, aes(x=Localidad, y=Lecture)) +
  facet_grid(Matriz2$Substance, scales ="free_y") +
  geom_dotplot (binwidth = 0.5)+
  theme_Publication() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.title.x=element_blank())
Boxplot1

#Jitter
Boxplot1 = ggplot(Matriz2, aes(x=Localidad, y=Lecture)) +
  facet_grid(Matriz2$Substance) +
  geom_jitter(width=0.001) +
  theme_Publication() +
  theme(strip.text = element_text(size=5), axis.text.x = element_text(angle = 45, 
        hjust = 1, size = 8), axis.title.x=element_blank())
Boxplot1 

                         
#Raster
Boxplot1 = ggplot(Matriz2, aes(x=Localidad, y=Substance)) +
  geom_raster(aes(fill=Lecture), interpolate = TRUE) +
  scale_fill_gradientn(colours = terrain.colors(100), breaks=seq(1, 25, by=1)) +
  theme_Publication() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.title.x=element_blank()
        , legend.key.size = unit(0.8, "line"), legend.key.width = unit(7, "line"))
Boxplot1

#Raster
Boxplot1 = ggplot(Matriz2, aes(x=Localidad, y=Substance)) +
  facet_grid(Matriz2$Substance) +
  geom_raster(aes(fill=Lecture), interpolate=FALSE)+
  scale_fill_gradientn(colours = terrain.colors(100), breaks=seq(1, 25, by=1)) +
  theme_Publication() +
  theme(strip.text = element_text(size=15, angle = 90), axis.text.x = element_text(angle = 90, hjust = 1, size = 8), axis.title.x=element_blank()
        , legend.key.size = unit(0.8, "line"), legend.key.width = unit(7, "line"))
Boxplot1

#Point range
Boxplot1 = ggplot(Matriz2, aes(x=Localidad, y=Lecture)) +
  facet_grid(Matriz2$Substance) +
  geom_pointrange(ymin=0, ymax=40) +
  theme_Publication() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.title.x=element_blank())
Boxplot1

#STACKED BAR GRAPH
Stacked = ggplot(Matriz2, aes(x=reorder(Substance, -Lecture), Lecture, fill=Localidad)) +
  geom_bar(stat="identity") +
  scale_fill_brewer(palette = "Dark2") +
  theme_Publication() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.title.x=element_blank())
Stacked

#Kruskal-Wallis test (non-parametric one-way ANOVA)
kruskal.test(Localidad ~ Fosfatos, data = Matriz)
#Mann-Whitney U text / Wilcoxon rank-sum test
pairwise.wilcox.test(Matriz$Calcio, Matriz$Localidad,
                     p.adjust.method = "BH")
#ANOVA
res.aov <- aov(Nitritos ~ Localidad, data = Matriz)
# Summary of the analysis
summary(res.aov)


write.csv(Matriz2, "Matriz2.csv")

ggsave ("ANIMCOME1h.pdf", plot = last_plot(), device="pdf", width = 9)

ggsave ("ANIMCOME1.pdf", plot = last_plot(), device="pdf")