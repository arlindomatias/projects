#Pacotes
library("corrr")
library("FactoMineR")
library("tidyverse")
library("ggplot2")
library("ggcorrplot")
library("factoextra")
library("dplyr")
library("corrplot")
library("Hmisc")
library("Rtsne")
library("plotly")
library("umap")
library(Cairo)
library(pheatmap)
library(ggstatsplot)

dev.new(width= 15, height=15, noRStudioGD = TRUE)

#Importação dos dados
PCA <- read.delim2("C:/Users/AsRock/Documents/PCA.txt")

#Converter dados em df
PCA <- as.data.frame(PCA)

#Checar valores ausentes
colSums(is.na(PCA))

#Remover variáveis não-numéricas e substituir ","
PCA_n <- PCA[,2:ncol(PCA)]

# Converter se necessário
PCA_n[] <- lapply(PCA_n, function(x) {
  if (is.character(x) | is.factor(x)) {
    as.numeric(as.character(x))
  } else {
    x
  }
})

#Normalizar
PCA_normalizado <- scale(PCA_n)

#PCA
data.pca <- princomp(PCA_normalizado)
summary(data.pca)

#Screeplot
fviz_eig(data.pca, addlabels = TRUE)

#Biplot
fviz_pca_var(data.pca, col.var = "cos2",
             axes = c(1,2),
             labelsize = 4,
             gradient.cols = c("blue", "red"), 
             repel = TRUE)

#Biplot completo
fviz_pca_biplot(data.pca, col.var = "cos2",
             axes = c(1,2),
             labelsize = 4,
             gradient.cols = c("blue", "red"), 
             repel = TRUE)

#Biplot Circular
cir.pca <- PCA(PCA_normalizado, scale.unit = TRUE)


#Contribuição para cada componente
cir.pca$var$coord  

#Qualidade da representação
fviz_cos2(data.pca, 
          choice = "var", 
          axes = 1:2)

#Individual
fviz_pca_ind(data.pca,
             geom = "point", 
             habillage = PCA$Group,
             #col.ind = PCA$Group,
             palette = c("#3399FF", "#CC0000","#00AFBB", "#E7B800", "#FC4E07", "purple"),
             addEllipses = TRUE,
             ellipse.level=0.95,
             #alpha.ind = 0.5,
             ellipse.type = "confidence", #confidence, convex, norm, euclid, t
             legend.title = "Groups"
)


#tSNE
digits <- PCA

digits %>%
  group_by(Group) %>%
  count()

data_sample <- digits %>%
  group_by(Group) %>%
  ungroup()

data_sample %>%
  group_by(Group) %>%
  count()

X <- data_sample %>% select(-Group)
y <- data_sample %>% select(Group)

tsne_results <- Rtsne(PCA_normalizado, dims = 2, perplexity = 20, verbose = TRUE, max_iter = 10000)

tsne_df <- data.frame(
  X = tsne_results$Y[, 1],
  Y = tsne_results$Y[, 2],
  Group = y
)

colors <- c("#E6194B", "#3CB44B", "#FFE119", "#4363D8", "#F58231", "#911EB4", "#46F0F0", "#F032E6", "#BCF60C", "#FABEBE")

ggplot(tsne_df, aes(x = X, y = Y, color = factor(Group))) +
  geom_point(size = 1.5) +
  scale_color_manual(values = colors) +
  labs(
    title = "t-SNE 2-Dimensional Digit Visualization",
    x = "t-SNE Dimension 1",
    y = "t-SNE Dimension 2"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 20)
  )

#UMAP

umap.data <- PCA_normalizado

umap.labels = PCA[, "Group"] 

umap.umap <- umap(umap.data, n_components = 10, random_state = 1000) 

layout <- umap.umap[["layout"]] 

layout <- data.frame(layout) 

final <- cbind(layout, PCA$Group) 


fig <- plot_ly(final, x = ~X1, y = ~X2, color = ~PCA$Group, colors = c('#636EFA','#EF553B','#00CC96'), type = 'scatter', mode = 'markers')%>%  
  
  layout(
    
    plot_bgcolor = "#e5ecf6",
    
    legend=list(title=list(text='Group')), 
    
    xaxis = list( 
      
      title = "UMAP1"),  
    
    yaxis = list( 
      
      title = "UMAP2")) 

fig

## MANOVA
# Criar Dataframe dos PCs
PC <- data.frame(cir.pca[["ind"]][["coord"]])
PC$Group <- PCA$Group

# Teste MANOVA global
manova_res <- manova(cbind(Dim.1, Dim.2, Dim.3, Dim.4, Dim.5) ~ Group,
                     data = PC)

# Valores P MANOVA
summary(manova_res, test = "Wilks")  # ou "Pillai", "Hotelling-Lawley"

# ANOVA para cada PC
summary.aov(manova_res)

# Fazer pós teste
aov_res <- aov(Dim.1 ~ Group, data = PC)
TukeyHSD(aov_res)


ggplot(PC, aes(x = Group, y = Dim.1, fill = Group)) +
  geom_violin() +
  geom_boxplot(width = 0.1) +
  theme_minimal()



