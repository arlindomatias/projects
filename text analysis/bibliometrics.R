library(bibliometrix)
library(kableExtra)
library(dplyr)
library(ggplot2)
library(reshape2)
library(pubmedR)
library(gridExtra)

## Importando dados
raw <- convert2df(file = "G:/Meu Drive/Doutorado/Resultados/Bioinformática/Scientometrics/chemobrainextendido.txt", dbsource="wos",format="plaintext")

## Converter resultados
results <- biblioAnalysis(raw, sep = ";")
summary <- summary(results, k = 10, pause = F)

## Tabelas 
# Tabela geral
knitr::kable(summary$MainInformationDF,
             caption = "Information Summary",
             align = "llccl",
             format = "html") %>% 
  kable_classic(full_width = F, position = "center")

# Tabela produção anual
knitr::kable(summary$AnnualProduction, 
             caption = 
               "Annualy Production for scientific Documents",
             align = "cc", 
             format = "html") %>%
  kable_classic(full_width = F, position = "center")

# Top autores
knitr::kable(summary$MostProdAuthors, 
             caption = "Top 10 Authors", 
             align = "lclc", 
             format = "html") %>% 
  kable_classic(full_width = F, position = "center")

# Top papers
knitr::kable(summary$MostCitedPapers[,1:2], 
             caption = "10 Most Cited Papers",  
             align = "ll",
             format = "html") %>%  
  kable_classic(full_width = F, position = "center")

# Top journals
knitr::kable(summary$MostRelSources, 
             caption = "Top 10 Journals", 
             align =       
               "lc", 
             format = "html") %>% 
  kable_classic(full_width = F, position = "center")

# Top Keywords
knitr::kable(summary$MostRelKeywords, 
             caption = "Top 10 Keywords", 
             align = "lclc", 
             format = "html") %>%
  add_footnote(c("DE: Keywords Extracted from Articles","ID: Keywords
     Extracted from References of Articles"), 
               notation="alphabet") %>%  
  kable_classic(position = "center")

# Top países
knitr::kable(summary$MostProdCountries, 
             caption = "Top Productive Countries", 
             align = "lclc", 
             format = "html") %>%
  kable_classic(position = "center")


## Gráficos

# Resultados gerais
plot(x=results, pause=F)

# Produção autores
authorProdOverTime(raw, k=10)

# Journal
topSO=sourceGrowth(raw, top=1, cdf=FALSE)
DF=melt(topSO, id='Year')
ggplot(DF,aes(Year,value, group=variable, color=variable))+geom_line()

# Lotka law
L=lotka(results)
lotkaTable=cbind(L$AuthorProd[,1],L$AuthorProd[,2],L$AuthorProd[,3],L$fitted)
knitr::kable(lotkaTable, 
             caption = "Frequency Of Authors Based on Lotka's Law", 
             digits = 3, 
             align = "cccc", 
             format = "html",
             col.names = c("Number of article", "Number of authors", "Frequency based on data", "Frequency based on Lotka's law")) %>%
  kable_classic(full_width = F, position = "center")

# Collaboration network
NetMatrix <- biblioNetwork(raw, analysis = "collaboration", 
                           network = "authors", 
                           sep = ";")
networkPlot(NetMatrix, n = 10, type = "auto", Title = "collaboration Network",labelsize=1, halo = TRUE) 

# Thematic map
tm1 = thematicMap(raw, field = "DE",n.labels = 2, ngrams = 1)
plot(tm1$map)

# Associations
threeFieldsPlot(raw)

# Thematic evolution
years=c(2015)

nexus <- thematicEvolution(raw,field="DE",years=years,n=100,
                           minFreq=3, ngrams = 1)
plotThematicEvolution(nexus$Nodes,nexus$Edges)

# Co-Word analysis
CS = conceptualStructure(raw, field = "DE", 
                         method = "CA", 
                         minDegree = 4, 
                         clust = "auto",
                         stemming = FALSE, 
                         labelsize = 8, 
                         documents = 801, 
                         graph = T)

grid.arrange(CS[[4]], CS[[5]], ncol = 2, nrow = 1)

# Author colaboration
NetMatrix <- biblioNetwork(raw, analysis = "collaboration", 
                           network = "authors", 
                           sep = ";")

net = networkPlot(NetMatrix, n = 50, 
                  Title = "Author collaboration", 
                  type = "auto",
                  size = 10, 
                  size.cex = T, 
                  edgesize = 3, 
                  labelsize = 0.6)

# Keyword
Netmatrix2 = biblioNetwork(raw, analysis = "co-occurrences", 
                           network = "keywords",
                           sep = ";")

net = networkPlot(Netmatrix2, normalize = "association", 
                  weighted = T, 
                  n = 50, 
                  Title = "Keyword Co-occurrences",
                  type = "fruchterman", 
                  size = T, 
                  edgesize = 5, 
                  labelsize = 0.7)
