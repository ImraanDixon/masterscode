## Load required packages
library(org.Mm.eg.db)
library(org.Rn.eg.db)
library(ggplot2)
library(tidyverse)
library(RColorBrewer)
library(pheatmap)
library(dendextend)

##----------------------------------------------------------------------------------

genesofinterest <- function(x = unique(c(keys(org.Mm.eg.db, "SYMBOL"), 
                                         keys(org.Rn.eg.db, "SYMBOL")))) {
  ## Load in DE data for each dataset at each time point within the latency period
  GSE47752 <- list(
    "24h" = list(
      "Nadler" = read.csv("GSE47752/Nadler 1 Day.csv", row.names = 1),
      "Whadman" = read.csv("GSE47752/Whadman 1 Day.csv", row.names = 1)
    ),
    "72h" = list(
      "Nadler" = read.csv("GSE47752/Nadler 3 Days.csv", row.names = 1),
      "Whadman" = read.csv("GSE47752/Whadman 3 Days.csv", row.names = 1)
    ),
    "240h" = list(
      "Nadler" = read.csv("GSE47752/Nadler 10 Days.csv", row.names = 1),
      "Whadman" = read.csv("GSE47752/Whadman 10 Days.csv", row.names = 1)
    )
  )
  GSE1834 <- list(
    "1h" = read.csv("GSE1834/1h.csv", row.names = 1),
    "6h" = read.csv("GSE1834/6h.csv", row.names = 1),
    "24h" = read.csv("GSE1834/24h.csv", row.names = 1),
    "72h" = read.csv("GSE1834/72h.csv", row.names = 1),
    "240h" = read.csv("GSE1834/144h.csv", row.names = 1) ## Erroneous file name
  )
  GSE1831 <- list(
    "1h" = read.csv("GSE1831/1h.csv", row.names = 1),
    "6h" = read.csv("GSE1831/6h.csv", row.names = 1),
    "24h" = read.csv("GSE1831/24h.csv", row.names = 1),
    "72h" = read.csv("GSE1831/72h.csv", row.names = 1),
    "240h" = read.csv("GSE1831/144h.csv", row.names = 1) ## Erroneous file name
  )
  GSE10923 <- list(
    "24h" = read.csv("GSE10923/differential expression.csv", row.names = 1)
  )
  GSE1833 <- list(
    "240h" = read.csv("GSE1833/NEN.csv", row.names = 1)
  )
  GSE136913 <- list(
    "168h" = read.csv("GSE136913/GSE136913_IsoformLevelDifExpResults.csv", row.names = 1)
  )
  GSE122228 <- list(
    "8h" = read.csv("GSE122228/8h.csv", row.names = 1),
    "336h" = read.csv("GSE122228/14d.csv", row.names = 1)
  )
  GSE88992 <- list(
    "6h" = read.csv("GSE88992/6h.csv", row.names = 1),
    "12h" = read.csv("GSE88992/12h.csv", row.names = 1),
    "24h" = read.csv("GSE88992/24h.csv", row.names = 1)
  )
  GSE73878 <- list(
    "168h" = read.csv("GSE73878/7d CL.csv", row.names = 1)
  )
  GSE213393 <- list(
    "72h" = read.csv("GSE213393/C3d.csv", row.names = 1),
    "168h" = read.csv("GSE213393/C7d.csv", row.names = 1),
    "336h" = read.csv("GSE213393/C14d.csv", row.names = 1)
  )
  GSE138100 <- list(
    "24h" = list(
      "Neuronal" = read.csv("GSE138100/Neuron.csv", row.names = 1),
      "Non-neuronal" = read.csv("GSE138100/Glial.csv", row.names = 1)
    )
  )
  GSE49030 <- list(
    "1h" = read.csv("GSE49030/1h.csv", row.names = 1),
    "4h" = list(
      " " = read.csv("GSE49030/4h.csv", row.names = 1),
      "Hybrid Strain" = read.csv("GSE49030/4h Alt.csv", row.names = 1)
    ),
    "8h" = read.csv("GSE49030/8h.csv", row.names = 1),
    "24h" = read.csv("GSE49030/24h.csv", row.names = 1)
  )
  GSE79129 <- list(
    "72h" = read.csv("GSE79129/72h Differential Expression.csv", row.names = 1)
  )
  
  ##------------------------------------------------------------------------------
  
  ## Initialise list of datasets over which to iterate and apply functions
  datasets <- paste0("GSE", c("47752", "1834", "1831", "1833", "10923", "136913"))
  Rn.datasetlist <- list()
  for (set in datasets) {
    tmpset <- eval(parse(text = set))
    for (subset in 1:length(tmpset)) {
      Rn.datasetlist[[set]][names(tmpset[subset])] <- tmpset[subset]
    }
  }
  
  datasets <- paste0("GSE", c("122228", "88992", "73878", "213393", "138100",
                              "49030", "79129"))
  Mm.datasetlist <- list()
  for (set in datasets) {
    tmpset <- eval(parse(text = set))
    for (subset in 1:length(tmpset)) {
      Mm.datasetlist[[set]][names(tmpset[subset])] <- tmpset[subset]
    }
  }
  
  ##------------------------------------------------------------------------------
  
  ## Extract a dataframe containing the genes of interest from all of the datasets
  interesting.genes <- data.frame(row.names = c("symbol", "logfc", "padj",
                                                "group", "time", "dataset"))
  genes <- x ## Default is all rat and mouse genes
  
  ## Rattus norvegicus
  ndataset <- 1:length(Rn.datasetlist)
  for(i in ndataset) {
    set <- Rn.datasetlist[[i]]
    print(names(Rn.datasetlist)[[i]])
    ntimepoints <- 1:length(set)
    for(t in ntimepoints) {
      time <- set[[t]]
      if(!is.data.frame(time)) {
        nsubsets <- 1:length(time)
        for(x in nsubsets) {
          data <- time[[x]]
          names(Rn.datasetlist[[i]][[t]][[x]])[names(data) == "SYMBOL" | names(data) == "symbol" | names(data) == "Symbol"] <- "symbol"
          names(Rn.datasetlist[[i]][[t]][[x]])[names(data) == "logFC" | names(data) == "log2FoldChange"] <- "logfc"
          names(Rn.datasetlist[[i]][[t]][[x]])[names(data) == "adj.P.Val" | names(data) == "padj"] <- "padj"
          Rn.datasetlist[[i]][[t]][[x]]$group <- names(time)[x]
          Rn.datasetlist[[i]][[t]][[x]]$time <- names(set)[t]
          Rn.datasetlist[[i]][[t]][[x]]$species <- "rat"
          Rn.datasetlist[[i]][[t]][[x]]$dataset <- paste(names(time)[x],
                                                         names(Rn.datasetlist)[i], " ")
          data <- Rn.datasetlist[[i]][[t]][[x]]
          DEset <- Rn.datasetlist[[i]][[t]][[x]][data$symbol %in% genes,
                                                 c("symbol", "logfc", "padj", "species",
                                                   "group", "time", "dataset")]
          interesting.genes <- rbind(interesting.genes, DEset)
        }
      } else {
        data <- time
        names(Rn.datasetlist[[i]][[t]])[names(data) == "SYMBOL" | names(data) == "symbol" | names(data) == "Symbol"] <- "symbol"
        names(Rn.datasetlist[[i]][[t]])[names(data) == "logFC" | names(data) == "log2FoldChange"] <- "logfc"
        names(Rn.datasetlist[[i]][[t]])[names(data) == "adj.P.Val" | names(data) == "padj"] <- "padj"
        Rn.datasetlist[[i]][[t]]$time <- names(set)[t]
        Rn.datasetlist[[i]][[t]]$species <- "rat"
        Rn.datasetlist[[i]][[t]]$dataset <- paste(names(Rn.datasetlist)[i], " ")
        Rn.datasetlist[[i]][[t]]$group <- NA
        data <- Rn.datasetlist[[i]][[t]]
        DEset <- Rn.datasetlist[[i]][[t]][data$symbol %in% genes,
                                          c("symbol", "logfc", "padj", "species",
                                            "group", "time", "dataset")]
        interesting.genes <- rbind(interesting.genes, DEset)
      }
    }
  }
  
  ## Mus musculus
  ndataset <- 1:length(Mm.datasetlist)
  for(i in ndataset) {
    set <- Mm.datasetlist[[i]]
    ntimepoints <- 1:length(set)
    for(t in ntimepoints) {
      time <- set[[t]]
      if(!is.data.frame(time)) {
        nsubsets <- 1:length(time)
        for(x in nsubsets) {
          data <- time[[x]]
          names(Mm.datasetlist[[i]][[t]][[x]])[names(data) == "SYMBOL" | names(data) == "symbol" | names(data) == "Symbol"] <- "symbol"
          names(Mm.datasetlist[[i]][[t]][[x]])[names(data) == "logFC" | names(data) == "log2FoldChange"] <- "logfc"
          names(Mm.datasetlist[[i]][[t]][[x]])[names(data) == "adj.P.Val" | names(data) == "padj"] <- "padj"
          Mm.datasetlist[[i]][[t]][[x]]$group <- names(time)[x]
          Mm.datasetlist[[i]][[t]][[x]]$time <- names(set)[t]
          Mm.datasetlist[[i]][[t]][[x]]$species <- "mouse"
          if (names(time)[x] != " ") {
            Mm.datasetlist[[i]][[t]][[x]]$dataset <- paste(names(time)[x],
                                                           names(Mm.datasetlist)[i], " ")
          } else {
            Mm.datasetlist[[i]][[t]][[x]]$dataset <- paste(names(Mm.datasetlist)[i], " ")
          }
          data <- Mm.datasetlist[[i]][[t]][[x]]
          DEset <- Mm.datasetlist[[i]][[t]][[x]][data$symbol %in% genes,
                                                 c("symbol", "logfc", "padj", "species",
                                                   "group", "time", "dataset")]
          interesting.genes <- rbind(interesting.genes, DEset)
        }
      } else {
        data <- time
        names(Mm.datasetlist[[i]][[t]])[names(data) == "SYMBOL" | names(data) == "symbol" | names(data) == "Symbol"] <- "symbol"
        names(Mm.datasetlist[[i]][[t]])[names(data) == "logFC" | names(data) == "log2FoldChange"] <- "logfc"
        names(Mm.datasetlist[[i]][[t]])[names(data) == "adj.P.Val" | names(data) == "padj"] <- "padj"
        Mm.datasetlist[[i]][[t]]$time <- names(set)[t]
        Mm.datasetlist[[i]][[t]]$species <- "mouse"
        Mm.datasetlist[[i]][[t]]$dataset <- paste(names(Mm.datasetlist)[i], " ")
        Mm.datasetlist[[i]][[t]]$group <- NA
        data <- Mm.datasetlist[[i]][[t]]
        DEset <- Mm.datasetlist[[i]][[t]][data$symbol %in% genes,
                                          c("symbol", "logfc", "padj", "species",
                                            "group", "time", "dataset")]
        interesting.genes <- rbind(interesting.genes, DEset)
      }
    }
  }
  print(interesting.genes)
  return(interesting.genes)
}

##------------------------------------------------------------------------------

## Plot genes involved in GABAergic signalling and chloride homeostasis
## Get all genes related to the relevant GO terms
m <- AnnotationDbi::select(org.Mm.eg.db, 
                               c("GO:0055064", "GO:0051932", "GO:0160001"), 
                               "SYMBOL", "GO")$SYMBOL ## Mouse
r <- AnnotationDbi::select(org.Rn.eg.db, 
                               c("GO:0055064", "GO:0051932", "GO:0160001"), 
                               "SYMBOL", "GO")$SYMBOL ## Rat
genes <- na.omit(unique(c(m, r)))

## Get log fold change of selected genes
interesting.genes <- genesofinterest(c(genes, "Slc6a1", "Slc6a11", "Slc6a12", "Slc6a13")) ## No GABA reuptake GO term for rats and mice on AmiGO, 
                                                                                          ## so instead added manually

## Order time labels in chronological order
interesting.genes$time <- factor(interesting.genes$time, levels = c("1h", "4h", "6h",
                                                                    "8h", "12h", "24h",
                                                                    "72h", "168h", 
                                                                    "240h", "336h"))
## Convert data to wide to be used for clustering of genes
interesting.genes$id <- paste0(interesting.genes$time,interesting.genes$dataset) ## A unique identifier of a gene so that occurrences in other time points or
                                                                                 ## or datasets do not collapse into the same entry in the list
wide <- pivot_wider(interesting.genes, id_cols = "id", names_from = "symbol", values_from = "logfc",
                    values_fill = 0, values_fn = function(x) {mean(x)})
matrix <- t(as.matrix(wide))
colnames(matrix) <- matrix[1,]
matrix <- matrix[-1,]
distances <- dist(matrix)
heatmap <- pheatmap(as.matrix(distances),
                    clustering_distance_rows = distances,
                    clustering_distance_cols = distances)

## Arrange the gene list such that it aligns with the clustering          
clusterorders <- partition_leaves(as.dendrogram(heatmap[[2]]))[[1]]
interesting.genes$symbol <- factor(interesting.genes$symbol,
                                   levels = clusterorders)

## Unit conversions for use with the labeller parameter of facet_grid()
hour.to.day <- list('1h' = "1h",
                    '4h' = "4h",
                    '6h' = "6h",
                    '8h' = "8h",
                    '12h' = "12h",
                    '24h' = "24h",
                    '72h' = "3d",
                    '168h' = "7d",
                    '240h' = "10d",
                    '336h' = "14d")
hours.to.days <- function(variable,value){
  return(hour.to.day[value])
}

## Plot expression
## scale_fill_gradientn() works on a scale of [0,1] and thus a conversion is needed to find log2(0) (no differential expression)
## on the new scale
min <- min(interesting.genes$logfc)
max <- max(interesting.genes$logfc)                    
## Formula for rescaling [a b] to [c d] is f(x) = (x-a)*(d-c)/(b-a) + c, hence
## the 0 value of the colour gradient is (-min)/(max-min): 
## [min, max] -> [0, 1]; x = 0 on original scale 
ggplot(interesting.genes, aes(x = dataset, y = symbol, fill = logfc)) + 
  geom_tile() + scale_fill_gradientn(values=c(1, (-min)/(max-min), 0), colours=c("magenta", "white", "yellow")) +
  facet_grid(~time, scales = "free", space = "free",  labeller = hours.to.days,
             shrink = T) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        plot.margin = unit(c(1,1,1,1), "cm"), axis.text = element_text(size = 10),
        axis.title.x = element_blank(),
        panel.background = element_rect(fill = 'white', colour = 'grey')) +
  ylab("Gene")

## Plot dendrogram
par(mar=c(2,2,2,5))
colors <- heatmap[[2]] %>%
  as.dendrogram() %>%
  set("branches_k_color", k=7) %>%
  set("labels_colors", k=7) %>% 
  labels_colors() 
edit <- c("#CC476B", "#228B00", "#9F7000", "#FF7F00", "#009681", "#0082CE", "#C03FBE") ## Replacement for colours
match <- cbind(unique(colors), edit) 
a <- unlist(lapply(colors, function(x) {match[,2][match[,1]==x]})) ## Each gene is assigned a new colour based on the corresponding replacement
heatmap[[2]] %>%
  as.dendrogram() %>%
  set("branches_k_col", value=unique(a), k = 7) %>%
  set("labels_colors", a) %>%  
  set("branches_lwd", 2) %>%  
  plot(horiz=T, lwd=2)

##------------------------------------------------------------------------------

## Plot pro- and anti-apoptotic genes
## Get all genes related to the relevant GO terms
m.pro <- AnnotationDbi::select(org.Mm.eg.db, 
                           c("GO:0043525"), 
                           "SYMBOL", "GO")$SYMBOL ## Mouse
m.anti <- AnnotationDbi::select(org.Mm.eg.db, 
                               c("GO:0043525"), 
                               "SYMBOL", "GO")$SYMBOL ## Mouse
r.pro <- AnnotationDbi::select(org.Rn.eg.db, 
                               c("GO:0043524"), 
                               "SYMBOL", "GO")$SYMBOL ## Rat
r.anti <- AnnotationDbi::select(org.Rn.eg.db, 
                                c("GO:0043524"), 
                                "SYMBOL", "GO")$SYMBOL ## Rat
genes <- unique(c(m.pro,r.pro,m.anti,r.anti))

## Get log fold change of selected genes
interesting.genes <- genesofinterest(genes)

## Assign genes a pro-apoptotic or anti-apoptotic label
interesting.genes$apo <- ""
interesting.genes$apo[interesting.genes$symbol %in% c(m.pro, r.pro)] <- "Pro-apoptotic"
interesting.genes$apo[interesting.genes$symbol %in% c(m.anti, r.anti)] <- "Anti-apoptotic"
interesting.genes$apo[interesting.genes$symbol %in% c(m.anti, r.anti) & interesting.genes$symbol %in% c(m.pro, r.pro)] <- "Both"

## Order time labels in chronological order
interesting.genes$time <- factor(interesting.genes$time, levels = c("1h", "4h", "6h",
                                                                    "8h", "12h", "24h",
                                                                    "72h", "168h", 
                                                                    "240h", "336h"))

## Convert data to wide to be used for clustering
interesting.genes$id <- paste0(interesting.genes$time,interesting.genes$dataset) ## A unique identifier of a gene so that occurrences in other time points or
                                                                                 ## or datasets do not collapse into the same entry in the list
wide <- pivot_wider(interesting.genes, id_cols = "id", names_from = "symbol", values_from = "logfc",
                    values_fill = 0, values_fn = function(x) {mean(x)})
matrix <- t(as.matrix(wide))
colnames(matrix) <- matrix[1,]
matrix <- matrix[-1,]
distances <- dist(matrix)
heatmap <- pheatmap(as.matrix(distances),
                    clustering_distance_rows = distances,
                    clustering_distance_cols = distances)

## Arrange the gene list such that it aligns with the clustering  
clusterorders <- partition_leaves(as.dendrogram(heatmap[[2]]))[[1]]
interesting.genes$symbol <- factor(interesting.genes$symbol,
                                   levels = clusterorders)

## Select only clusters of genes significantly differentially expressed across multiple datasets
interesting.genes <- interesting.genes[interesting.genes$symbol %in% 
                                         clusterorders[c(1:24,232:272)],]

## Unit conversions for use with the labeller parameter of facet_grid()
## Also a dictionary for the apoptosis labels required for the function to work
hour.to.day <- list('1h' = "1h",
                    '4h' = "4h",
                    '6h' = "6h",
                    '8h' = "8h",
                    '12h' = "12h",
                    '24h' = "24h",
                    '72h' = "3d",
                    '168h' = "7d",
                    '240h' = "10d",
                    '336h' = "14d",
                    'Pro-apoptotic' = "Pro-apoptotic",
                    'Anti-apoptotic' = "Anti-apoptotic",
                    'Both' = "Both")
hours.to.days <- function(variable,value){
  return(hour.to.day[value])
}

## Plot expression
## scale_fill_gradientn() works on a scale of [0,1] and thus a conversion is needed to find log2(0) (no differential expression)
## on the new scale
min <- min(interesting.genes$logfc)
max <- max(interesting.genes$logfc)
## Formula for rescaling [a b] to [c d] is f(x) = (x-a)*(d-c)/(b-a) + c, hence
## the 0 value of the colour gradient is (-min)/(max-min): 
## [min, max] -> [0, 1]; x = 0 on original scale 
ggplot(interesting.genes, aes(x = dataset, y = symbol, fill = logfc)) + 
  geom_tile() + scale_fill_gradientn(values=c(1, (-min)/(max-min), 0), colours=c("magenta", "white", "yellow")) +
  facet_grid(apo~time, scales = "free", space = "free",  labeller = hours.to.days,
             shrink = T) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        plot.margin = unit(c(1,1,1,1), "cm"), axis.text = element_text(size = 10),
        axis.title.x = element_blank(),
        panel.background = element_rect(fill = 'white', colour = 'grey')) +
  ylab("Gene")
