## Load required packages
library(clusterProfiler)
library(ggplot2)
library(GO.db)
library(org.Mm.eg.db)
library(org.Rn.eg.db)
library(tidyverse)

##----------------------------------------------------------------------------------

## Mus musculus
## Load in DE data for each dataset at each time point within the latency period
GSE122228 <- list(
  "8h" = read.csv("GSE122228/8h.csv", row.names = 1),
  "14d" = read.csv("GSE122228/14d.csv", row.names = 1)
)
GSE88992 <- list(
  "6h" = read.csv("GSE88992/6h.csv", row.names = 1),
  "12h" = read.csv("GSE88992/12h.csv", row.names = 1),
  "24h" = read.csv("GSE88992/24h.csv", row.names = 1)
)
GSE73878 <- list(
  "7d" = list(
    "Contralateral" = read.csv("GSE73878/7d CL.csv", row.names = 1),
    "Ipsilateral" = read.csv("GSE73878/7d IL.csv", row.names = 1)
  )
)
GSE213393 <- list(
  "3d" = list(
    "Contralateral" = read.csv("GSE213393/C3d.csv", row.names = 1),
    "Ipsilateral" = read.csv("GSE213393/I3d.csv", row.names = 1)
    ),
  "7d" = list(
    "Contralateral" = read.csv("GSE213393/C7d.csv", row.names = 1),
    "Ipsilateral" = read.csv("GSE213393/I7d.csv", row.names = 1)
  ),
  "14d" = list(
    "Contralateral" = read.csv("GSE213393/C14d.csv", row.names = 1),
    "Ipsilateral" = read.csv("GSE213393/I14d.csv", row.names = 1)
  )
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
  "3d" = read.csv("GSE79129/72h Differential Expression.csv", row.names = 1)
)

##------------------------------------------------------------------------------

## Initialise list of datasets over which to iterate and apply functions
datasets <- paste0("GSE", c("122228", "88992", "73878", "213393", "138100",
                            "49030", "79129"))
datasetlist <- list()
for (set in datasets) {
  tmpset <- eval(parse(text = set))
  for (subset in 1:length(tmpset)) {
    datasetlist[[set]][names(tmpset[subset])] <- tmpset[subset]
  }
}

##------------------------------------------------------------------------------

## Add Entrez IDs and split upregulated and downregulated genes to prepare for
## downstream GO analysis
upregulated <- data.frame(row.names = c("entrez", "logfc", "group", "time", "dataset"))
downregulated <- data.frame(row.names = c("entrez", "logfc", "group", "time", "dataset"))
ndataset <- 1:length(datasetlist)
for(i in ndataset) {
  set <- datasetlist[[i]]
  ntimepoints <- 1:length(set)
  for(t in ntimepoints) {
    time <- set[[t]]
    if(!is.data.frame(time)) {
      nsubsets <- 1:length(time)
      for(x in nsubsets) {
        data <- time[[x]]
        names(data)[names(data) == "SYMBOL" | names(data) == "symbol" | names(data) == "Symbol"] <- "symbol"
        names(datasetlist[[i]][[t]][[x]])[names(data) == "logFC" | names(data) == "log2FoldChange"] <- "logfc"
        match <- AnnotationDbi::select(org.Mm.eg.db,
                                       keys = data$symbol,
                                       columns = c("ENTREZID", "SYMBOL"),
                                       keytype = "SYMBOL")
        datasetlist[[i]][[t]][[x]]$entrez <- match$ENTREZID[match(data$symbol, match$SYMBOL)]
        datasetlist[[i]][[t]][[x]]$group <- names(time)[x]
        datasetlist[[i]][[t]][[x]]$time <- names(set)[t]
        if (names(time)[x] != " ") {
          datasetlist[[i]][[t]][[x]]$dataset <- paste(names(time)[x],
                                                      names(datasetlist)[i], " ")
        } else {
          datasetlist[[i]][[t]][[x]]$dataset <- paste(names(datasetlist)[i], " ")
        }
        DEset <- datasetlist[[i]][[t]][[x]][,c("entrez", "logfc", "group", "time", "dataset")]
        upregulated <- rbind(upregulated, DEset[DEset$logfc > 0.584,])
        downregulated <- rbind(downregulated, DEset[DEset$logfc < -0.584,])
      }
    } else {
      data <- time
      names(data)[names(data) == "SYMBOL" | names(data) == "symbol" | names(data) == "Symbol"] <- "symbol"
      names(datasetlist[[i]][[t]])[names(data) == "logFC" | names(data) == "log2FoldChange"] <- "logfc"
      match <- AnnotationDbi::select(org.Mm.eg.db,
                                     keys = data$symbol,
                                     columns = c("ENTREZID", "SYMBOL"),
                                     keytype = "SYMBOL")
      datasetlist[[i]][[t]]$entrez <- match$ENTREZID[match(data$symbol, match$SYMBOL)]
      datasetlist[[i]][[t]]$time <- names(set)[t]
      datasetlist[[i]][[t]]$dataset <- paste(names(datasetlist)[i], " ")
      datasetlist[[i]][[t]]$group <- NA
      DEset <- datasetlist[[i]][[t]][,c("entrez", "logfc", "group", "time", "dataset")]
      upregulated <- rbind(upregulated, DEset[DEset$logfc > 0.584,])
      downregulated <- rbind(downregulated, DEset[DEset$logfc < -0.584,])
    }
  }
}
DEgenes.mouse <- rbind(upregulated, downregulated)

##------------------------------------------------------------------------------

## Run the GO enrichment functions, getting rid of similar GO terms under a more
## general term and plot
upGO <- compareCluster(entrez~time+dataset, fun = "enrichGO",
                       data = upregulated, ont = "BP", pvalueCutoff = 0.05,
                       OrgDb = org.Mm.eg.db, readable = T)

downGO <- compareCluster(entrez~time+dataset, fun = "enrichGO",
                         data = downregulated, ont = "BP", pvalueCutoff = 0.05,
                         OrgDb = org.Mm.eg.db, readable = T)

go <- c(unlist(get("GO:0002376", GOBPCHILDREN)), ## Immune system process
        unlist(get("GO:0016032", GOBPCHILDREN)), ## Viral process
        unlist(get("GO:0001816", GOBPCHILDREN))) ## Cytokine production
prevgo <- c() 
while(length(prevgo)!=length(go)) {
  prevgo <- unique(unname(go))
  go <- c(go, unlist(lapply(unlist(go), function(x) {get(x, GOBPCHILDREN)})))
  go <- unique(unname(na.omit(go)))
  print(length(prevgo))
  print(length(go))
} ## This process takes a while
immune <- unique(unname(go))

## Filter out non-immune terms
GOtermdf <- data.frame("ID" = upGO@compareClusterResult$ID, 
                       "name" = upGO@compareClusterResult$Description)
nonimmune <- GOtermdf$ID[! GOtermdf$ID %in% immune]
upGO.immune <- dropGO(upGO, term = nonimmune)
downGO.immune <- dropGO(downGO, term = nonimmune)

## Simplify
immuneup <- clusterProfiler::simplify(upGO.immune)
immunedown <- clusterProfiler::simplify(downGO.immune)

## Create the dotplot function
clusteredGO <- c("immuneup", "immunedown")
names(clusteredGO) <- c("Enriched Immune-related GO terms in upregulated genes (Mouse)", 
                        "Enriched Immune-related GO terms in downregulated genes (Mouse)")
wrapper <- function(x, ...) ## For fitting the GO term descriptions on the y-axis
{
  lapply(x, function(x){paste(strwrap(x, ...), collapse = "\n")})
}
plot.go <- function(go, title) {
  dotplot(go, x="dataset", label_format = function(x){
    ifelse(nchar(x) > 59, 
           wrapper(x, width = 47) %>% strtrim(width = 60) %>% paste0("..."), ## Cut off very long descriptions
           wrapper(x, width = 47))}) + 
    facet_grid(~factor(time, levels = c("1h", "4h", "6h", "8h", 
                                        "12h", "24h", "3d", "7d", "14d")),
               scales = "free", space = "free") + 
    ggtitle(wrapper(title, width = 47)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          plot.margin = unit(c(1,1,1,1), "cm"), axis.text = element_text(size = 10),
          axis.title.x = element_blank())    
}

## Plot enrichment dotplots
lapply(seq_along(clusteredGO), function(x,y,i) {plot.go(eval(parse(text=x[i])), y[i])},
       x=clusteredGO, y=names(clusteredGO))

##-----------------------------------------------------------------------------------

## Rattus norvegicus
## Load in DE data for each dataset at each time point within the latency period
GSE47752 <- list(
  "24h" = list(
    "Nadler" = read.csv("GSE47752/Nadler 1 Day.csv", row.names = 1),
    "Whadman" = read.csv("GSE47752/Whadman 1 Day.csv", row.names = 1)
  ),
  "3d" = list(
    "Nadler" = read.csv("GSE47752/Nadler 3 Days.csv", row.names = 1),
    "Whadman" = read.csv("GSE47752/Whadman 3 Days.csv", row.names = 1)
  ),
  "10d" = list(
    "Nadler" = read.csv("GSE47752/Nadler 10 Days.csv", row.names = 1),
    "Whadman" = read.csv("GSE47752/Whadman 10 Days.csv", row.names = 1)
  )
)
GSE1834 <- list(
  "1h" = read.csv("GSE1834/1h.csv", row.names = 1),
  "6h" = read.csv("GSE1834/6h.csv", row.names = 1),
  "24h" = read.csv("GSE1834/24h.csv", row.names = 1),
  "3d" = read.csv("GSE1834/72h.csv", row.names = 1),
  "10d" = read.csv("GSE1834/144h.csv", row.names = 1)
)
GSE1831 <- list(
  "1h" = read.csv("GSE1831/1h.csv", row.names = 1),
  "6h" = read.csv("GSE1831/6h.csv", row.names = 1),
  "24h" = read.csv("GSE1831/24h.csv", row.names = 1),
  "3d" = read.csv("GSE1831/72h.csv", row.names = 1),
  "10d" = read.csv("GSE1831/144h.csv", row.names = 1)
)
GSE10923 <- list(
  "24h" = read.csv("GSE10923/differential expression.csv", row.names = 1)
)
GSE1833 <- list(
  "10d" = read.csv("GSE1833/NEN.csv", row.names = 1)
  )
GSE136913 <- list(
  "7d" = read.csv("GSE136913/GSE136913_IsoformLevelDifExpResults.csv", row.names = 1)
)


##------------------------------------------------------------------------------

## Initialise list of datasets over which to iterate and apply functions
datasets <- paste0("GSE", c("47752", "1834", "1831", "1833", "10923", "136913"))
datasetlist <- list()
for (set in datasets) {
  tmpset <- eval(parse(text = set))
  for (subset in 1:length(tmpset)) {
    datasetlist[[set]][names(tmpset[subset])] <- tmpset[subset]
  }
}

##------------------------------------------------------------------------------

## Add Entrez IDs and split upregulated and downregulated genes to prepare for
## downstream GO analysis
upregulated <- data.frame(row.names = c("entrez", "logfc", "group", "time", "dataset"))
downregulated <- data.frame(row.names = c("entrez", "logfc", "group", "time", "dataset"))
ndataset <- 1:length(datasetlist)
for(i in ndataset) {
  set <- datasetlist[[i]]
  ntimepoints <- 1:length(set)
  for(t in ntimepoints) {
    time <- set[[t]]
    if(!is.data.frame(time)) {
      nsubsets <- 1:length(time)
      for(x in nsubsets) {
        data <- time[[x]]
        names(data)[names(data) == "SYMBOL" | names(data) == "symbol" | names(data) == "Symbol"] <- "symbol"
        names(datasetlist[[i]][[t]][[x]])[names(data) == "logFC" | names(data) == "log2FoldChange"] <- "logfc"
        match <- AnnotationDbi::select(org.Rn.eg.db,
                                       keys = data$symbol,
                                       columns = c("ENTREZID", "SYMBOL"),
                                       keytype = "SYMBOL")
        datasetlist[[i]][[t]][[x]]$entrez <- match$ENTREZID[match(data$symbol, match$SYMBOL)]
        datasetlist[[i]][[t]][[x]]$group <- names(time)[x]
        datasetlist[[i]][[t]][[x]]$time <- names(set)[t]
        datasetlist[[i]][[t]][[x]]$dataset <- paste(names(time)[x],
                                                    names(datasetlist)[i], " ")
        DEset <- datasetlist[[i]][[t]][[x]][,c("entrez", "logfc", "group", "time", "dataset")]
        upregulated <- rbind(upregulated, DEset[DEset$logfc > 0.584,])
        downregulated <- rbind(downregulated, DEset[DEset$logfc < -0.584,])
      }
    } else {
      data <- time
      names(data)[names(data) == "SYMBOL" | names(data) == "symbol" | names(data) == "Symbol"] <- "symbol"
      names(datasetlist[[i]][[t]])[names(data) == "logFC" | names(data) == "log2FoldChange"] <- "logfc"
      match <- AnnotationDbi::select(org.Rn.eg.db,
                                     keys = data$symbol,
                                     columns = c("ENTREZID", "SYMBOL"),
                                     keytype = "SYMBOL")
      datasetlist[[i]][[t]]$entrez <- match$ENTREZID[match(data$symbol, match$SYMBOL)]
      datasetlist[[i]][[t]]$time <- names(set)[t]
      datasetlist[[i]][[t]]$dataset <- paste(names(datasetlist)[i], " ")
      datasetlist[[i]][[t]]$group <- NA
      DEset <- datasetlist[[i]][[t]][,c("entrez", "logfc", "group", "time", "dataset")]
      upregulated <- rbind(upregulated, DEset[DEset$logfc > 0.584,])
      downregulated <- rbind(downregulated, DEset[DEset$logfc < -0.584,])
    }
  }
}
DEgenes.rat <- rbind(upregulated, downregulated)

##------------------------------------------------------------------------------

## Run the GO enrichment functions, getting rid of similar GO terms under a more
## general term and plot
upGO <- compareCluster(entrez~time+dataset, fun = "enrichGO",
                       data = upregulated, ont = "BP", pvalueCutoff = 0.05,
                       OrgDb = org.Rn.eg.db, readable = T)

downGO <- compareCluster(entrez~time+dataset, fun = "enrichGO",
                         data = downregulated, ont = "BP", pvalueCutoff = 0.05,
                         OrgDb = org.Rn.eg.db, readable = T)

go <- c(unlist(get("GO:0002376", GOBPCHILDREN)), ## Immune system process
        unlist(get("GO:0016032", GOBPCHILDREN)), ## Viral process
        unlist(get("GO:0001816", GOBPCHILDREN))) ## Cytokine production
prevgo <- c() 
while(length(prevgo)!=length(go)) {
  prevgo <- unique(unname(go))
  go <- c(go, unlist(lapply(unlist(go), function(x) {get(x, GOBPCHILDREN)})))
  go <- unique(unname(na.omit(go)))
  print(length(prevgo))
  print(length(go))
} ## This process takes a while
immune <- unique(unname(go))

## Filter out non-immune terms
GOtermdf <- data.frame("ID" = upGO@compareClusterResult$ID, 
                       "name" = upGO@compareClusterResult$Description)
nonimmune <- GOtermdf$ID[! GOtermdf$ID %in% immune]
upGO.immune <- dropGO(upGO, term = nonimmune)
downGO.immune <- dropGO(downGO, term = nonimmune)

## Simplify
immuneup <- clusterProfiler::simplify(upGO.immune)
immunedown <- clusterProfiler::simplify(downGO.immune)

## Create the dotplot function
clusteredGO <- c("immuneup", "immunedown")
names(clusteredGO) <- c("Enriched Immune-related GO terms in upregulated genes (Rat)", 
                        "Enriched Immune-related GO terms in downregulated genes (Rat)")
wrapper <- function(x, ...) ## For fitting the GO term descriptions on the y-axis
{
  lapply(x, function(x){paste(strwrap(x, ...), collapse = "\n")})
}
plot.go <- function(go, title) {
  dotplot(go, x="dataset", label_format = function(x){
    ifelse(nchar(x) > 59, 
           wrapper(x, width = 47) %>% strtrim(width = 60) %>% paste0("..."), ## Cut off very long descriptions
           wrapper(x, width = 47))}) + 
    facet_grid(~factor(time, levels = c("1h", "4h", "6h", "8h", 
                                        "12h", "24h", "3d", "7d", "14d")),
               scales = "free", space = "free") + 
    ggtitle(wrapper(title, width = 47)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          plot.margin = unit(c(1,1,1,1), "cm"), axis.text = element_text(size = 10),
          axis.title.x = element_blank())    
}

## Plot enrichment dotplots
lapply(seq_along(clusteredGO), function(x,y,i) {plot.go(eval(parse(text=x[i])), y[i])},
       x=clusteredGO, y=names(clusteredGO))

##---------------------------------------------------------------------------------------------------------

## STEM clusters (Mouse)

## Load data
data <- read.table("Mouse clusters.txt", header = T, sep = "\t", ) ## STEM output
data$Symbol <- str_to_title(data$Symbol)
rownames(data) <- data$Symbol
data <- data[,-(c(1:3,5))]

## Attach EntrezIDs to the genes
clusters <- data.frame(Cluster = data[,1], row.names = rownames(data))
match <- AnnotationDbi::select(org.Mm.eg.db,
                               keys = rownames(clusters),
                               columns = c("ENTREZID", "SYMBOL"),
                               keytype = "SYMBOL")
clusters$entrez <- match$ENTREZID[match(rownames(clusters), match$SYMBOL)]
clusters <- clusters[clusters$Cluster %in% c(9,10,38,1,23,43,30,8,25),]
clusters$Cluster <- factor(clusters$Cluster, levels = c(8,9,23,25,38,30,10,43,1))

## Run GO term enrichment analysis
clustergo <- compareCluster(entrez~Cluster, fun = "enrichGO", ont = "BP", 
                            pvalueCutoff = 0.05, data = na.omit(clusters),
                            OrgDb = org.Mm.eg.db, readable = T)
simplifiedclusters <- clusterProfiler::simplify(clustergo)

## Plot enrichment dotplot
wrapper <- function(x, ...) ## For fitting the GO term descriptions on the y-axis
{
  lapply(x, function(x){paste(strwrap(x, ...), collapse = "\n")})
}
y <- dotplot(simplifiedclusters, label_format = function(x){
  ifelse(nchar(x) > 59, 
        wrapper(x, width = 47) %>% strtrim(width = 60) %>% paste0("..."), ## Cut off very long descriptions
        wrapper(x, width = 47))}) 
y

##--------------------------------------------------------------------------------------

## STEM clusters (Rat)

## Load data
data <- read.table("Rat clusters.txt", header = T, sep = "\t", ) ## STEM output
data$Symbol <- str_to_title(data$Symbol)
rownames(data) <- data$Symbol
data <- data[,-(c(1:3,5))]

## Attach EntrezIDs to genes
clusters <- data.frame(Cluster = data[,1], row.names = rownames(data))
match <- AnnotationDbi::select(org.Rn.eg.db,
                               keys = rownames(clusters),
                               columns = c("ENTREZID", "SYMBOL"),
                               keytype = "SYMBOL")
clusters$entrez <- match$ENTREZID[match(rownames(clusters), match$SYMBOL)]
clusters <- clusters[clusters$Cluster %in% c(8,10,24,27,29,47),]
clusters$Cluster <- factor(clusters$Cluster, levels = c(8,10,24,27,29,47))

## Run GO term enrichment analysis
clustergo <- compareCluster(entrez~Cluster, fun = "enrichGO", ont = "BP", 
                            pvalueCutoff = 0.05, data = na.omit(clusters),
                            OrgDb = org.Rn.eg.db, readable = T)
simplifiedclusters <- clusterProfiler::simplify(clustergo)

## Plot enrichment dotplot
wrapper <- function(x, ...) ## For fitting the GO term descriptions on the y-axis
{
  lapply(x, function(x){paste(strwrap(x, ...), collapse = "\n")}) 
}
y <- dotplot(simplifiedclusters, label_format = function(x){
  ifelse(nchar(x) > 59, 
        wrapper(x, width = 47) %>% strtrim(width = 60) %>% paste0("..."), ## Cut off very long descriptions
        wrapper(x, width = 47))}) 
y
