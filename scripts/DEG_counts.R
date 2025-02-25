## Load libraries
library(dplyr)
library(ggplot2)
library(GO.db)
library(AnnotationDbi)
library(org.Mm.eg.db)
library(org.Rn.eg.db)

##---------------------------------------------------------------------------

## Load DEGs from files
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

##-----------------------------------------------------------------------------------------------

## Plotting the number of DEGs
## Count the number of DEGs
DEgenes <- rbind(DEgenes.mouse, DEgenes.rat)
DEgenes$unique <- paste(DEgenes$entrez, DEgenes$dataset, DEgenes$time) ## A unique identifier so that only duplicate genes at the same time in the same
                                                                       ## dataset are deleted
DEgenes <- DEgenes[!duplicated(DEgenes$unique),] ## Remove duplicate genes
gene.counts <- dplyr::rename(dplyr::count(DEgenes[,4:5], time, dataset), Genes = n)

## Add padding to separate early time points from later time points visually
## and order time lables in chronological order
padding <- data.frame(time = c(" ", "  ", "   "),
                      dataset = c("padding 1", "padding 2", "padding 3"),
                      Genes = c(0,0,0))
gene.counts <- rbind(gene.counts,padding)
gene.counts$time <- factor(gene.counts$time, levels = c("1h", "4h", "6h",
                                                        "8h", "12h", "24h", 
                                                        " ", "  ", "   ",
                                                        "3d", "7d",
                                                        "10d", "14d"))
## Arrange datasets                                                     
gene.counts$dataset <- factor(gene.counts$dataset, 
                              levels = paste(c("Ipsilateral GSE73878", "Contralateral GSE73878",
                                               "Ipsilateral GSE213393", "Contralateral GSE213393",
                                               "GSE1831", "GSE1834", 
                                               "GSE10923","GSE79129", "GSE88992",
                                               "GSE49030", "Hybrid Strain GSE49030",
                                               "Neuronal GSE138100", "Non-neuronal GSE138100",
                                               "Nadler GSE47752", "Whadman GSE47752", 
                                               "GSE122228", "GSE1833"), " "))
                                               
## Plot number of DEGs
gene.counts$Genes[40:42] <- 6500 ## Dummy numbers for padding
ggplot(gene.counts, aes(fill=dataset, y=Genes, x=time)) + 
  theme_bw() + xlab("Time") + labs(color="Dataset") +
  geom_line(aes(x = time, y = Genes, color = dataset, group = dataset), 
            size = 1) + geom_point(aes(x = time, y = Genes, color = dataset, group = dataset), 
                                   size = 2, show.legend = F)

##----------------------------------------------------------------------------------------

## Plot number of immune-related DEGs
## Get all genes related to immune functions 
go <- c(unlist(get("GO:0002376", GOBPCHILDREN)), 
        unlist(get("GO:0016032", GOBPCHILDREN)),
        unlist(get("GO:0001816", GOBPCHILDREN))) 
prevgo <- c() 
while(length(prevgo)!=length(go)) {
  prevgo <- unique(unname(go))
  go <- c(go, unlist(lapply(unlist(go), function(x) get(x, GOBPCHILDREN))))
  go <- unique(unname(na.omit(go)))
} ## This routine takes a while
g <- unique(c(AnnotationDbi::select(org.Mm.eg.db, 
                                    c(go), 
                                    "ENTREZID", "GO")$ENTREZID,
              AnnotationDbi::select(org.Rn.eg.db, 
                                    c(go), 
                                    "ENTREZID", "GO")$ENTREZID))
DEgenes <- rbind(DEgenes.mouse, DEgenes.rat)
DEgenes <- DEgenes[DEgenes$entrez %in% g,]

DEgenes$unique <- paste(DEgenes$entrez, DEgenes$dataset, DEgenes$time) ## A unique identifier so that only duplicate genes at the same time in the same
                                                                       ## dataset are deleted
DEgenes <- DEgenes[!duplicated(DEgenes$unique),] ## Remove duplicate genes
gene.counts <- dplyr::rename(dplyr::count(DEgenes[,4:5], time, dataset), Genes = n)

## Add padding to separate early time points from later time points visually
## and order time lables in chronological order
padding <- data.frame(time = c(" ", "  ", "   "),
                      dataset = c("padding 1", "padding 2", "padding 3"),
                      Genes = c(0,0,0))
gene.counts <- rbind(gene.counts,padding)
gene.counts$time <- factor(gene.counts$time, levels = c("1h", "4h", "6h",
                                                        "8h", "12h", "24h", 
                                                        " ", "  ", "   ",
                                                        "3d", "7d",
                                                        "10d", "14d"))
## Arrange datasets                                                     
gene.counts$dataset <- factor(gene.counts$dataset, 
                              levels = paste(c("Ipsilateral GSE73878", "Contralateral GSE73878",
                                               "Ipsilateral GSE213393", "Contralateral GSE213393",
                                               "GSE1831", "GSE1834", 
                                               "GSE10923","GSE79129", "GSE88992",
                                               "GSE49030", "Hybrid Strain GSE49030",
                                               "Neuronal GSE138100", "Non-neuronal GSE138100",
                                               "Nadler GSE47752", "Whadman GSE47752", 
                                               "GSE122228", "GSE1833"), " "))
                                               
## Plot number of DEGs
gene.counts$Genes[40:42] <- 900 ## Dummy numbers for padding
ggplot(gene.counts, aes(fill=dataset, y=Genes, x=time)) + 
  theme_bw() + xlab("Time") + labs(fill="Dataset") +
  geom_bar(position=position_dodge2(width = 0.9, preserve = "single"), 
           stat="identity")

##------------------------------------------------------------------------------------------------

## Counting DEG overlaps
DEgenes <- rbind(DEgenes.mouse, DEgenes.rat)
DEgenes$unique <- paste(DEgenes$entrez, DEgenes$dataset, DEgenes$time) ## A unique identifier so that only duplicate genes at the same time in the same
## dataset are deleted
DEgenes <- DEgenes[!duplicated(DEgenes$unique),] ## Remove duplicated genes
down <- DEgenes[DEgenes$logfc<0,]
up <- DEgenes[DEgenes$logfc>0,]

reg <- c("up", "down")

for (regulation in reg) {
  df <- eval(parse(text=regulation))
  for (time in unique(DEgenes$time)) {
    time.point <- time
    ndatasets <- length(unique(DEgenes$dataset[DEgenes$time==time]))/2
                        if (ndatasets < 1) next
                        if (ndatasets == 1) ndatasets <- 2
                        ndatasets <- ceiling(ndatasets)
                        
                        ## Filter and group the data
                        result <- df %>%
                          filter(time == time.point) %>%
                          group_by(entrez) %>%
                          summarise(dataset_count = n_distinct(dataset)) %>%
                          filter(dataset_count >= ndatasets) %>%
                          nrow()
                        
                        unique_genes_count <- df %>%
                          filter(time == time.point) %>%
                          summarise(unique_genes = n_distinct(entrez)) %>%
                          pull(unique_genes)
                        
                        ## Print the result
                        print(paste(regulation, time.point))
                              print(unique_genes_count)
                              print(result)
                        cat("\n")
    ## The percentage of overlapping DEGs were graphed with Microsoft Excel
  }
}
