## Load required packages
library(affy)
library(dplyr)
library(GEOquery)
library(DESeq2)
library(SummarizedExperiment)
library(rgu34acdf)
library(rgu34a.db)
library(mouse4302.db)
library(mouse4302cdf)
library(rat2302cdf)
library(rat2302.db)

## Data was downloaded using getGEOSuppFiles()
##-----------------------------------------------------------------------------------------------

## GSE1831
gse <- getGEO("GSE1831")
names <- pData(gse[[1]])$title
cels <- list.files("GSE1831/data/", pattern = "cel") ## Data stored locally
raw.data <- ReadAffy(verbose = FALSE, 
                     filenames = paste("GSE1831/data", cels, sep = "/"),
                     cdfname = "rgu34acdf") ## Load microarray data
data.rma.norm <- affy::rma(raw.data)
gse <- exprs(data.rma.norm)
gse <- data.frame(gse)

## Make sample labels readable
colnames(gse) <- names
gse <- relocate(gse, grep("C", colnames(gse)), .before = colnames(gse)[1])
colnames(gse) <- gsub("brain, hippocampus: ", "", colnames(gse)) ## It's all hippocampus
colnames(gse) <- gsub("-", "", colnames(gse))
colnames(gse) <- gsub("h.*", "", colnames(gse))

## Run PCA
vargene <- apply(gse, 1, var)
names(vargene) <- rownames(gse)
vars <- sort(vargene,decreasing=T)
topgenes <- names(vars)[1:1000]
filtered <- gse[rownames(gse) %in% topgenes,]
PC <- prcomp(t(filtered))
varexpl <- (PC$sdev^2)/sum(PC$sdev^2)*100

## PCA plot
dev.new(width=20, height=20)
colfunc <- colorRampPalette(c("red", "yellow"))
colours <- c(rep("blue", 15), rep(colfunc(5)[1:5], each = 3))
colors <- c(
  "#000080",  # Navy Blue
  "#1E90FF",  # Dodger Blue (lighter blue)
  "#87CEEB",  # Sky Blue
  "#e1c4ff",  # Lilac
  "#9B30FF",  # Purple
  "#FF3333",  # Red
  "#FF9933",  # Light orange
  "#FFFF00"   # Bright yellow
)
df <- data.frame(Colour = colors,
                 Time = c("1h", "4h", "6h", "8h",
                          "12h", "24h", "72h","240h"))
colours <- c(rep("green", 15), rep(df$Colour[c(1,3,6,7,8)], each = 3)) ## Green = Controls
plot(PC$x[,1],PC$x[,2],cex=1.2,bg=colours,
     xlab=paste0("PC1 (", round(varexpl[1], 2), "%)"),
     ylab=paste0("PC2 (", round(varexpl[2], 2), "%)"),
     main="GSE1831",las=1, pch=21
)

##----------------------------------------------------------------------------------------------

## GSE1834
gse <- getGEO("GSE1834")
names <- pData(gse[[1]])$title
cels <- list.files("GSE1834/data/", pattern = "cel") ## Data stored locally
raw.data <- ReadAffy(verbose = FALSE, 
                     filenames = paste("GSE1834/data", cels, sep = "/"),
                     cdfname = "rgu34acdf") ## Load microarray data
data.rma.norm <- affy::rma(raw.data)
gse <- exprs(data.rma.norm)
gse <- data.frame(gse)

## Make sample names readable
colnames(gse) <- names
gse <- relocate(gse, grep("C", colnames(gse)), .before = colnames(gse)[4])
colnames(gse) <- gsub("brain, hippocampus: ", "", colnames(gse)) ## It's all hippocampus
colnames(gse) <- gsub("-.*", "", colnames(gse))
colnames(gse) <- gsub("h.*", "", colnames(gse))

## Run PCA
vargene <- apply(gse, 1, var)
names(vargene) <- rownames(gse)
vars <- sort(vargene,decreasing=T)
topgenes <- names(vars)[1:1000]
filtered <- gse[rownames(gse) %in% topgenes,]
PC <- prcomp(t(filtered))
varexpl <- (PC$sdev^2)/sum(PC$sdev^2)*100

## Plot PCA
dev.new(width=20, height=20)
colfunc <- colorRampPalette(c("red", "yellow"))
colors <- c(
  "#000080",  # Navy Blue
  "#1E90FF",  # Dodger Blue (lighter blue)
  "#87CEEB",  # Sky Blue
  "#e1c4ff",  # Lilac
  "#9B30FF",  # Purple
  "#FF3333",  # Red
  "#FF9933",  # Light orange
  "#FFFF00"   # Bright yellow
)
df <- data.frame(Colour = colors,
                 Time = c("1h", "4h", "6h", "8h",
                          "12h", "24h", "72h","240h"))
colours <- c(rep("green", 15), rep(df$Colour[c(1,3,6,7,8)], each = 3)) ## Green = Controls
plot(PC$x[,1],PC$x[,2],cex=1.2,bg=colours,
     xlab=paste0("PC1 (", round(varexpl[1], 2), "%)"),
     ylab=paste0("PC2 (", round(varexpl[2], 2), "%)"),
     main="GSE1834",las=1, pch=21
)

##----------------------------------------------------------------------------------------------

## GSE49030
gse <- getGEO("GSE49030")
names <- pData(gse[[1]])$title
cels = list.files("GSE49030/data/", pattern = "CEL") ## Data stored locally
raw.data <- ReadAffy(verbose = FALSE, 
                     filenames = paste("GSE49030/data", cels, sep = "/"),
                     cdfname = "mouse4302cdf") ## Load microarray data
data.rma.norm <- affy::rma(raw.data)
gse <- exprs(data.rma.norm)
gse <- data.frame(gse)

## Make sample names readable
colnames(gse) <- names
gse <- relocate(gse, grep("Control", colnames(gse)), .before = colnames(gse)[1]) ## Order columns to make downstream analysis simpler
colnames(gse) <- gsub("Control", "C", colnames(gse))
colnames(gse) <- gsub("Kainic acid ", "KA", colnames(gse))
colnames(gse) <- gsub(" other.*", ".2", colnames(gse))
gse <- relocate(gse, grep("KA1hr", colnames(gse)), .before = "KA4hrs 1") ## Order columns to make downstream analysis simpler
colnames(gse) <- gsub("r...", "", colnames(gse))
colnames(gse) <- gsub("r.*", "", colnames(gse))

## Run PCA
vargene <- apply(gse, 1, var)
names(vargene) <- rownames(gse)
vars <- sort(vargene,decreasing=T)
topgenes <- names(vars)[1:1000]
filtered <- gse[rownames(gse) %in% topgenes,]
PC <- prcomp(t(filtered))
varexpl <- (PC$sdev^2)/sum(PC$sdev^2)*100

## Plot PCA
dev.new(width=20, height=20)
colors <- c(
  "#000080",  # Navy Blue
  "#1E90FF",  # Dodger Blue (lighter blue)
  "#87CEEB",  # Sky Blue
  "#e1c4ff",  # Lilac
  "#9B30FF",  # Purple
  "#FF3333",  # Red
  "#FF9933",  # Light orange
  "#FFFF00"   # Bright yellow
)
df <- data.frame(Colour = colors,
                 Time = c("1h", "4h", "6h", "8h",
                          "12h", "24h", "72h","240h"))
colours <- c(rep("green", 9), 
             rep(df$Colour[c(1,2,4,6,2)], each = 3)) ## Green = Controls
plot(PC$x[,1],PC$x[,2],cex=1.2,bg=colours,
     xlab=paste0("PC1 (", round(varexpl[1], 2), "%)"),
     ylab=paste0("PC2 (", round(varexpl[2], 2), "%)"),
     pch=21,
     main="GSE49030",las=1
)

##-------------------------------------------------------------------------------

## GSE88992
gse <- getGEO("GSE88992")
names <- pData(gse[[1]])$title
cels = list.files("GSE88992/data/", pattern = "CEL") ## Data stored locally
raw.data <- ReadAffy(verbose = FALSE, 
                     filenames = paste("GSE88992/data", cels, sep = "/"),
                     cdfname = "mouse4302cdf") ## Load microarray data
data.rma.norm <- affy::rma(raw.data)
gse <- exprs(data.rma.norm)
gse <- data.frame(gse)

## Make sample names readable
colnames(gse) <- names
gse <- relocate(gse, grep("saline", colnames(gse)), .before = colnames(gse)[1])
colnames(gse) <- gsub("Hippocampus_ ", "", colnames(gse)) ## It's all hippocampus
colnames(gse) <- gsub("kainate_", "KA", colnames(gse))
colnames(gse) <- gsub("saline_", "C", colnames(gse)) ## Saline controls
colnames(gse) <- gsub("_rep.*", "", colnames(gse))

## Run PCA
vargene <- apply(gse, 1, var)
names(vargene) <- rownames(gse)
vars <- sort(vargene,decreasing=T)
topgenes <- names(vars)[1:1000]
filtered <- gse[rownames(gse) %in% topgenes,]
PC <- prcomp(t(filtered))
varexpl <- (PC$sdev^2)/sum(PC$sdev^2)*100

## Plot PCA
dev.new(width=20, height=20)
colors <- c(
  "#000080",  # Navy Blue
  "#1E90FF",  # Dodger Blue (lighter blue)
  "#87CEEB",  # Sky Blue
  "#e1c4ff",  # Lilac
  "#9B30FF",  # Purple
  "#FF3333",  # Red
  "#FF9933",  # Light orange
  "#FFFF00"  # Bright yellow
)
df <- data.frame(Colour = colors,
                 Time = c("1h", "4h", "6h", "8h",
                          "12h", "24h", "72h","240h"))
colours <- c(rep("green", 9), rep(c(df$Colour[c(3,5)]), each = 3), 
             rep(df$Colour[6], 2)) ## Green = Controls
plot(PC$x[,1],PC$x[,2],cex=1.2,bg=colours,
     xlab=paste0("PC1 (", round(varexpl[1], 2), "%)"),
     ylab=paste0("PC2 (", round(varexpl[2], 2), "%)"),
     main="GSE88992",las=1, pch=21
)

##------------------------------------------------------------------------------------

## GSE73878
gse <- getGEO("GSE73878")
names <- pData(gse[[1]])$title
gse <- exprs(gse[[1]])
gse <- data.frame(gse)

## Make sample names readable
colnames(gse) <- names
colnames(gse) <- gsub("L", ".", colnames(gse))
colnames(gse) <- gsub("_", "", colnames(gse))
colnames(gse) <- gsub("Sham", "C", colnames(gse))
colnames(gse) <- gsub("rep", "d", colnames(gse))
colnames(gse) <- gsub("7d", "07d", colnames(gse))
gse <- gse[,order(colnames(gse))] ## Order columns to make downstream analysis simpler
gse <- relocate(gse, grep("C.", colnames(gse)), .after = colnames(gse)[80])
gse <- relocate(gse, grep(".C", colnames(gse)), .before = colnames(gse)[1])

## Run PCA
vargene <- apply(gse, 1, var)
names(vargene) <- rownames(gse)
vars <- sort(vargene,decreasing=T)
topgenes <- names(vars)[1:1000]
filtered <- gse[rownames(gse) %in% topgenes,]
PC <- prcomp(t(filtered))
varexpl <- (PC$sdev^2)/sum(PC$sdev^2)*100

dev.new(width=20, height=20)
colours <- c(rep("blue", 28), 
             rep(c(rep("red", 8), rep("orange", 9), rep("yellow", 9)), 2))
plot(PC$x[,1],PC$x[,2],cex=1.2,bg=colours,
     xlab=paste0("PC1 (", round(varexpl[1], 2), "%)"),
     ylab=paste0("PC2 (", round(varexpl[2], 2), "%)"),
     main="Expression profiles beyond the latency period",
     las=1, pch=20
)
text(PC$x[,1], PC$x[,2], labels = colnames(gse), pos = 4) ## Check outlier sample names

## Remove outliers and redo PCA
gse3 <- gse[,!(colnames(gse) %in% c("C.C28d3", "C.KA60d3", "C.KA28d4", "I.C28d5"))]
vargene <- apply(gse3, 1, var)
names(vargene) <- rownames(gse3)
vars <- sort(vargene,decreasing=T)
topgenes <- names(vars)[1:1000]
filtered <- gse3[rownames(gse3) %in% topgenes,]
PC <- prcomp(t(filtered))
varexpl <- (PC$sdev^2)/sum(PC$sdev^2)*100

## Plot PCA
colors <- c(
  "#000080",  # Navy Blue
  "#00FFFF",  # Sky Blue
  "#e1c4ff",  # Lilac
  "#9B30FF",  # Purple
  "#FF3333",  # Red
  "#FF9933",  # Light orange
  "#FFFF00"   # Bright yellow
)
df <- data.frame(Colour = colors,
                 Time = c("1d", "3d", "7d",
                          "10d", "14d", "28d", "60d"))
colours <- c(rep("green", 13),
             rep("green", 13),
             rep(df$Colour[3], 8), 
             rep(df$Colour[6], 9), 
             rep(df$Colour[7], 9),
             rep(df$Colour[3], 8),
             rep(df$Colour[6], 8),
             rep(df$Colour[7], 8)) ## Green = Controls
shapes <- colnames(gse3)
shapes[grep("C\\.", shapes)] <- 21 ## Contralateral
shapes[grep("I\\.", shapes)] <- 22 ## Ipislateral
plot(PC$x[,1],PC$x[,2],cex=1.2,bg=colours,
     xlab=paste0("PC1 (", round(varexpl[1], 2), "%)"),
     ylab=paste0("PC2 (", round(varexpl[2], 2), "%)"),
     main="7 Days to 60 Days",
     las=1, pch=21
)

## Ipsilateral samples only
colours <- c(rep("green", 13),
             rep("green", 13),
             rep(df$Colour[3], 8), 
             rep(df$Colour[6], 9), 
             rep(df$Colour[7], 9),
             rep(df$Colour[3], 8),
             rep(df$Colour[6], 8),
             rep(df$Colour[7], 8)) ## Green = Controls
shapes <- colnames(gse3)
shapes[grep("C\\.", shapes)] <- 21 ## Contralateral
shapes[grep("I\\.", shapes)] <- 22 ## Ipsilateral
colours <- data.frame(cbind(colours, shapes))
colours$colours[colours$shapes=="21"] <- alpha("white", 0) ## Remove contralateral data points by increasing transparency to max
lines <- colours$colours
lines[lines!="#FFFFFF00"] <- "black" ## Only draw lines around non-transparent data points
par(mar = c(5, 4, 6, 2))
plot(PC$x[,1],PC$x[,2],cex=1.2,bg=colours$colours,
     xlab="", col=lines,
     ylab=paste0("PC2 (", round(varexpl[2], 2), "%)"),
     main="", 
     las=1, pch=22, xaxt="n"
) 
axis(3)

## Contralateral samples only
colours <- c(rep("green", 13),
             rep("green", 13),
             rep(df$Colour[3], 8), 
             rep(df$Colour[6], 9), 
             rep(df$Colour[7], 9),
             rep(df$Colour[3], 8),
             rep(df$Colour[6], 8),
             rep(df$Colour[7], 8)) ## Green = Controls
shapes <- colnames(gse3)
shapes[grep("C\\.", shapes)] <- 21 ## Contralateral
shapes[grep("I\\.", shapes)] <- 22 ## Ipsilateral
colours <- data.frame(cbind(colours, shapes))
colours$colours[colours$shapes=="22"] <- alpha("white", 0) ## Remove ipsilateral data points by increasing transparency to max
lines <- colours$colours
lines[lines!="#FFFFFF00"] <- "black" ## Only draw lines around non-transparent data points
par(mar = c(5, 4, 6, 2))
plot(PC$x[,1],PC$x[,2],cex=1.2,bg=colours$colours,
     xlab=paste0("PC1 (", round(varexpl[1], 2), "%)"), col=lines,
     ylab=paste0("PC2 (", round(varexpl[2], 2), "%)"),
     main="", 
     las=1, pch=25
) 
mtext(paste0("PC1 (", round(varexpl[1], 2), "%)"), side = 3, line = 3)
mtext("7 Days to 60 Days", side = 3, line = 4.3, font = 2, cex = 1.2)

##-------------------------------------------------------------------------------------------------------------

## GSE213393
gse <- read.csv("GSE213393/GSE213393_geneCounts.tsv", sep = "\t") ## Load stored data
rownames(gse) <- gse$Geneid
gse <- gse[,-1]

## Create a summarizedExperiment object from the count data
colnames(gse) <- gsub("Sample_", "", colnames(gse))
colData <- data.frame(location = rep(c("contralateral", "ipsilateral"), 24),
                      treatment = rep(c(rep("kainate", 4), rep("ctrl", 4)), 6),
                      time = rep(c("3d", "7d", "14d"), each = 16),
                      accession = paste0("GSM6585", 453:500),
                      animal = gsub("DC", "", gsub("DI", "", colnames(gse))),
                      sex = rep(c("male", "female"), 6, each = 2),
                      row.names = colnames(gse))
colData$location <- as.factor(colData$location)
colData$time <- as.factor(colData$time)
colData$animal <- as.factor(colData$animal)
colData$accession <- as.factor(colData$accession)
colData$treatment <- as.factor(colData$treatment)
colData$sex <- as.factor(colData$sex)
se <- SummarizedExperiment(as.matrix(gse), colData = colData)

## Filter out low counts
groupsize <- 3
mincount <- 10
keep <- rowSums(assay(se) >= mincount) >= groupsize
se <- se[keep,]

## Run PCA
counts <- rlog(assay(se)) ## Normalise
vargene <- apply(counts, 1, var)
names(vargene) <- rownames(counts)
vars <- sort(vargene,decreasing=T)
topgenes <- names(vars)[1:1000]
filtered <- counts[rownames(counts) %in% topgenes,]
PC <- prcomp(t(filtered))
varexpl <- (PC$sdev^2)/sum(PC$sdev^2)*100

dev.new()
colors <- c(
  "#000080",  # Navy Blue
  "#00FFFF",  # Sky Blue
  "#e1c4ff",  # Lilac
  "#9B30FF",  # Purple
  "#FF3333",  # Red
  "#FF9933",  # Light orange
  "#FFFF00"   # Bright yellow
)
df <- data.frame(Colour = colors,
                 Time = c("1d", "3d", "7d",
                          "10d", "14d", "28d", "60d"))
colours <- c(rep(c(df$Colour[2], "green"), 2, each = 4), 
             rep(c(df$Colour[3], "green"), 2, each = 4),
             rep(c(df$Colour[5], "green"), 2, each = 4)) ## Green = Controls
plot(PC$x[,1],PC$x[,2],cex=1.2,bg=colours,
     xlab=paste0("PC1 (", round(varexpl[1], 2), "%)"),
     ylab=paste0("PC2 (", round(varexpl[2], 2), "%)"),
     pch=21,main="3 Days to 14 Days",las=1
)

## Contralateral samples only
colours <- c(rep(c(rep(c(df$Colour[2], "#FFFFFF00"), 2), 
                   rep(c("green", "#FFFFFF00"), 2)), 2), 
             rep(c(rep(c(df$Colour[3], "#FFFFFF00"), 2), 
                   rep(c("green", "#FFFFFF00"), 2)), 2),
             rep(c(rep(c(df$Colour[5], "#FFFFFF00"), 2), 
                   rep(c("green", "#FFFFFF00"), 2)), 2)) ## Green = Controls
par(mar = c(5, 4, 6, 2))
plot(PC$x[,1],PC$x[,2],cex=1.2,bg=colours,
     xlab=paste0("PC1 (", round(varexpl[1], 2), "%)"),
     ylab=paste0("PC2 (", round(varexpl[2], 2), "%)"),
     pch=25,main="",las=1, col = c("black", "#FFFFFF00"),
)

## Ipsilateral samples only
colours <- c(rep(c(rep(c("#FFFFFF00", df$Colour[2]), 2), 
                   rep(c("#FFFFFF00", "green"), 2)), 2), 
             rep(c(rep(c("#FFFFFF00", df$Colour[3]), 2), 
                   rep(c("#FFFFFF00", "green"), 2)), 2),
             rep(c(rep(c("#FFFFFF00", df$Colour[5]), 2), 
                   rep(c("#FFFFFF00", "green"), 2)), 2)) ## Green = Controls
par(mar = c(5, 4, 6, 2))
plot(PC$x[,1],PC$x[,2],cex=1.2,bg=colours,
     xlab="", xaxt = "n",
     ylab=paste0("PC2 (", round(varexpl[2], 2), "%)"),
     pch=22,main="",las=1, col = c("#FFFFFF00", "black")
)
axis(3)
mtext(paste0("PC1 (", round(varexpl[1], 2), "%)"), side = 3, line = 3)
mtext("3 Days to 14 Days", side = 3, line = 4.3, font = 2, cex = 1.2)

##-------------------------------------------------------------------------------------------

## GSE47752
gse <- getGEO("GSE47752")
names <- pData(gse[[1]])$title
cels = list.files("GSE47752/data/", pattern = "CEL") ## Data stored locally
raw.data <- ReadAffy(verbose = FALSE, 
                     filenames = paste("GSE47752/data", cels, sep = "/"),
                     cdfname = "rat2302cdf") ## Load microarray data
data.rma.norm <- affy::rma(raw.data)
gse <- exprs(data.rma.norm)
gse <- data.frame(gse)

## Make labels shorter and readable
colnames(gse) <- names
gse <- gse[,grep("kainate", colnames(gse))]
foo <- function(x){
  strsplit(x, "-")
}
a <- lapply(colnames(gse), foo)
for(i in 1:length(a)) {
  colnames(gse)[i] <- paste0(a[[i]][[1]][[3]], 
                             a[[i]][[1]][[4]],
                             a[[i]][[1]][[5]])
}
colnames(gse) <- gsub("Nadler", "N.", colnames(gse))
colnames(gse) <- gsub("Wadman", "W.", colnames(gse))
colnames(gse) <- gsub("day ", "", colnames(gse))
colnames(gse) <- gsub("control", "C", colnames(gse))
colnames(gse) <- gsub("rat ", ".", colnames(gse))
colnames(gse) <- gsub("day", "", colnames(gse))
gse <- gse[,-7]
gse$N.3.5 <- apply(gse[,grep("N.3.5", colnames(gse))], 1, mean) ## Multiple recordings from the same sample
gse <- gse[,order(colnames(gse))] ## Order columns for easier downstream analysis
gse <- gse[,-18]

## Run PCA
vargene <- apply(gse, 1, var)
names(vargene) <- rownames(gse)
vars <- sort(vargene,decreasing=T)
topgenes <- names(vars)[1:1000]
filtered <- gse[rownames(gse) %in% topgenes,]
PC <- prcomp(t(filtered))
varexpl <- (PC$sdev^2)/sum(PC$sdev^2)*100

## Plot PCA
dev.new(width=20, height=20)
colors <- c(
  "#000080",  # Navy Blue
  "#00FFFF",  # Sky Blue
  "#e1c4ff",  # Lilac
  "#9B30FF",  # Purple
  "#FF3333",  # Red
  "#FF9933",  # Light orange
  "#FFFF00"   # Bright yellow
)
df <- data.frame(Colour = colors,
                 Time = c("1d", "3d", "7d",
                          "10d", "14d", "28d", "60d"))
colours <- rep(c(df$Colour[c(1,2,4)],"green"), 2, each = 6) ## Green = Controls
plot(PC$x[,2],PC$x[,1],cex=1.2,bg=colours,
     xlab=paste0("PC2 (", round(varexpl[2], 2), "%)"),
     ylab=paste0("PC1 (", round(varexpl[1], 2), "%)"),
     main="1 Day to 10 Days",las=1, pch=rep(21,each = 48)
)

##-----------------------------------------------------------------------------------------------

## Plot legends
plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1) ## Empty plot
colors <- c(
  "#000080",  # Navy Blue
  "#1E90FF",  # Dodger Blue (lighter blue)
  "#87CEEB",  # Sky Blue
  "#e1c4ff",  # Lilac
  "#9B30FF",  # Purple
  "#FF3333",  # Red
  "#FF9933",  # Light orange
  "#FFFF00"  # Bright yellow
)
legend("topleft", pt.bg=c("green",unique(colors[1:6])), legend = c("Control",
                                      "1h","4h","6h", "8h",
                                      "12h", "24h"),
       
       pch = 21, cex=1, col = "black")
legend("top", pt.bg=c("green",unique(colors[c(1,3,6:8)])), legend = c("Control",
                                                                   "1h","6h","24h",
                                                                   "3d", "10d"),
       
       pch = 21, cex=1, col = "black")

colors <- c(
  "#000080",  # Navy Blue
  "#00FFFF",  # Sky Blue
  "#e1c4ff",  # Lilac
  "#9B30FF",  # Purple
  "#FF3333",  # Red
  "#FF9933",  # Light orange
  "#FFFF00"  # Bright yellow
)
legend("bottomleft", pt.bg=c("green",unique(colors)), 
       legend = c("Control", 
                  "1d", "3d", "7d",
                  "10d", "14d", "28d", "60d"),
       
       pch = 21, cex=1, col = "black")
legend("bottom", pt.bg="white",
       legend = c("Ipsilateral", "Contralateral"),
       pch = c(22,25), cex=1)
