## Load required packages
library(affy)
library(dplyr)
library(GEOquery)
library(DESeq2)
library(SummarizedExperiment)
library(limma)
library(rgu34acdf)
library(rgu34a.db)
library(mouse4302.db)
library(mouse4302cdf)
library(rat2302cdf)
library(rat2302.db)
library(org.Mm.eg.db)
library(AnnotationDbi)

## Data was downloaded using getGEOSuppFiles()
## First, run contrasts to determine differential expression
## Then, gene expression can be plotted using plot_expression.R
##-------------------------------------------------------------------------------------

## GSE1831
gse <- getGEO("GSE1831")
names <- pData(gse[[1]])$title
cels = list.files("GSE1831/data/", pattern = "cel") ## Data stored locally
raw.data <- ReadAffy(verbose = FALSE, 
                     filenames = paste("GSE1831/data", cels, sep = "/"),
                     cdfname = "rgu34acdf") ## Load microarray data
data.rma.norm <- affy::rma(raw.data)
data <- exprs(data.rma.norm)

## Filter low intensity probes
medians <- rowMedians(data)
hist_res <- hist(medians, 100, col = "cornsilk1", freq = FALSE, 
                 main = "Histogram of the median intensities", 
                 border = "antiquewhite4",
                 xlab = "Median intensities")
filter <- apply(data, 1, function(x){sum(x > 6) >= 3}) ## Visual cutoff based on hist_res
data <- subset(data.rma.norm, filter)

## Filter probes that match to multiple genes
match <- AnnotationDbi::select(rgu34a.db,
                               keys = featureNames(data),
                               columns = c("SYMBOL", "ENSEMBL"),
                               keytype = "PROBEID")
match <- subset(match, !is.na(ENSEMBL))
match_grouped <- group_by(match, PROBEID)
match_summarized <- dplyr::summarize(match_grouped, no_of_matches = n_distinct(SYMBOL))
match_filtered <- filter(match_summarized, no_of_matches > 1)
ids_to_exclude <- (featureNames(data) %in% match_filtered$PROBEID)
data <- subset(data, !ids_to_exclude)
match_summarized <- dplyr::summarize(match_grouped, no_of_matches = n_distinct(ENSEMBL))
match_filtered <- filter(match_summarized, no_of_matches > 1)
ids_to_exclude <- (featureNames(data) %in% match_filtered$PROBEID)
data <- subset(data, !ids_to_exclude)
fData(data)$PROBEID <- rownames(fData(data))
fData(data) <- left_join(fData(data), match)
rownames(fData(data)) <- fData(data)$PROBEID
fData(data) <- fData(data)[,-1]

## Prepare names of factors
kainate <- names
kainate <- gsub(".*KA", "KA", kainate)
kainate <- gsub("KA.*", "KA", kainate)
kainate <- gsub(".*Ctr", "Ctr", kainate)
kainate <- gsub("Ctr.*", "Ctr", kainate)

time <- names
time <- gsub(".*KA.", "", time)
time <- gsub(".*Ctr.", "", time)
time <- gsub("-.*", "", time)

## Run contrasts
factor <- factor(paste(kainate,time,sep="."))
design <- model.matrix(~ 0 + factor)
colnames(design) <- levels(factor)
fit <- lmFit(data, design)
cont.matrix <- makeContrasts(`1h` = KA.1h-Ctr.1h,
                             `6h` = KA.6h-Ctr.6h,
                             `24h` = KA.24h-Ctr.24h,
                             `72h` = KA.72h-Ctr.72h,
                             `240h` = KA.240h-Ctr.240h,
                             levels = design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)
`1h` <- topTable(fit2, p.value = 0.05, coef = 1, number = Inf)
`6h` <- topTable(fit2, p.value = 0.05, coef = 2, number = Inf)
`24h` <- topTable(fit2, p.value = 0.05, coef = 3, number = Inf)
`72h` <- topTable(fit2, p.value = 0.05, coef = 4, number = Inf)
`240h` <- topTable(fit2, p.value = 0.05, coef = 5, number = Inf)

## Write DE results to file
write.csv(`1h`, "GSE1831/1h.csv")
write.csv(`6h`, "GSE1831/6h.csv")
write.csv(`24h`, "GSE1831/24h.csv")
write.csv(`72h`, "GSE1831/72h.csv")
write.csv(`240h`, "GSE1831/240h.csv")

##-------------------------------------------------------------------------

## GSE1833
gse <- getGEO("GSE1833")
names <- pData(gse[[1]])$title
cels = list.files("GSE1833/data/", pattern = "cel") ## Data stored locally
raw.data <- ReadAffy(verbose = FALSE, 
                     filenames = paste("GSE1833/data", cels, sep = "/"),
                     cdfname = "rae230acdf") ## Load microarray data
data.rma.norm <- affy::rma(raw.data)
data <- exprs(data.rma.norm)

## Filter low intensity probes
medians <- rowMedians(data)
hist_res <- hist(medians, 100, col = "cornsilk1", freq = FALSE, 
                 main = "Histogram of the median intensities", 
                 border = "antiquewhite4",
                 xlab = "Median intensities")
filter <- apply(data, 1, function(x){sum(x > 6) >= 3}) ## Visual cutoff based on hist_res
data <- subset(data.rma.norm, filter)

## Filter probes that match multiple genes
match <- AnnotationDbi::select(rae230a.db,
                               keys = featureNames(data),
                               columns = c("SYMBOL", "ENSEMBL"),
                               keytype = "PROBEID")
match <- subset(match, !is.na(ENSEMBL))
match_grouped <- group_by(match, PROBEID)
match_summarized <- dplyr::summarize(match_grouped, no_of_matches = n_distinct(SYMBOL))
match_filtered <- filter(match_summarized, no_of_matches > 1)
ids_to_exclude <- (featureNames(data) %in% match_filtered$PROBEID)
data <- subset(data, !ids_to_exclude)
match_summarized <- dplyr::summarize(match_grouped, no_of_matches = n_distinct(ENSEMBL))
match_filtered <- filter(match_summarized, no_of_matches > 1)
ids_to_exclude <- (featureNames(data) %in% match_filtered$PROBEID)
data <- subset(data, !ids_to_exclude)
fData(data)$PROBEID <- rownames(fData(data))
fData(data) <- left_join(fData(data), match)
rownames(fData(data)) <- fData(data)$PROBEID
fData(data) <- fData(data)[,-1]

## Prepare names of factors
kainate <- names
kainate <- gsub(".*KA", "KA", kainate)
kainate <- gsub("KA.*", "KA", kainate)
kainate <- gsub(".*C", "C", kainate)
kainate <- gsub("C.*", "C", kainate)

enrichment <- names
enrichment <- gsub(".*EN.*", "EN", enrichment)
enrichment[-grep("EN", enrichment)] <- "NEN"

## Run contrasts
factor <- factor(paste(kainate,enrichment,sep="."))
design <- model.matrix(~ 0 + factor)
colnames(design) <- levels(factor)
fit <- lmFit(data, design)
cont.matrix <- makeContrasts(KA.NEN-C.NEN,
                             KA.EN-C.EN,
                             levels = design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)
NEN <- topTable(fit2, p.value = 0.05, coef = 1, number = Inf)
EN <- topTable(fit2, p.value = 0.05, coef = 2, number = Inf)

## Write DE results to file
write.csv(NEN, "GSE1833/NEN.csv")
write.csv(EN, "GSE1833/Enriched.csv")

##------------------------------------------------------------------------------------

## GSE1834
gse <- getGEO("GSE1834")
names <- pData(gse[[1]])$title
cels = list.files("GSE1834/data/", pattern = "cel") ## Data stored locally
raw.data <- ReadAffy(verbose = FALSE, 
                     filenames = paste("GSE1834/data", cels, sep = "/"),
                     cdfname = "rgu34acdf") ## Load microarray data
data.rma.norm <- affy::rma(raw.data)
data <- exprs(data.rma.norm)

## Filter low intensity probes
medians <- rowMedians(data)
hist_res <- hist(medians, 100, col = "cornsilk1", freq = FALSE, 
                 main = "Histogram of the median intensities", 
                 border = "antiquewhite4",
                 xlab = "Median intensities")
filter <- apply(data, 1, function(x){sum(x > 6) >= 3}) ## Visual cutoff based on hist_res
data <- subset(data.rma.norm, filter)

## Filter probes that match multiple genes
match <- AnnotationDbi::select(rgu34a.db,
                               keys = featureNames(data),
                               columns = c("SYMBOL", "ENSEMBL"),
                               keytype = "PROBEID")
match <- subset(match, !is.na(ENSEMBL))
match_grouped <- group_by(match, PROBEID)
match_summarized <- dplyr::summarize(match_grouped, no_of_matches = n_distinct(SYMBOL))
match_filtered <- filter(match_summarized, no_of_matches > 1)
ids_to_exclude <- (featureNames(data) %in% match_filtered$PROBEID)
data <- subset(data, !ids_to_exclude)
match_summarized <- dplyr::summarize(match_grouped, no_of_matches = n_distinct(ENSEMBL))
match_filtered <- filter(match_summarized, no_of_matches > 1)
ids_to_exclude <- (featureNames(data) %in% match_filtered$PROBEID)
data <- subset(data, !ids_to_exclude)
fData(data)$PROBEID <- rownames(fData(data))
fData(data) <- left_join(fData(data), match)
rownames(fData(data)) <- fData(data)$PROBEID
fData(data) <- fData(data)[,-1]

## Prepare names of factors
kainate <- names
kainate <- gsub(".*KA", "KA", kainate)
kainate <- gsub("KA.*", "KA", kainate)
kainate <- gsub(".*C", "C", kainate)
kainate <- gsub("C.*", "C", kainate)

time <- names
time <- gsub(".*KA", "", time)
time <- gsub(".*C", "", time)
time <- gsub("-.*", "", time)

## Run contrasts
factor <- factor(paste(kainate,time,sep="."))
design <- model.matrix(~ 0 + factor)
colnames(design) <- levels(factor)
fit <- lmFit(data, design)
cont.matrix <- makeContrasts(`1h` = KA.1-C.1,
                             `6h` = KA.6-C.6,
                             `24h` = KA.24-C.24,
                             `72h` = KA.72-C.72,
                             `240h` = KA.240-C.240,
                             levels = design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)
`1h` <- topTable(fit2, p.value = 0.05, coef = 1, number = Inf)
`6h` <- topTable(fit2, p.value = 0.05, coef = 2, number = Inf)
`24h` <- topTable(fit2, p.value = 0.05, coef = 3, number = Inf)
`72h` <- topTable(fit2, p.value = 0.05, coef = 4, number = Inf)
`240h` <- topTable(fit2, p.value = 0.05, coef = 5, number = Inf)

## Write DE results to file
write.csv(`1h`, "GSE1834/1h.csv")
write.csv(`6h`, "GSE1834/6h.csv")
write.csv(`24h`, "GSE1834/24h.csv")
write.csv(`72h`, "GSE1834/72h.csv")
write.csv(`240h`, "GSE1834/240h.csv")

##--------------------------------------------------------------------------------

## GSE10923
gse <- getGEO("GSE10923")
names <- pData(gse[[1]])$title
cels = list.files("GSE10923/data/", pattern = "CEL") ## Data stored locally
raw.data <- ReadAffy(verbose = FALSE, 
                     filenames = paste("GSE10923/data", cels, sep = "/"),
                     cdfname = "rat2302cdf") ## Load microarray data
data.rma.norm <- affy::rma(raw.data)
data <- exprs(data.rma.norm)

## Filtre low intensity probes
medians <- rowMedians(data)
hist_res <- hist(medians, 100, col = "cornsilk1", freq = FALSE, 
                 main = "Histogram of the median intensities", 
                 border = "antiquewhite4",
                 xlab = "Median intensities")
filter <- apply(data, 1, function(x){sum(x > 3) >= 3})
data <- subset(data.rma.norm, filter)

## Filter probes that match multiple genes
match <- AnnotationDbi::select(rat2302.db,
                               keys = featureNames(data),
                               columns = c("SYMBOL", "ENSEMBL"),
                               keytype = "PROBEID")
match <- subset(match, !is.na(ENSEMBL))
match_grouped <- group_by(match, PROBEID)
match_summarized <- dplyr::summarize(match_grouped, no_of_matches = n_distinct(SYMBOL))
match_filtered <- filter(match_summarized, no_of_matches > 1)
ids_to_exclude <- (featureNames(data) %in% match_filtered$PROBEID)
data <- subset(data, !ids_to_exclude)
match_summarized <- dplyr::summarize(match_grouped, no_of_matches = n_distinct(ENSEMBL))
match_filtered <- filter(match_summarized, no_of_matches > 1)
ids_to_exclude <- (featureNames(data) %in% match_filtered$PROBEID)
data <- subset(data, !ids_to_exclude)
fData(data)$PROBEID <- rownames(fData(data))
fData(data) <- left_join(fData(data), match)
rownames(fData(data)) <- fData(data)$PROBEID
fData(data) <- fData(data)[,-1]

## Prepare names of factors
kainate <- names
kainate <- gsub("Kainate.*", "KA", kainate)
kainate[-grep("KA", kainate)] <- "C"

peptide <- names
peptide <- gsub(".*Neuroprotective peptide.*", "NP", peptide, ignore.case = T)
peptide[-grep("NP", peptide)] <- "ALONE"

## Run contrasts
factor <- factor(paste(kainate,peptide,sep="."))
design <- model.matrix(~ 0 + kainate)
colnames(design) <- levels(factor)
fit <- lmFit(data, design)
cont.matrix <- makeContrasts(kainateKA-kainateC,
                             levels = design) ## We only care about the results without the neuroprotective peptide
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)
kainate <- topTable(fit2, p.value = 0.05, coef = 1, number = Inf)

## Write DE results to file
write.csv(kainate, "GSE10923/differential expression.csv")

##-----------------------------------------------------------------------------------------

## GSE47752
cels = list.files("GSE47752/data/", pattern = "CEL") ## Data stored locally
raw.data <- ReadAffy(verbose = FALSE, 
                     filenames = paste("GSE47752/data", cels, sep = "/"),
                     cdfname = "rat2302cdf") ## Load microarray data
filter <- colnames(raw.data)[grepl("WW", rownames(raw.data@phenoData@data)) |
                             grepl("VN", rownames(raw.data@phenoData@data))] ## Only keep kainate-related data
raw.filter <- raw.data[,filter]
filter <- colnames(raw.filter)[grepl("rehyb", rownames(raw.filter@phenoData@data))]
raw.filter <- raw.filter[,rownames(raw.filter@phenoData@data)!=filter]
data.rma.norm <- affy::rma(raw.filter)
names <- rownames(pData(data.rma.norm))
gse <- exprs(data.rma.norm)
gse <- data.frame(gse)

## Filter low intensity probes
medians <- rowMedians(data)
hist_res <- hist(medians, 100, col = "cornsilk1", freq = FALSE, 
                 main = "Histogram of the median intensities", 
                 border = "antiquewhite4",
                 xlab = "Median intensities")
filter <- apply(data, 1, function(x){sum(x > 4.2) >= 3}) ## Visual cutoff based on hist_res
data <- subset(data.rma.norm, filter)

## Filter probes that match to multiple genes
match <- AnnotationDbi::select(rat2302.db,
                               keys = featureNames(data),
                               columns = c("SYMBOL", "ENSEMBL"),
                               keytype = "PROBEID")
match <- subset(match, !is.na(ENSEMBL))
match_grouped <- group_by(match, PROBEID)
match_summarized <- dplyr::summarize(match_grouped, no_of_matches = n_distinct(SYMBOL))
match_filtered <- filter(match_summarized, no_of_matches > 1)
ids_to_exclude <- (featureNames(data) %in% match_filtered$PROBEID)
data <- subset(data, !ids_to_exclude)
match_summarized <- dplyr::summarize(match_grouped, no_of_matches = n_distinct(ENSEMBL))
match_filtered <- filter(match_summarized, no_of_matches > 1)
ids_to_exclude <- (featureNames(data) %in% match_filtered$PROBEID)
data <- subset(data, !ids_to_exclude)
fData(data)$PROBEID <- rownames(fData(data))
fData(data) <- left_join(fData(data), match)
rownames(fData(data)) <- fData(data)$PROBEID
fData(data) <- fData(data)[,-1]

## Prepare names of factors
kainate <- names
kainate <- gsub(".*KA.*", "KA", kainate)
kainate <- gsub(".*C.*", "C", kainate)

time <- names
foo <- function(x){
  strsplit(x, "-")
}
a <- lapply(time, foo)
for(i in 1:length(a)) {
  time[i] <- paste0(a[[i]][[1]][[4]],
                    a[[i]][[1]][[5]],
                    a[[i]][[1]][[6]])
}
time[-grep("Day", time, ignore.case = T)] <- "0d" ## Controls don't have a day label, so we assign them "0d"
time <- gsub(".*WW", "", time)
time <- gsub(".*VN", "", time)
time <- gsub("Day.*", "d", time, ignore.case = T)
time <- gsub("_", "", time)

strain <- names
strain <- gsub(".*WW.*", "Whadman", strain)
strain <- gsub(".*VN.*", "Nadler", strain)

## Run contrasts
factor <- factor(paste(kainate,time,strain,sep="."))
design <- model.matrix(~ 0 + factor)
colnames(design) <- levels(factor)
fit <- lmFit(data, design)
cont.matrix <- makeContrasts(VN1d = KA.1d.Nadler-C.0d.Nadler,
                             VN3d = KA.3d.Nadler-C.0d.Nadler,
                             VN10d = KA.10d.Nadler-C.0d.Nadler,
                             WW1d = KA.1d.Whadman-C.0d.Whadman,
                             WW3d = KA.3d.Whadman-C.0d.Whadman,
                             WW10d = KA.10d.Whadman-C.0d.Whadman,
                             levels = design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)
VN1d <- topTable(fit2, p.value = 0.05, coef = 1, number = Inf)
VN3d <- topTable(fit2, p.value = 0.05, coef = 2, number = Inf)
VN10d <- topTable(fit2, p.value = 0.05, coef = 3, number = Inf)
WW1d <- topTable(fit2, p.value = 0.05, coef = 4, number = Inf)
WW3d <- topTable(fit2, p.value = 0.05, coef = 5, number = Inf)
WW10d <- topTable(fit2, p.value = 0.05, coef = 6, number = Inf)

## Write DE results to file
write.csv(VN1d, "GSE47752/Nadler 1 Day.csv")
write.csv(VN3d, "GSE47752/Nadler 3 Days.csv")
write.csv(VN10d, "GSE47752/Nadler 10 Days.csv")
write.csv(WW1d, "GSE47752/Whadman 1 Day.csv")
write.csv(WW3d, "GSE47752/Whadman 3 Days.csv")
write.csv(WW10d, "GSE47752/Whadman 10 Days.csv")

##---------------------------------------------------------------------------------------------

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

## Filter low intensity probes
medians <- rowMedians(data)
hist_res <- hist(medians, 100, col = "cornsilk1", freq = FALSE, 
                 main = "Histogram of the median intensities", 
                 border = "antiquewhite4",
                 xlab = "Median intensities")
filter <- apply(data, 1, function(x){sum(x > 2.5) >= 3})
data <- subset(data.rma.norm, filter)

## Filter probes that match multiple genes
match <- AnnotationDbi::select(mouse4302.db,
                               keys = featureNames(data),
                               columns = c("SYMBOL", "ENSEMBL"),
                               keytype = "PROBEID")
match <- subset(match, !is.na(ENSEMBL))
match_grouped <- group_by(match, PROBEID)
match_summarized <- dplyr::summarize(match_grouped, no_of_matches = n_distinct(SYMBOL))
match_filtered <- filter(match_summarized, no_of_matches > 1)
ids_to_exclude <- (featureNames(data) %in% match_filtered$PROBEID)
data <- subset(data, !ids_to_exclude)
match_summarized <- dplyr::summarize(match_grouped, no_of_matches = n_distinct(ENSEMBL))
match_filtered <- filter(match_summarized, no_of_matches > 1)
ids_to_exclude <- (featureNames(data) %in% match_filtered$PROBEID)
data <- subset(data, !ids_to_exclude)
fData(data)$PROBEID <- rownames(fData(data))
fData(data) <- left_join(fData(data), match)
rownames(fData(data)) <- fData(data)$PROBEID
fData(data) <- fData(data)[,-1]

## Prepare names of factors
kainate <- names
kainate <- gsub(".*Kainic acid.*", "KA", kainate)
kainate <- gsub(".*Control.*", "C", kainate)


time <- names
time <- gsub("Control.*", "0h 0h 0h", time) ## Controls don't have a time lable, so "0h" is assigned
foo <- function(x){
  strsplit(x, " ")
}
a <- lapply(time, foo)
for(i in 1:length(a)) {
  time[i] <- paste0(a[[i]][[1]][[3]])
}
time <- gsub("hr.*", "h", time)

strain <- names
strain <- gsub(".*other.*", "2", strain)
strain[strain!="2"] <- "1"

## Run contrasts
factor <- factor(paste(kainate,time,strain,sep="."))
design <- model.matrix(~ 0 + factor)
colnames(design) <- levels(factor)
fit <- lmFit(data, design)
cont.matrix <- makeContrasts(`1h` = KA.1h.1-C.0h.1,
                             `4h` = KA.4h.1-C.0h.1,
                             `8h` = KA.8h.1-C.0h.1,
                             `24h` = KA.24h.1-C.0h.1,
                             `4h Alt` = KA.4h.2-C.0h.2,
                             levels = design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)
`1h` <- topTable(fit2, p.value = 0.05, coef = 1, number = Inf)
`4h` <- topTable(fit2, p.value = 0.05, coef = 2, number = Inf)
`8h` <- topTable(fit2, p.value = 0.05, coef = 3, number = Inf)
`24h` <- topTable(fit2, p.value = 0.05, coef = 4, number = Inf)
`4h Alt` <- topTable(fit2, p.value = 0.05, coef = 5, number = Inf)

## Write DE results to file
write.csv(`1h`, "GSE49030/1h.csv")
write.csv(`4h`, "GSE49030/4h.csv")
write.csv(`8h`, "GSE49030/8h.csv")
write.csv(`24h`, "GSE49030/24h.csv")
write.csv(`4h Alt`, "GSE49030/4h Alt.csv")

##-------------------------------------------------------------------------------------------

## GSE73878
gse <- getGEO("GSE73878")
names <- pData(gse[[1]])$title
data <- exprs(gse[[1]])

## Filter low intensity probes
medians <- rowMedians(data)
hist_res <- hist(medians, 100, col = "cornsilk1", freq = FALSE, 
                 main = "Histogram of the median intensities", 
                 border = "antiquewhite4",
                 xlab = "Median intensities")
filter <- apply(data, 1, function(x){sum(x > 4) >= 3})
data <- subset(gse[[1]], filter)

## Filter probes that match multiple genes
match <- read.csv("GPL6885-11608 Illumina Microbead Array.txt", sep = "\t") ## Downloaded from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL6885
match <- match[,c(1,6)]
colnames(match) <- c("PROBEID", "SYMBOL")
match_grouped <- group_by(match, PROBEID)
match_summarized <- dplyr::summarize(match_grouped, no_of_matches = n_distinct(SYMBOL))
match_filtered <- filter(match_summarized, no_of_matches > 1)
ids_to_exclude <- (featureNames(data) %in% match_filtered$PROBEID)
data <- subset(data, !ids_to_exclude)
fData(data)$PROBEID <- rownames(fData(data))
fData(data) <- left_join(fData(data), match)
rownames(fData(data)) <- fData(data)$PROBEID
fData(data) <- fData(data)[,c(8,12)]

## Prepare names of factors
kainate <- names
foo <- function(x){
  strsplit(x, "_")
}
a <- lapply(kainate, foo)
for(i in 1:length(a)) {
  kainate[i] <- paste0(a[[i]][[1]][[2]])
}
kainate <- gsub("Sham", "CON", kainate)


time <- names
foo <- function(x){
  strsplit(x, "_")
}
a <- lapply(time, foo)
for(i in 1:length(a)) {
  time[i] <- paste0(a[[i]][[1]][[3]])
}
time <- paste0(time, "d")

location <- names
foo <- function(x){
  strsplit(x, "_")
}
a <- lapply(location, foo)
for(i in 1:length(a)) {
  location[i] <- paste0(a[[i]][[1]][[1]])
}

## Run contrasts
factor <- factor(paste(kainate,time,location,sep="."))
design <- model.matrix(~ 0 + factor)
colnames(design) <- levels(factor)
fit <- lmFit(data, design)
cont.matrix <- makeContrasts(IL.7d = KA.7d.IL-CON.7d.IL,
                             IL.28d = KA.28d.IL-CON.28d.IL,
                             IL.60d = KA.60d.IL-CON.60d.IL,
                             CL.7d = KA.7d.CL-CON.7d.CL,
                             CL.28d = KA.28d.CL-CON.28d.CL,
                             CL.60d = KA.60d.CL-CON.60d.CL,
                             levels = design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)
IL.7d <- topTable(fit2, p.value = 0.05, coef = 1, number = Inf)
IL.28d <- topTable(fit2, p.value = 0.05, coef = 2, number = Inf)
IL.60d <- topTable(fit2, p.value = 0.05, coef = 3, number = Inf)
CL.7d <- topTable(fit2, p.value = 0.05, coef = 4, number = Inf)
CL.28d <- topTable(fit2, p.value = 0.05, coef = 5, number = Inf)
CL.60d <- topTable(fit2, p.value = 0.05, coef = 6, number = Inf)

## Write DE results to file
write.csv(IL.7d, "GSE73878/7d IL.csv")
write.csv(IL.28d, "GSE73878/28d IL.csv")
write.csv(IL.60d, "GSE73878/60d IL.csv")
write.csv(CL.7d, "GSE73878/7d CL.csv")
write.csv(CL.28d, "GSE73878/28d CL.csv")
write.csv(CL.60d, "GSE73878/60d CL.csv")

##---------------------------------------------------------------------------------------

## GSE79129
gse <- getGEO("GSE79129")
names <- pData(gse[[1]])$title
data <- exprs(gse[[1]])

## Filter low intensity probes
medians <- rowMedians(data)
hist_res <- hist(medians, 100, col = "cornsilk1", freq = FALSE, 
                 main = "Histogram of the median intensities", 
                 border = "antiquewhite4",
                 xlab = "Median intensities")
filter <- apply(data, 1, function(x){sum(x > 14) >= 3})
data <- subset(gse[[1]], filter)
fData(data) <- fData(data)[,c(8,12)]

## Prepare names of factors
kainate <- rep(c("CON", "KA"), each = 3)

## Run contrasts
factor <- factor(paste(kainate,sep="."))
design <- model.matrix(~ 0 + factor)
colnames(design) <- levels(factor)
fit <- lmFit(data, design)
cont.matrix <- makeContrasts(KA = KA-CON,
                             levels = design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)
KA <- topTable(fit2, p.value = 0.05, coef = 1, number = Inf)

## Write DE results to file
write.csv(KA, "GSE79129/72h Differential expression.csv")

##----------------------------------------------------------------------------------

## GSE88992
gse <- getGEO("GSE88992")
names <- pData(gse[[1]])$title
cels = list.files("GSE88992/data/", pattern = "CEL") ## Data stored locally
raw.data <- ReadAffy(verbose = FALSE, 
                     filenames = paste("GSE88992/data", cels, sep = "/"),
                     cdfname = "mouse4302cdf") ## Load microarray data
data.rma.norm <- affy::rma(raw.data)
data <- exprs(data.rma.norm)

## Filter low intensity probes
medians <- rowMedians(data)
hist_res <- hist(medians, 100, col = "cornsilk1", freq = FALSE, 
                 main = "Histogram of the median intensities", 
                 border = "antiquewhite4",
                 xlab = "Median intensities")
filter <- apply(data, 1, function(x){sum(x > 2.5) >= 3})
data <- subset(data.rma.norm, filter)

## Filter probes that match to multiple genes
match <- AnnotationDbi::select(mouse4302.db,
                               keys = featureNames(data),
                               columns = c("SYMBOL", "ENSEMBL"),
                               keytype = "PROBEID")
match <- subset(match, !is.na(ENSEMBL))
match_grouped <- group_by(match, PROBEID)
match_summarized <- dplyr::summarize(match_grouped, no_of_matches = n_distinct(SYMBOL))
match_filtered <- filter(match_summarized, no_of_matches > 1)
ids_to_exclude <- (featureNames(data) %in% match_filtered$PROBEID)
data <- subset(data, !ids_to_exclude)
match_summarized <- dplyr::summarize(match_grouped, no_of_matches = n_distinct(ENSEMBL))
match_filtered <- filter(match_summarized, no_of_matches > 1)
ids_to_exclude <- (featureNames(data) %in% match_filtered$PROBEID)
data <- subset(data, !ids_to_exclude)
fData(data)$PROBEID <- rownames(fData(data))
fData(data) <- left_join(fData(data), match)
rownames(fData(data)) <- fData(data)$PROBEID
fData(data) <- fData(data)[,-1]

## Prepare names of factors
kainate <- names
kainate <- gsub(".*kainate.*", "KA", kainate)
kainate <- gsub(".*saline.*", "C", kainate)


time <- names
foo <- function(x){
  strsplit(x, "_")
}
a <- lapply(time, foo)
for(i in 1:length(a)) {
  time[i] <- paste0(a[[i]][[1]][[3]])
}

## Run contrasts
factor <- factor(paste(kainate,time,sep="."))
design <- model.matrix(~ 0 + factor)
colnames(design) <- levels(factor)
fit <- lmFit(data, design)
cont.matrix <- makeContrasts(`8h` = KA.6h-C.6h,
                             `12h` = KA.12h-C.12h,
                             `24h` = KA.24h-C.24h,
                             levels = design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)
`6h` <- topTable(fit2, p.value = 0.05, coef = 1, number = Inf)
`12h` <- topTable(fit2, p.value = 0.05, coef = 2, number = Inf)
`24h` <- topTable(fit2, p.value = 0.05, coef = 3, number = Inf)

## Write contrasts to file
write.csv(`6h`, "GSE88992/6h.csv")
write.csv(`12h`, "GSE88992/12h.csv")
write.csv(`24h`, "GSE88992/24h.csv")

##------------------------------------------------------------------------------------

## GSE122228
gse <- getGEO("GSE122228")
names <- pData(gse[[1]])$title
cels = list.files("GSE122228/data/", pattern = "cel") ## Data stored locally
raw.data <- ReadAffy(verbose = FALSE, 
                     filenames = paste("GSE122228/data", cels, sep = "/"),
                     cdfname = "htmg430pmcdf") ## Load microarray data
data.rma.norm <- affy::rma(raw.data)
data <- exprs(data.rma.norm)

## Filter low intensity probes
medians <- rowMedians(data)
hist_res <- hist(medians, 100, col = "cornsilk1", freq = FALSE, 
                 main = "Histogram of the median intensities", 
                 border = "antiquewhite4",
                 xlab = "Median intensities")
filter <- apply(data, 1, function(x){sum(x > 3.5) >= 3})
data <- subset(data.rma.norm, filter)

## Filter probes that match multiple genes
match <- AnnotationDbi::select(htmg430pm.db,
                               keys = featureNames(data),
                               columns = c("SYMBOL", "ENSEMBL"),
                               keytype = "PROBEID")
match <- subset(match, !is.na(ENSEMBL))
match_grouped <- group_by(match, PROBEID)
match_summarized <- dplyr::summarize(match_grouped, no_of_matches = n_distinct(SYMBOL))
match_filtered <- filter(match_summarized, no_of_matches > 1)
ids_to_exclude <- (featureNames(data) %in% match_filtered$PROBEID)
data <- subset(data, !ids_to_exclude)
match_summarized <- dplyr::summarize(match_grouped, no_of_matches = n_distinct(ENSEMBL))
match_filtered <- filter(match_summarized, no_of_matches > 1)
ids_to_exclude <- (featureNames(data) %in% match_filtered$PROBEID)
data <- subset(data, !ids_to_exclude)
fData(data)$PROBEID <- rownames(fData(data))
fData(data) <- left_join(fData(data), match)
rownames(fData(data)) <- fData(data)$PROBEID
fData(data) <- fData(data)[,-1]

## Prepare names of factors
kainate <- names
kainate <- gsub(".*KA.*", "KA", kainate)
kainate <- gsub(".*CTR.*", "C", kainate)


time <- names
time <- gsub("Input ", "", time)
time <- gsub("s.*", "", time)
time <- gsub(" hour", "h", time)
time <- gsub(" day", "d", time)

## Run contrasts
factor <- factor(paste(kainate,time,sep="."))
design <- model.matrix(~ 0 + factor)
colnames(design) <- levels(factor)
fit <- lmFit(data, design)
cont.matrix <- makeContrasts(`8h` = KA.8h-C.8h,
                             `14d` = KA.14d-C.14d,
                             levels = design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)
`8h` <- topTable(fit2, p.value = 0.05, coef = 1, number = Inf)
`14d` <- topTable(fit2, p.value = 0.05, coef = 2, number = Inf)

## Write DE results to file
write.csv(`8h`, "GSE122228/8h.csv")
write.csv(`14d`, "GSE122228/14d.csv")

##----------------------------------------------------------------------------------

## GSE138100
## Load count data
gse <- read.csv("GSE138100/GSE138100_counts_gencode.M16.txt", sep = "\t") ## Data stored locally
rownames(gse) <- gse$GeneID
gse <- gse[,-1]

## Create a summarizedExperiment object from the count data
colData <- data.frame(treatment = gsub("K.*", "KA", gsub("S.*", "CON", colnames(gse))),
                      accession = paste0("GSM", 4099316:4099327),
                      cell.type = gsub(".*p.*", "Neuron", gsub(".*n.*", "Glial", colnames(gse))),
                      row.names = colnames(gse))
colData$accession <- as.factor(colData$accession)
colData$treatment <- as.factor(colData$treatment)
colData$cell.type <- as.factor(colData$cell.typ)
se <- SummarizedExperiment(as.matrix(gse), colData = colData)

## Filter out low counts
groupsize <- 3
mincount <- 10
keep <- rowSums(assay(se) >= mincount) >= groupsize
se <- se[keep,]

## Create a DESeqDataSet
dds <- DESeqDataSet(se, ~ treatment + cell.type)
dds$group <- factor(paste(dds$treatment, dds$cell.type, sep="."))
design(dds) <- ~ group
dds <- DESeq(dds)

## Differential Gene Expression and then annotate with AnnotationDbi
reslist <- list()
reslist[["Glial"]] <- results(dds, contrast = c("group", "KA.Glial", "CON.Glial"))
reslist[["Neuron"]] <- results(dds, contrast = c("group", "KA.Neuron", "CON.Neuron"))

for(i in 1:length(reslist)) {
  mm <- org.Mm.eg.db
  keys <- gsub("\\..*", "", rownames(reslist[[i]]))
  match <- AnnotationDbi::select(mm,
                                 keys = keys,
                                 columns = c("ENSEMBL", "SYMBOL"),
                                 keytype = "ENSEMBL")
  reslist[[i]]$symbol <- match$SYMBOL[match(keys, match$ENSEMBL)] 
  reslist[[i]] <- reslist[[i]][!is.na(reslist[[i]]$padj),]
  sigres <- reslist[[i]][reslist[[i]]$padj<0.05,] ## Only significantly differentially expressed genes
  sigres <- sigres[order(sigres$log2FoldChange, decreasing = T),]
  reslist[[i]] <- sigres

  ## Write DE result to file
  write.csv(reslist[[i]], paste0("GSE138100/", names(reslist)[i], ".csv"))
}

##----------------------------------------------------------------------------------------------

## GSE213393
## Load count data
gse <- read.csv("GSE213393/GSE213393_geneCounts.tsv", sep = "\t") ## Data stored locally
rownames(gse) <- gse$Geneid
gse <- gse[,-1]

## Create a summarizedExperiment object from the count data
colnames(gse) <- gsub("Sample_", "", colnames(gse))
colData <- data.frame(location = rep(c("C", "I"), 24),
                      treatment = rep(c(rep("KA", 4), rep("CON", 4)), 6),
                      time = rep(c("72h", "168h", "336h"), each = 16),
                      accession = paste0("GSM6585", 453:500),
                      animal = gsub("DC", "", gsub("DI", "", colnames(gse))),
                      sex = rep(c("M", "F"), 6, each = 2),
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

## Create a DESeqDataSet
dds <- DESeqDataSet(se, ~ treatment + location + time)
dds$group <- factor(paste(dds$treatment, dds$time, dds$location, sep="."))
design(dds) <- ~ group
dds <- DESeq(dds)

## Differential Gene Expression and then annotate with AnnotationDbi
reslist <- list()
reslist[["I3d"]] <- results(dds, contrast = c("group", "KA.72h.I", "CON.72h.I"))
reslist[["I7d"]] <- results(dds, contrast = c("group", "KA.168h.I", "CON.168h.I"))
reslist[["I14d"]] <- results(dds, contrast = c("group", "KA.336h.I", "CON.336h.I"))
reslist[["C3d"]] <- results(dds, contrast = c("group", "KA.72h.C", "CON.72h.C"))
reslist[["C7d"]] <- results(dds, contrast = c("group", "KA.168h.C", "CON.168h.C"))
reslist[["C14d"]] <- results(dds, contrast = c("group", "KA.336h.C", "CON.336h.C"))

for(i in 1:length(reslist)) {
  mm <- org.Mm.eg.db
  keys <- gsub("\\..*", "", rownames(reslist[[i]]))
  match <- AnnotationDbi::select(mm,
                                 keys = keys,
                                 columns = c("ENSEMBL", "SYMBOL"),
                                 keytype = "ENSEMBL")
  reslist[[i]]$symbol <- match$SYMBOL[match(keys, match$ENSEMBL)] 
  reslist[[i]] <- reslist[[i]][!is.na(reslist[[i]]$padj),]
  sigres <- reslist[[i]][reslist[[i]]$padj<0.05,] ## Only significantly differentially expressed genes
  sigres <- sigres[order(sigres$log2FoldChange, decreasing = T),]
  reslist[[i]] <- sigres

  ## Write DE results to file
  write.csv(reslist[[i]], paste0("GSE213393/", names(reslist)[i], ".csv"))
}
