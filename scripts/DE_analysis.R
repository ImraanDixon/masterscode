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
library(AnnotationDbi)

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
