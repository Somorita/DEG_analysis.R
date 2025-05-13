# DEG_analysis.R

# Load the count data
urld <- "https://www.ncbi.nlm.nih.gov/geo/download/?format=file&type=rnaseq_counts"
path <- paste(urld, "acc=GSE224537", "file=GSE224537_raw_counts_GRCh38.p13_NCBI.tsv.gz", sep="&")
tbl <- as.matrix(data.table::fread(path, header=T, colClasses="integer"), rownames="GeneID")

# Load annotation
apath <- paste(urld, "type=rnaseq_counts", "file=Human.GRCh38.p13.annot.tsv.gz", sep="&")
annot <- data.table::fread(apath, header=T, quote="", stringsAsFactors=F, data.table=F)
rownames(annot) <- annot$GeneID

# Subset tbl to only Control and DMOG samples
tbl_dmog <- tbl[, 1:4]

# Define sample groups manually
sml <- c("Control", "Control", "DMOG", "DMOG")
gs <- factor(sml)
coldata <- data

# Filter low expressed genes
keep <- rowSums(tbl_dmog >= 10) >= min(table(gs))
tbl_dmog <- tbl_dmog[keep, ]

# Create DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = tbl_dmog, colData = coldata, design = ~ condition)

#run Deseq
dds_dmog <- DESeq(dds, test = "Wald")

# Get DEGs
res_dmog <- results(dds_dmog, contrast = c("condition", "DMOG", "Control"), alpha = 0.05,
                    pAdjustMethod ="fdr", lfcThreshold = 1)
res_dmog <- as.data.frame(res_dmog)
res_dmog$GeneID <- rownames(res_dmog)
res_dmog <- merge(res_dmog, annot[, c("GeneID", "Symbol")], by = "GeneID")
write.table(res_dmog, file = "DEGs_DMOG_vs_Control.txt", sep = "\t", quote = FALSE, 
            row.names = FALSE)
