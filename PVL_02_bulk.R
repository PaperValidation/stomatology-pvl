#!/usr/bin/env Rscript
source(file.path(getwd(), "_env_paths.R"))  # paths from .env (single source of truth)
## ============================================================
## PVL Phase 2: Bulk transcriptomic DEG + pathway enrichment
##   Cohort A (discovery): GSE26549  - OPL microarray, progression outcome
##   Cohort B (validation): GSE227919 - RNA-seq gradient + progression
## ============================================================
Sys.setenv(http_proxy  = "socks5://127.0.0.1:8428",
           https_proxy = "socks5://127.0.0.1:8428",
           ftp_proxy   = "socks5://127.0.0.1:8428")
options(download.file.method = "libcurl")

suppressMessages({
  library(limma); library(edgeR); library(ggplot2)
  library(enrichR); library(pheatmap); library(GEOquery)
  library(data.table); library(dplyr); library(tidyr); library(stringr)
  library(RColorBrewer); library(ComplexHeatmap); library(circlize)
  library(httr)
})
# enrichR talks to the Enrichr API -> route through the SOCKS5 proxy
tryCatch(set_config(use_proxy("127.0.0.1", 8428, scheme = "socks5")),
         error = function(e) invisible(NULL))

dir.create(RES, showWarnings = FALSE)
cat("=== PVL Phase 2 start ===\n")

## ---------- helpers ----------
parse_gse_matrix <- function(gz){
  con <- gzfile(gz); lines <- readLines(con); close(con)
  begin <- which(grepl("!series_matrix_table_begin", lines))
  end   <- which(grepl("!series_matrix_table_end",   lines))
  block <- lines[(begin+1):(end-1)]
  hdr   <- gsub('"', '', strsplit(block[1], "\t")[[1]])
  dat   <- read.delim(text = block[-1], sep = "\t", header = FALSE,
                      stringsAsFactors = FALSE, check.names = FALSE)
  probes <- dat[,1]; mat <- as.matrix(dat[,-1]); rownames(mat) <- probes
  colnames(mat) <- hdr[-1]; mat
}

volcano_plot <- function(deg, file, title, top=10){
  deg$neglogP <- -log10(pmax(deg$P.Value, 1e-300))
  deg$lab <- ifelse(deg$neglogP > 10 & abs(deg$logFC) > 1, deg$gene, "")
  p <- ggplot(deg, aes(x=logFC, y=neglogP, color=signif)) +
    geom_point(size=0.7, alpha=0.6) + theme_minimal() +
    scale_color_manual(values=c("grey","red","blue")) +
    labs(title=title, x="log2 FC", y="-log10 P") +
    geom_text(aes(label=lab), size=2.5, vjust=-0.5, check_overlap=TRUE) +
    theme(legend.position="none")
  ggsave(file, p, width=7, height=5)
}

run_enrich <- function(gup, gdn, tag){
  dbs <- c("GO_Biological_Process_2021","KEGG_2021_Human","MSigDB_Hallmark_2020")
  for (dirn in c("up","down")) {
    g <- if (dirn=="up") gup else gdn
    if (length(g) < 5) { cat("  skip enrich", tag, dirn, "(n<5)\n"); next }
    g <- unique(g[!is.na(g) & g!=""])
    cat("  enrichR", tag, dirn, "n=", length(g), "\n")
    res <- tryCatch(enrichr(g, dbs), error=function(e){ cat("  enrichr ERR:", e$message, "\n"); NULL })
    if (is.null(res)) next
    for (db in names(res)) {
      df <- res[[db]]
      if (is.null(df) || nrow(df)==0) next
      safename <- gsub("[ ,/]", "_", db)
      write.csv(df, file.path(RES, sprintf("enrich_%s_%s_%s.csv", safename, tag, dirn)), row.names=FALSE)
      tt <- head(df[order(-df$Combined.Score), ], 15)
      tt$Term <- factor(tt$Term, levels=rev(tt$Term))
      p <- ggplot(tt, aes(x=Combined.Score, y=Term, fill=-log10(.data$Adjusted.P.value))) +
        geom_bar(stat="identity") + theme_minimal() +
        labs(title=sprintf("%s | %s | %s", db, tag, dirn), x="Combined Score", y="") +
        scale_fill_gradient(low="blue", high="red")
      ggsave(file.path(RES, sprintf("enrich_%s_%s_%s.png", safename, tag, dirn)), p, width=9, height=6)
    }
  }
}

## ============================================================
## Cohort A: GSE26549  (HuGene 1.0 ST, progression outcome)
## ============================================================
cat("\n--- Cohort A: GSE26549 ---\n")
matA <- parse_gse_matrix(file.path(RAW, "GSE26549_series_matrix.txt.gz"))
clinA <- fread(file.path(MET, "GSE26549_clinical.csv"))
clinA_base <- clinA[`time of biopsy` == "biopsy at baseline"]
cat("baseline samples:", nrow(clinA_base),
    "| progressing:", sum(clinA_base$outcome=="progressing"),
    "| non_progressing:", sum(clinA_base$outcome=="non_progressing"), "\n")

keep <- intersect(colnames(matA), clinA_base$GSM)
matA_b <- matA[, keep, drop=FALSE]
clinA_b <- clinA_base[match(keep, GSM)]
matA_b <- matA_b[, match(clinA_b$GSM, colnames(matA_b))]

# probe -> gene symbol via local GPL6244.annot.gz (getGPL removed in recent GEOquery)
annot <- file.path(RAW, "GPL6244.annot.gz")
alines <- readLines(gzfile(annot))
abegin <- which(grepl("!platform_table_begin", alines))
ahdr <- strsplit(alines[abegin+1], "\t")[[1]]
abody <- read.delim(text = alines[(abegin+2):length(alines)], sep = "\t",
                    header = FALSE, quote = "", stringsAsFactors = FALSE,
                    check.names = FALSE, fill = TRUE)
colnames(abody) <- ahdr
gpltab <- abody[, c("ID","Gene symbol")]
gpltab <- gpltab[gpltab$`Gene symbol` != "" & gpltab$`Gene symbol` != "---", ]
gpltab$sym <- sapply(strsplit(gpltab$`Gene symbol`, " /// "), function(x) x[1])
gpltab <- gpltab[gpltab$sym != "", ]
matA_sym <- matA_b[rownames(matA_b) %in% gpltab$ID, , drop=FALSE]
gpltab2 <- gpltab[match(rownames(matA_sym), gpltab$ID), ]
# collapse duplicates: max expression per symbol
collapsed <- sapply(unique(gpltab2$sym), function(s){
  idx <- which(gpltab2$sym == s)
  if(length(idx)==1) matA_sym[idx,] else apply(matA_sym[idx,,drop=FALSE], 2, max)
})
collapsed <- t(collapsed)
cat("gene-level matrix:", nrow(collapsed), "genes x", ncol(collapsed), "samples\n")

grp <- factor(clinA_b$outcome, levels=c("non_progressing","progressing"))
design <- model.matrix(~ grp)
fit <- lmFit(collapsed, design)
fit <- eBayes(fit, trend=TRUE)
topA <- topTable(fit, coef=2, number=Inf, sort.by="P")
topA$gene <- rownames(topA)
topA$signif <- ifelse(topA$adj.P.Val<0.05 & abs(topA$logFC)>1, "yes", "no")
write.csv(topA, file.path(RES, "GSE26549_progression_DEG.csv"), row.names=FALSE)
volcano_plot(topA, file.path(RES,"GSE26549_progression_volcano.png"),
             "GSE26549 progression (progressing vs non)")
degA_up <- topA$gene[topA$adj.P.Val<0.05 & topA$logFC>1]
degA_down<- topA$gene[topA$adj.P.Val<0.05 & topA$logFC< -1]
cat("GSE26549 DEG up/down:", length(degA_up), "/", length(degA_down), "\n")
run_enrich(degA_up, degA_down, "GSE26549_progression")

## ============================================================
## Cohort B: GSE227919 (RNA-seq, log2 expr, 4-class gradient)
## ============================================================
cat("\n--- Cohort B: GSE227919 ---\n")
raw <- readLines(gzfile(file.path(RAW, "GSE227919_2023_20_03_expr_mat_raw.txt.gz")))
hdr <- strsplit(raw[1], " ")[[1]]
body <- read.table(text=raw[-1], sep=" ", header=FALSE, stringsAsFactors=FALSE,
                   check.names=FALSE)
geneids <- body[,1]; exprB <- as.matrix(body[,-1])
rownames(exprB) <- make.unique(as.character(geneids))
colnames(exprB) <- hdr
clinB <- fread(file.path(MET, "GSE227919_clinical.csv"))
keepB <- intersect(hdr, clinB$Sample_ID)
exprB <- exprB[, keepB, drop=FALSE]
clinB <- clinB[match(keepB, Sample_ID)]
exprB <- exprB[, match(clinB$Sample_ID, colnames(exprB))]
cat("exprB:", nrow(exprB), "genes x", ncol(exprB), "samples\n")
cat("Class:", toString(table(clinB$Class)), "\n")

# contrast 1: Cancer vs Control
cl1 <- clinB$Class %in% c("Control","Cancer")
design1 <- model.matrix(~ 0 + Class, data=clinB[cl1,])
colnames(design1) <- gsub("Class","",colnames(design1))
fit1 <- lmFit(exprB[,cl1], design1); fit1 <- eBayes(fit1, trend=TRUE)
c1 <- makeContrasts(Cancer - Control, levels=design1)
fit1c <- contrasts.fit(fit1, c1); fit1c <- eBayes(fit1c)
topB_cc <- topTable(fit1c, number=Inf); topB_cc$gene <- rownames(topB_cc)
topB_cc$signif <- ifelse(topB_cc$adj.P.Val<0.05 & abs(topB_cc$logFC)>1,"yes","no")
write.csv(topB_cc, file.path(RES,"GSE227919_Cancer_vs_Control_DEG.csv"), row.names=FALSE)
volcano_plot(topB_cc, file.path(RES,"GSE227919_Cancer_vs_Control_volcano.png"),
             "GSE227919 Cancer vs Control")
up_cc <- topB_cc$gene[topB_cc$adj.P.Val<0.05 & topB_cc$logFC>1]
dn_cc <- topB_cc$gene[topB_cc$adj.P.Val<0.05 & topB_cc$logFC< -1]
cat("GSE227919 Cancer vs Control up/down:", length(up_cc),"/",length(dn_cc),"\n")
run_enrich(up_cc, dn_cc, "GSE227919_CancerVsControl")

# contrast 2: progression among PML (Progressed vs Stable)
prog <- clinB[Progression_status %in% c("Stable","Progressed-SCC","Progressed-Dys")]
prog$prog_bin <- ifelse(grepl("Progressed", prog$Progression_status), "Progressed","Stable")
exprP <- exprB[, colnames(exprB) %in% prog$Sample_ID]
prog <- prog[match(colnames(exprP), Sample_ID)]
designP <- model.matrix(~ prog_bin, data=prog)
fitP <- lmFit(exprP, designP); fitP <- eBayes(fitP, trend=TRUE)
topP <- topTable(fitP, coef=2, number=Inf); topP$gene <- rownames(topP)
topP$signif <- ifelse(topP$adj.P.Val<0.05 & abs(topP$logFC)>1,"yes","no")
write.csv(topP, file.path(RES,"GSE227919_progression_DEG.csv"), row.names=FALSE)
volcano_plot(topP, file.path(RES,"GSE227919_progression_volcano.png"),
             "GSE227919 progression (Progressed vs Stable)")
up_p <- topP$gene[topP$adj.P.Val<0.05 & topP$logFC>1]
dn_p <- topP$gene[topP$adj.P.Val<0.05 & topP$logFC< -1]
cat("GSE227919 progression up/down:", length(up_p),"/",length(dn_p),"\n")
run_enrich(up_p, dn_p, "GSE227919_progression")

## gradient trend (Control<HkNR<Dysplasia<Cancer) among ordered classes
ord <- as.character(clinB$Class)
scores <- as.numeric(recode(ord, Control="0", HkNR="1", Dysplasia="2", Cancer="3"))
keepOrd <- !is.na(scores)
designT <- model.matrix(~ scores[keepOrd])
fitT <- lmFit(exprB[,keepOrd], designT); fitT <- eBayes(fitT, trend=TRUE)
topT <- topTable(fitT, coef=2, number=Inf); topT$gene <- rownames(topT)
write.csv(topT, file.path(RES,"GSE227919_gradient_trend_DEG.csv"), row.names=FALSE)
cat("GSE227919 gradient trend significant:", sum(topT$adj.P.Val<0.05), "\n")

cat("\n=== PVL Phase 2 DONE ===\n")
