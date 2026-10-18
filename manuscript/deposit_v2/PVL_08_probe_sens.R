source(file.path(getwd(), "_env_paths.R"))  # paths from .env (single source of truth)
suppressMessages({library(data.table); library(limma)})

## ---------- parse GSE26549 matrix ----------
parse_gse_matrix <- function(gz){
  lines <- readLines(gz)
  beg <- which(grepl("!series_matrix_table_begin", lines))[1] + 1
  end <- which(grepl("!series_matrix_table_end",   lines))[1]
  block <- lines[beg:end]
  hdr <- gsub('"','', strsplit(block[1], "\t")[[1]])
  m <- fread(text = block[-1], sep = "\t", header = FALSE, fill = TRUE)
  rn <- as.character(m[[1]]); m <- as.matrix(m[, -1, with = FALSE])
  rownames(m) <- rn; colnames(m) <- hdr[-1]; m
}
expr49 <- parse_gse_matrix(file.path(RAW, "GSE26549_series_matrix.txt.gz"))

## ---------- annotation GPL6244 ----------
annot_lines <- readLines(file.path(RAW, "GPL6244.annot.gz"))
cand <- grep("Gene symbol", annot_lines, ignore.case = TRUE)
hdr_i <- cand[!grepl("^#", annot_lines[cand])][1]
ann <- fread(text = annot_lines[hdr_i:length(annot_lines)], sep = "\t", header = TRUE, fill = TRUE)
sym_col <- which(grepl("Gene symbol", names(ann), ignore.case = TRUE))[1]
sym <- sapply(strsplit(as.character(ann[[sym_col]]), " /// "), function(x) x[1])
gi <- data.frame(probe = as.character(ann[[1]]), gene = sym, stringsAsFactors = FALSE)
gi <- gi[!is.na(gi$gene) & gi$gene != "", ]
mt <- expr49[rownames(expr49) %in% gi$probe, , drop = FALSE]
gi <- gi[match(rownames(mt), gi$probe), ]

collapse <- function(mt, gi, fun){
  genes <- unique(gi$gene)
  out <- sapply(genes, function(g){
    idx <- which(gi$gene == g)
    if (length(idx) == 1) mt[idx, ] else apply(mt[idx, , drop = FALSE], 2, fun)
  })
  out <- as.matrix(out)
  colnames(out) <- genes
  rownames(out) <- colnames(mt)   # sample ids
  out                            # samples x genes
}
x49_max  <- collapse(mt, gi, max)
x49_mean <- collapse(mt, gi, mean)
cat("collapsed dims (max):", dim(x49_max), " (mean):", dim(x49_mean), "\n")

## ---------- baseline progression DEG under each method ----------
clin <- fread(file.path(MET, "GSE26549_clinical.csv"))
baseline <- clin[`time of biopsy` == "biopsy at baseline"]
y <- as.integer(baseline$outcome == "progressing"); names(y) <- as.character(baseline$GSM)

run_deg <- function(mat){
  X <- t(mat[as.character(baseline$GSM), , drop = FALSE])   # genes x samples
  design <- model.matrix(~ y)
  fit <- lmFit(X, design); fit <- eBayes(fit)
  tt <- topTable(fit, coef = 2, number = Inf, sort.by = "none")
  tt$gene <- rownames(tt); tt
}
deg_max  <- run_deg(x49_max)
deg_mean <- run_deg(x49_mean)

sig_max  <- unique(deg_max[deg_max$adj.P.Val < 0.05 & abs(deg_max$logFC) > 0.5, ]$gene)
sig_mean <- unique(deg_mean[deg_mean$adj.P.Val < 0.05 & abs(deg_mean$logFC) > 0.5, ]$gene)
cat("DEG (adj.p<0.05,|logFC|>0.5): max =", length(sig_max), " mean =", length(sig_mean), "\n")

## overlap of DEG sets
inter <- intersect(sig_max, sig_mean)
jac <- length(inter) / length(union(sig_max, sig_mean))
cat("DEG overlap: intersect =", length(inter), " Jaccard =", round(jac, 3), "\n")

## logFC concordance on shared genes
shared <- intersect(rownames(deg_max), rownames(deg_mean))
r <- cor(deg_max[shared, "logFC"], deg_mean[shared, "logFC"], method = "spearman")
cat("logFC Spearman (max vs mean, shared genes):", round(r, 3), "\n")

## ---------- external validation under each method (GSE227919 severity gradient) ----------
lines <- readLines(file.path(RAW, "GSE227919_2023_20_03_expr_mat_raw.txt.gz"))
hdr <- strsplit(trimws(lines[1]), "\\s+")[[1]]
fm <- fread(text = lines[-1], sep = " ", header = FALSE, fill = TRUE)
genes <- as.character(fm[[1]]); exm <- as.matrix(fm[, -1, with = FALSE])
storage.mode(exm) <- "double"; rownames(exm) <- genes; colnames(exm) <- hdr
clin9 <- fread(file.path(MET, "GSE227919_clinical.csv"))
ord <- setNames(c(0,1,2,3), c("Control","HkNR","Dysplasia","Cancer"))
clin9$ord <- ord[clin9$Class]
y9 <- clin9$ord; names(y9) <- clin9$Sample_ID
exm <- exm[, as.character(clin9$Sample_ID), drop = FALSE]

grad <- function(sig){
  cf <- intersect(sig, rownames(exm))
  w <- setNames(deg_max[cf, "logFC"], cf)
  Z <- scale(t(exm[cf, , drop = FALSE]))
  score <- as.numeric(Z %*% w)
  cor.test(score, y9, method = "spearman")$estimate
}
rho_max  <- grad(sig_max)
rho_mean <- grad(sig_mean)
cat("external severity-gradient rho: max-method signature =", round(rho_max, 3),
    " mean-method signature =", round(rho_mean, 3), "\n")

write.csv(data.frame(
  metric = c("DEG_max","DEG_mean","overlap","Jaccard","logFC_spearman",
             "ext_rho_max","ext_rho_mean"),
  value = c(length(sig_max), length(sig_mean), length(inter), round(jac,3),
            round(r,3), round(rho_max,3), round(rho_mean,3))),
  file.path(RES, "PVL8_probe_sensitivity.csv"), row.names = FALSE)
cat("DONE Phase8\n")
