#!/usr/bin/env Rscript
## ============================================================
## PVL Phase 5: LINCS L1000 drug repurposing (signed Connectivity-Map)
## Principle: a candidate "malignancy-reversing" compound should
##   DOWN-regulate progression-UP genes  AND  UP-regulate progression-DOWN genes.
## Implemented via enrichR directional LINCS libraries (Kuleshov 2016; Subramanian 2017).
## NOTE: outputs are *transcriptome-matched candidate compounds*, NOT treatment advice.
## ============================================================
Sys.setenv(http_proxy  = "socks5://127.0.0.1:8428",
           https_proxy = "socks5://127.0.0.1:8428")
suppressMessages({ library(enrichR); library(data.table); library(dplyr); library(ggplot2) })

## 输入数据目录 = .env 的 DATASET_PATH（WSL 下为 /mnt/g/...）；分析结果仍留本地 ~/PVL/results
DATASET_PATH <- "/mnt/g/CMCC Sync/Medical Data/Stomatology/PVL"
RAW  <- file.path(DATASET_PATH, "raw")
RES <- "~/PVL/results"
dir.create(RES, showWarnings = FALSE)
cat("=== PVL Phase 5 (LINCS) start ===\n")

## ---------- signed discovery signature (validated 86-gene set) ----------
deg <- fread(file.path(RES, "GSE26549_progression_DEG.csv"))
sig <- unique(deg[adj.P.Val < 0.05 & abs(logFC) > 0.5]$gene); sig <- sig[!is.na(sig)]
## genes present in GSE227919 too (the cross-platform-validated signature)
lines <- readLines(file.path(RAW, "GSE227919_2023_20_03_expr_mat_raw.txt.gz"))
hdr <- strsplit(trimws(lines[1]), "\\s+")[[1]]
exm_genes <- as.character(fread(text = lines[-1], sep = " ", header = FALSE, fill = TRUE)[[1]])
common_feat <- intersect(sig, exm_genes)
w <- setNames(deg$logFC, deg$gene)[common_feat]; w <- w[!is.na(w)]
up   <- names(w)[w > 0]      # high in progression / malignancy
down <- names(w)[w < 0]      # low in progression (protective / normalized)
cat("signed signature: up =", length(up), " down =", length(down), "\n")

## ---------- helper: run an enrichR query, return tidy table ----------
qry <- function(genes, lib){
  if (length(genes) < 3) return(data.frame())
  res <- tryCatch(enrichr(genes, lib), error = function(e) NULL)
  if (is.null(res) || !(lib %in% names(res))) return(data.frame())
  d <- res[[lib]]
  if (!nrow(d)) return(data.frame())
  d <- as.data.frame(d)
  d$Drug <- sapply(strsplit(d$Term, " "), function(x){
    f <- x[length(x)]
    p <- strsplit(f, "-")[[1]]
    if (length(p) >= 3) paste(p[2:(length(p) - 1)], collapse = "-") else f
  })
  d$query_lib <- lib
  d[, c("Term","Drug","Overlap","P.value","Adjusted.P.value","Combined.Score","Genes","query_lib")]
}

## ---------- signed Connectivity-Map ----------
## protective: drug DOWN-regulates up-genes  (lib = Chem_Pert_down, query=up)
##             drug UP-regulates down-genes (lib = Chem_Pert_up,   query=down)
prot_down <- qry(up,   "LINCS_L1000_Chem_Pert_down")   # suppresses progression genes
prot_up   <- qry(down, "LINCS_L1000_Chem_Pert_up")     # rescues protective genes
## promoting (to flag/avoid): drug UP-regulates up-genes ; DOWN-regulates down-genes
promo_up   <- qry(up,   "LINCS_L1000_Chem_Pert_up")
promo_down <- qry(down, "LINCS_L1000_Chem_Pert_down")

## confirmatory consensus (curated) layer: drug UP-signature overlaps protective(down) genes
cons_prot <- qry(down, "LINCS_L1000_Chem_Pert_Consensus_Sigs")
cons_promo<- qry(up,   "LINCS_L1000_Chem_Pert_Consensus_Sigs")

neglog <- function(p) -log10(pmax(p, 1e-300))

## ---------- build protective candidate table (require BOTH directions) ----------
## A reversal compound must DOWN-regulate progression-UP genes (prot_down) AND
## UP-regulate progression-DOWN genes (prot_up). Intersection = high-confidence.
build_both <- function(a, b){
  if (!nrow(a) || !nrow(b)) return(data.frame())
  aggA <- a[a$Adjusted.P.value < 0.05, ] %>% group_by(Drug) %>%
    summarise(scoreA = sum(neglog(Adjusted.P.value)),
              termsA = paste(unique(Term), collapse = " | "),
              genesA = paste(unique(unlist(strsplit(Genes, ";"))), collapse = ";"),
              nA = n(), .groups = "drop")
  aggB <- b[b$Adjusted.P.value < 0.05, ] %>% group_by(Drug) %>%
    summarise(scoreB = sum(neglog(Adjusted.P.value)),
              termsB = paste(unique(Term), collapse = " | "),
              genesB = paste(unique(unlist(strsplit(Genes, ";"))), collapse = ";"),
              nB = n(), .groups = "drop")
  m <- inner_join(aggA, aggB, by = "Drug")
  m$reversal_score   <- m$scoreA + m$scoreB
  m$up_targets_down  <- m$genesA   # drug DOWN-regulates these progression-UP genes
  m$down_targets_up  <- m$genesB   # drug UP-regulates these progression-DOWN genes
  m[order(-m$reversal_score), ]
}
cand <- build_both(prot_down, prot_up)
promo <- build_both(promo_up, promo_down)

## flag candidates confirmed in consensus layer (normalized name match;
## directional vs consensus LINCS libraries use partly different nomenclature)
if (nrow(cand)){
  nn <- function(s) tolower(gsub("[^a-z0-9]+", "", s))
  cons_drugs_n <- nn(unique(cons_prot$Drug))
  cand$consensus_confirmed <- sapply(cand$Drug, function(d)
    any(grepl(nn(d), cons_drugs_n, fixed = TRUE)))
}

cat("protective reversal candidates (BOTH directions, adj.p<0.05):", nrow(cand), "\n")
cat("  consensus-confirmed:", if (nrow(cand)) sum(cand$consensus_confirmed) else 0, "\n")
cat("promoting (progression-like, BOTH directions):", nrow(promo), "\n")
cat("debug consensus: cons_prot nrow =", nrow(cons_prot),
    " cons_promo nrow =", nrow(cons_promo), "\n")

## ---------- outputs ----------
if (nrow(cand)) write.csv(cand, file.path(RES, "PVL5_lincs_reversal_candidates.csv"), row.names = FALSE)
if (nrow(promo)) write.csv(promo, file.path(RES, "PVL5_lincs_promoting.csv"), row.names = FALSE)
write.csv(data.frame(direction = c("up_in_progression","down_in_progression"),
                     n_genes = c(length(up), length(down)),
                     genes = c(paste(up, collapse = ";"), paste(down, collapse = ";"))),
          file.path(RES, "PVL5_signature_signed.csv"), row.names = FALSE)

## top-20 reversal bar plot
if (nrow(cand)){
  top <- head(cand, 20)
  png(file.path(RES, "PVL5_reversal_top20.png"), width = 900, height = 700, res = 120)
  print(ggplot(top, aes(x = reorder(Drug, reversal_score), y = reversal_score, fill = consensus_confirmed)) +
          geom_bar(stat = "identity") + coord_flip() +
          labs(title = "Top LINCS reversal candidates (PVL progression signature)",
               x = "Compound", y = "Reversal score (-log10 adj.P, summed)") +
          theme_minimal() + scale_fill_manual(values = c("grey70","steelblue")))
  dev.off()
  cat("top 5 reversal candidates:\n")
  print(head(cand[, c("Drug","reversal_score","consensus_confirmed","nA","nB")], 5))
}

cat("\n=== PVL Phase 5 DONE ===\n")
