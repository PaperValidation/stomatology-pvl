suppressMessages({library(data.table); library(enrichR)})
RES <- "~/PVL/results"

## ---------- load discovery signature (101 genes) ----------
deg  <- fread(file.path(RES, "GSE26549_progression_DEG.csv"))
sig  <- unique(deg[adj.P.Val < 0.05 & abs(logFC) > 0.5]$gene); sig <- sig[!is.na(sig)]
cat("signature genes for enrichment:", length(sig), "\n")

## ---------- enrichment libraries ----------
libs <- c("GO_Biological_Process_2021", "KEGG_2021_Human",
          "MSigDB_Hallmark_2020", "ChEA_2022", "ENCODE_TF_ChIP-seq_2015")
enr <- enrichR::enrichr(sig, libs)
saveRDS(enr, file.path(RES, "PVL7_enrichment_raw.rds"))

## functional modules (GO / KEGG / Hallmark): top terms
go  <- enr[["GO_Biological_Process_2021"]]; go  <- go[go$Adjusted.P.value < 0.05, ]
kegg<- enr[["KEGG_2021_Human"]];          kegg<- kegg[kegg$Adjusted.P.value < 0.05, ]
hal <- enr[["MSigDB_Hallmark_2020"]];     hal <- hal[hal$Adjusted.P.value < 0.05, ]
cat("GO terms (adj.p<0.05):", nrow(go), " KEGG:", nrow(kegg), " Hallmark:", nrow(hal), "\n")

## upstream TF (ChEA / ENCODE)
chea <- enr[["ChEA_2022"]];               chea <- chea[chea$Adjusted.P.value < 0.05, ]
enc  <- enr[["ENCODE_TF_ChIP-seq_2015"]]; enc  <- enc[enc$Adjusted.P.value < 0.05, ]
cat("ChEA_2022 TFs (adj.p<0.05):", nrow(chea), " ENCODE TF ChIP (adj.p<0.05):", nrow(enc), "\n")

## flag PVL/HNSCC-relevant master regulators
relevant <- c("TP63","NOTCH1","TP53","MYC","EGFR","STAT3","NFKB1","RELA","JUN","FOS",
              "ETS1","SPI1","CEBPA","CEBPB","SOX2","KLF4","GRHL3","IRF1","IRF7","STAT1")
chea_hits <- chea[chea$Term %in% relevant, ]
enc_hits  <- enc[enc$Term %in% relevant, ]
cat("PVL-relevant TFs in ChEA_2022:", if(nrow(chea_hits)) paste(chea_hits$Term, collapse=",") else "NONE", "\n")
cat("PVL-relevant TFs in ENCODE ChIP:", if(nrow(enc_hits)) paste(enc_hits$Term, collapse=",") else "NONE", "\n")

## ---------- outputs ----------
## FDR-filtered (stringent)
write.csv(head(go[, 1:6], 30),   file.path(RES, "PVL7_GO_fdr.csv"), row.names = FALSE)
write.csv(head(kegg[, 1:6], 30),  file.path(RES, "PVL7_KEGG_fdr.csv"), row.names = FALSE)
write.csv(head(hal[, 1:6], 30),   file.path(RES, "PVL7_Hallmark_fdr.csv"), row.names = FALSE)
## top by combined score (ranking, regardless of FDR) -- informative for small gene set
top_cs <- function(d, n = 15) head(d[order(-d$Combined.Score), c("Term","Adjusted.P.value","Combined.Score","Genes")], n)
write.csv(top_cs(go),   file.path(RES, "PVL7_GO_top.csv"), row.names = FALSE)
write.csv(top_cs(kegg),  file.path(RES, "PVL7_KEGG_top.csv"), row.names = FALSE)
write.csv(top_cs(hal),   file.path(RES, "PVL7_Hallmark_top.csv"), row.names = FALSE)
write.csv(top_cs(chea),  file.path(RES, "PVL7_ChEA_TF_top.csv"), row.names = FALSE)
write.csv(top_cs(enc),   file.path(RES, "PVL7_ENCODE_TF_top.csv"), row.names = FALSE)

summ <- data.frame(
  metric = c("sig_n","GO_adjp05","KEGG_adjp05","Hallmark_adjp05",
             "ChEA_TF_adjp05","ENCODE_TF_adjp05",
             "relevant_TF_ChEA","relevant_TF_ENCODE"),
  value = c(length(sig), nrow(go), nrow(kegg), nrow(hal),
            nrow(chea), nrow(enc),
            if(nrow(chea_hits)) nrow(chea_hits) else 0,
            if(nrow(enc_hits)) nrow(enc_hits) else 0))
write.csv(summ, file.path(RES, "PVL7_summary.csv"), row.names = FALSE)
cat("DONE Phase7\n")
