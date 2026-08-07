#!/usr/bin/env Rscript
## Focused Bioconductor install for Phase 2 enrichment (clusterProfiler family).
Sys.setenv(http_proxy  = "socks5://127.0.0.1:8428",
           https_proxy = "socks5://127.0.0.1:8428")
options(download.file.method = "libcurl", Ncpus = 4)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("clusterProfiler","DOSE","fgsea","enrichplot","org.Hs.eg.db"),
                     ask = FALSE, update = FALSE, dependencies = TRUE)
cat("\n==> BIOC PHASE2 INSTALL DONE\n")
