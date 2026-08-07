#!/usr/bin/env Rscript
## Install Bioconductor enrichment stack as PPM BINARY packages via pak (avoids fgsea boost compile error).
Sys.setenv(http_proxy  = "socks5://127.0.0.1:8428",
           https_proxy = "socks5://127.0.0.1:8428")
suppressMessages(library(pak))
repos <- c(
  CRAN = "https://packagemanager.posit.co/cran/__linux__/noble/latest",
  Bioc = "https://packagemanager.posit.co/bioc/__linux__/noble/latest"
)
options(repos = repos)
pkgs <- c("clusterProfiler","DOSE","fgsea","enrichplot","ggtree","treeio","org.Hs.eg.db",
           "ggpubr","rstatix","car","lme4","nloptr","Deriv","doBy","pbkrtest")
cat("==> pak installing (binary) :", paste(pkgs, collapse=", "), "\n")
pak::pkg_install(pkgs, ask = FALSE)
cat("\n==> PAK BIOC INSTALL DONE\n")
