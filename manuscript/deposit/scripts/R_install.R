## PVL project: install R/Bioconductor stack in WSL (behind SOCKS5 proxy, default repos)
Sys.setenv(http_proxy  = "socks5://127.0.0.1:8428",
           https_proxy = "socks5://127.0.0.1:8428",
           ftp_proxy   = "socks5://127.0.0.1:8428")
options(download.file.method = "libcurl")
## Use DEFAULT repos (cran.r-project.org + bioconductor.org) reachable via proxy
options(repos = c(CRAN = "https://cran.r-project.org"))
## Do NOT set BioC_mirror to a China mirror (unreachable via this proxy)

cat("=== test connectivity ===\n")
cat(R.version.string, "\n")
ok <- tryCatch({ download.file("https://cran.r-project.org/src/contrib/PACKAGES",
                               tempfile(), quiet=TRUE); "CRAN OK" },
               error=function(e) paste("CRAN FAIL:", e$message))
cat(ok, "\n")

cat("=== installing BiocManager ===\n")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
cat("BiocManager version:", as.character(packageVersion("BiocManager")), "\n")

cran_pkgs <- c("tidyverse","data.table","R.utils","glmnet","patchwork","circlize",
               "plyr","jsonlite","openxlsx","pheatmap","viridis","ggpubr","furrr",
               "future","cowplot","ggrepel","purrr","stringr","readr","RColorBrewer")
bioc_pkgs <- c("Biobase","GEOquery","limma","edgeR","DESeq2","clusterProfiler",
               "org.Hs.eg.db","Seurat","ComplexHeatmap","WGCNA","SingleR",
               "celldex","CellChat","harmony","scCustomize","Rtsne","umap")

cat("=== installing packages (Ncpus=4) ===\n")
BiocManager::install(c(cran_pkgs, bioc_pkgs), ask = FALSE, update = FALSE, Ncpus = 4)
cat("=== INSTALL SCRIPT DONE ===\n")
