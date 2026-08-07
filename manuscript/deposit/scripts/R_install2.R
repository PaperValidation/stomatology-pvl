#!/usr/bin/env Rscript
# Recovery install: use Posit Public Package Manager (PPM) Linux binaries for Ubuntu 24.04 (noble)
# to avoid source compilation (nloptr/Deriv/lme4 etc. failed before due to missing system libs).
options(repos = c(
  CRAN        = "https://packagemanager.posit.co/cran/__linux__/noble/latest",
  Bioconductor = "https://packagemanager.posit.co/bioc/__linux__/noble/latest"
))
Sys.setenv(http_proxy  = "socks5://127.0.0.1:8428",
           https_proxy = "socks5://127.0.0.1:8428")
options(download.file.method = "libcurl",
        pkgType = "binary",
        Ncpus = 4,
        HTTPUserAgent = sprintf("R/%s R (%s)", getRversion(),
                                paste(getRversion(), R.version$platform, R.version$arch, R.version$os)))

cat("==> repos:\n"); print(getOption("repos"))

# 1) CRAN packages that failed (mostly compiled deps)
cran_missing <- c("ggpubr","rstatix","car","lme4","nloptr","Deriv","doBy","pbkrtest","ggtree","treeio")
cat("\n==> installing CRAN packages (binary)...\n")
install.packages(cran_missing, dependencies = TRUE, type = "binary")

# 2) Bioconductor enrichment + annotation
bioc_missing <- c("clusterProfiler","DOSE","fgsea","enrichplot")
cat("\n==> installing Bioconductor packages (binary)...\n")
install.packages(bioc_missing, dependencies = TRUE, type = "binary")

# 3) remotes for CellChat (GitHub)
cat("\n==> installing remotes (binary)...\n")
install.packages("remotes", type = "binary")

# 4) CellChat from GitHub (needs R>=3.6; jinworks/CellChat)
cat("\n==> installing CellChat from GitHub (source)...\n")
Sys.unsetenv("http_proxy","https_proxy")  # github often reachable; keep proxy for safety
Sys.setenv(http_proxy  = "socks5://127.0.0.1:8428",
           https_proxy = "socks5://127.0.0.1:8428")
remotes::install_github("jinworks/CellChat", dependencies = TRUE, upgrade = "never")

cat("\n==> RECOVERY INSTALL DONE\n")
