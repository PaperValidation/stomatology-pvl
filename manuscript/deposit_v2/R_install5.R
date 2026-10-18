#!/usr/bin/env Rscript
## Install remaining CRAN missing packages as PPM BINARY (no compilation).
Sys.setenv(http_proxy  = "socks5://127.0.0.1:8428",
           https_proxy = "socks5://127.0.0.1:8428")
suppressMessages(library(pak))
options(repos = c(CRAN = "https://packagemanager.posit.co/cran/__linux__/noble/latest"))
pkgs <- c("ggpubr","rstatix","car","lme4","nloptr","Deriv","doBy","pbkrtest")
cat("==> pak (binary) CRAN:", paste(pkgs, collapse=", "), "\n")
pak::pkg_install(pkgs, ask = FALSE)
cat("\n==> PAK CRAN DONE\n")
