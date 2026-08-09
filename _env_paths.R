# =============================================================================
# _env_paths.R  --  single source of truth for dataset paths
# Sourced by all PVL analysis scripts (run from the project root in WSL).
# Reads DATASET_PATH and DATASET_PATH_TEMP from the project .env and converts
# Windows drive paths (G:\, E:\) to WSL mount paths (/mnt/g, /mnt/e).
# =============================================================================

.read_env_value <- function(key) {
  envp <- file.path(getwd(), ".env")
  if (!file.exists(envp)) {
    envp <- "/mnt/c/Users/tidus/PycharmProjects/Stomatology-PVL/.env"  # fallback
  }
  lines <- readLines(envp, warn = FALSE)
  pat <- paste0("^", key, "\\s*=")
  for (l in lines) {
    if (grepl(pat, l)) {
      val <- sub(paste0("^", key, "\\s*=\\s*"), "", l)
      val <- gsub("^['\"]|['\"]$", "", val)   # strip surrounding quotes
      val <- sub("\\\\$", "", val)            # strip trailing backslash
      return(val)
    }
  }
  stop(paste("key not found in .env:", key))
}

.win2wsl <- function(p) {
  m <- regexpr("^[A-Za-z]:", p)
  if (m > 0) {
    drv <- tolower(substr(p, 1, 1))
    p <- paste0("/mnt/", drv, substr(p, 3, nchar(p)))
  }
  p <- gsub("\\\\", "/", p)
  p <- gsub("/+", "/", p)
  p
}

# Canonical data (raw + metadata) -> G:  (DATASET_PATH)
DATASET_PATH      <- .win2wsl(.read_env_value("DATASET_PATH"))
# Temp / intermediate results          -> E:  (DATASET_PATH_TEMP)
DATASET_PATH_TEMP <- .win2wsl(.read_env_value("DATASET_PATH_TEMP"))

RAW <- file.path(DATASET_PATH, "raw")
MET <- file.path(DATASET_PATH, "metadata")
RES <- DATASET_PATH_TEMP   # analysis outputs now live on the temp volume (E:)
