#!/usr/bin/env Rscript

# Install GitHub packages for methylation analysis

packages_to_install <- list(
  list(repo = "hovestadtlab/conumee2", package = "conumee2", subdir = "conumee2"),
  list(repo = "achilleasNP/IlluminaHumanMethylationEPICanno.ilm10b5.hg38",
       package = "IlluminaHumanMethylationEPICanno.ilm10b5.hg38"),
  list(repo = "mwsill/RFpurify", package = "RFpurify"),
  list(repo = "Xiaoqizheng/InfiniumPurify", package = "InfiniumPurify")
)

if (!require("devtools", quietly = TRUE)) {
  install.packages("devtools", repos = "https://cloud.r-project.org")
}

for (pkg in packages_to_install) {
  if (!require(pkg$package, character.only = TRUE, quietly = TRUE)) {
    message(paste("Installing", pkg$package, "from GitHub:", pkg$repo))
    args <- list(pkg$repo, quiet = TRUE)
    if (!is.null(pkg$subdir)) {
      args$subdir <- pkg$subdir
    }
    do.call(devtools::install_github, args)
  } else {
    message(paste(pkg$package, "already installed"))
  }
}

message("GitHub packages installation complete")
