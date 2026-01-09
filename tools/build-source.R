#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
package_root <- if (length(args) > 0) args[[1]] else "."

tar_options <- "--format=posix --no-same-owner --numeric-owner --owner=0 --group=0"
if (nzchar(Sys.getenv("TAR_OPTIONS"))) {
  tar_options <- paste(Sys.getenv("TAR_OPTIONS"), tar_options)
}
Sys.setenv(TAR_OPTIONS = tar_options)

r_bin <- file.path(R.home("bin"), "R")
status <- system2(r_bin, c("CMD", "build", package_root))

if (!identical(status, 0L)) {
  stop("R CMD build failed with status: ", status)
}
