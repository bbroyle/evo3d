.onAttach <- function(libname, pkgname) {
  if (!requireNamespace("msa", quietly = TRUE)) {
    packageStartupMessage("Optional package 'msa' not found. For alignment features, install it via BiocManager::install('msa')")
  }
}