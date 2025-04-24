.onAttach <- function(libname, pkgname) {
  packageStartupMessage("Note: 'msa' is required for alignment functionality. Install via BiocManager::install('msa')")
}