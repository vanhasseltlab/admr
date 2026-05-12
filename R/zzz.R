.onAttach <- function(libname, pkgname) {
  packageStartupMessage(
    "WARNING: 'admr' is deprecated and no longer maintained.\n",
    "Please use 'admixr2' instead: https://github.com/LeidenPharmacology/admixr2\n",
    "Install with: devtools::install_github(\"LeidenPharmacology/admixr2\")"
  )
}
