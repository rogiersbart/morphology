.onAttach <- function(libname, pkgname) {
  std::err("! Package {.pkg morphology} is still in its experimental lifecycle stage.")
  std::err("! Use at your own risk, and submit issues here:")
  std::err("! {.url https://github.com/rogiersbart/morphology/issues}")
  invisible()
}
