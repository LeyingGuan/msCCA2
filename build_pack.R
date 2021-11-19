require(Rcpp)
require(RcppArmadillo)
require(tools)
#RcppArmadillo.package.skeleton("msCCA")
compileAttributes(pkgdir = ".", verbose = TRUE)
package_native_routine_registration_skeleton(".", character_only = FALSE)

pkgbuild::compile_dll(path = ".")
devtools::document(pkg = ".")
devtools::check(pkg = ".")
devtools::document(pkg = ".")

R CMD build msCCA
R CMD install msCCA