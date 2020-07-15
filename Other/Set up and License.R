library(installr)
updateR()

install.packages("Rcpp", repos="https://rcppcore.github.io/drat")

getOption("repos")
options(repos = "https://cran.rstudio.com")
packrat::snapshot()

packageVersion("rsconnect")

library(available)
available("dmf")

library(usethis)

usethis::use_agpl3_license("dynamicfit.app")