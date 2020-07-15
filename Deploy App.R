library(rsconnect)

rsconnect::setAccountInfo(name='melissagwolf', token='CDF073365AC304DFB687B4A44BDF9790', secret='WW8+SCx0PvWyZ+bDo1uks2L4ywfGt9aD0Xd33Jcf')

rsconnect::deployApp('C:/users/missg/OneDrive/Documents/Dynamic Model Fit/dynamic/WebsiteCopy.Rmd')

packageVersion("lavaan")

devtools::install_version("lavaan", version = "0.6-5", repos = "http://cran.us.r-project.org")