#devtools setup
# devtools::create('dcevalr')

#packages to use: imports
devtools::use_package('igraph')
devtools::use_package('foreach')
devtools::use_package('plyr')
devtools::use_package('stringr')
devtools::use_package('reshape2')
devtools::use_package('methods')

#packages to use: suggests
devtools::use_package('EBcoexpress', type = 'Suggests') #requireNamespace("EBcoexpress", quietly = TRUE)
devtools::use_package('EBarrays', type = 'Suggests') #requireNamespace("EBcoexpress", quietly = TRUE)
devtools::use_package('GeneNet', type = 'Suggests')
devtools::use_package('COSINE', type = 'Suggests')
devtools::use_package('mclust', type = 'Suggests')
devtools::use_package('minqa', type = 'Suggests')
devtools::use_package('testthat', type = 'Suggests')
devtools::use_package('SummarizedExperiment', type = 'Suggests')
devtools::use_package('Biobase', type = 'Suggests')
devtools::use_package('BiocStyle', type = 'Suggests')

#packages to use: enhances
devtools::use_package('parallel',  type = 'Enhances')
devtools::use_package('doSNOW',  type = 'Enhances')
devtools::use_package('doParallel',  type = 'Enhances')

#documentation OR Ctrl + Shift + Alt + R to add skeleton and Ctrl + Shift + B
devtools::document() # Ctrl + Shift + D
devtools::test() # Ctrl + Shift + T

# #vignettes: to setup vignettes
# devtools::use_vignette('dcevalr_vignette')

# #testing: to setup testing
# devtools::use_testthat()

# #to setup the data folder
# devtools::use_data(sim102)

#formatting
# formatR::tidy_dir('R')
# lintr::lint_package()
