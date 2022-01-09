using = function(...) {
  libs = unlist(list(...))
  req = unlist(lapply(libs, require, character.only = TRUE))
  need = libs[req == FALSE]
  if (length(need) > 0) {
    install.packages(need, repos = "http://cran.us.r-project.org")
    lapply(need, require, character.only = TRUE)
  }
}

using(
  "rmarkdown",
  "conflicted",
  "plyr",
  "doParallel",
  "reshape",
  "vars",
  "data.table",
  "Matrix",
  "optimr",
  "trustOptim",
  "systemfit",
  "haven",
  "ggrepel",
  "expint",
  "latex2exp",
  "knitr",
  "flextable",
  "lmtest",
  "broom",
  "purrr",
  "doFuture",
  "tidyverse"
)
conflict_prefer("filter", "dplyr")
conflict_prefer("mutate", "dplyr")
conflict_prefer("arrange", "dplyr")
conflict_prefer("select", "dplyr")