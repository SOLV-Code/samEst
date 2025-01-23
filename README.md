# samEst <img src="man/figures/logo.png" align="right" height="138" />
Collection of estimation models to estimate stock-recruitment parameters for salmon populations


## Install instructions
 samEst depends on the most recent version of Rstan that is not yet available on CRAN. To install it, run:

```{r} 
remove.packages(c("StanHeaders", "rstan"))
install.packages("StanHeaders", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
install.packages("rstan", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
```

Then install samEst:

```{r} 
install.packages("remotes") 
remotes::install_git('https://github.com/Pacific-salmon-assess/samEst')
```
