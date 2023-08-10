# FPSpatioTemp

This repository contains all code for the paper:

__Higham, M.__, Dumelle, M., Hammond, C., Ver Hoef, J. M., & Wells, J. (2023). An application of spatio-temporal modeling to finite population abundance prediction. _Journal of Agricultural, Biological and Environmental Statistics_, 1 - 25. [Link to Abstract](https://link.springer.com/article/10.1007/s13253-023-00565-y).

The data, provided by the Alaska Department of Fish and Game, that is used in the paper can be found at `/data/moose_14_20.rda`.

The vignette in `/vignettes` shows how to fit the model to a particular data set with the functions in `/R`.

To install from Git, run

```{r}
remotes::install_github("highamm/FPSpatioTemp")
```