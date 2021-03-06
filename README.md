# rcITR

**rcITR** is an R package for estimating risk-controlled individualized treatment rules (rcITR). An rcITR seeks to maximize the expected benefit (or reward) while controlling risk at a pre-specified threshold. Two modeling approaches are available:

- rcDT (Risk-Controlled Decision Tree): Single decision tree model
- rcRF (Risk-Controlled Random Forest): Ensemble of rcDT models 

Detailed documentation can be found <a href="https://kdoub5ha.github.io/rcITR/rcITR-vignette" target="_top">here</a>

# Installation
```{r, echo = TRUE, eval = FALSE, warning = FALSE, error = FALSE, message = FALSE}
devtools::install_github("kdoub5ha/rcITR")
```
