<!-- README.md is generated from README.Rmd. Please edit that file -->

# OSplines

<!-- badges: start -->
<!-- badges: end -->

The goal of OSplines is to efficiently implement model-based smoothing with the integrated Wiener's process, within a variety of Bayesian hierarchical models.
## Installation

You can install the development version of this package from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("Smoothing-IWP/OSplines")
```

## Example

This is a basic example which shows you how to use `OSplines` to fit and
analyze some models, we consider the following data set of COVID-19 mortality in Canada, which is available in the package:

``` r
options(warn=-1)
library(tidyverse)
library(OSplines)
head(covid_canada)
```

We can fit a model with $\text{IWP}_3(\sigma)$ prior using the function `model_fit`:
```{r warning=FALSE}
fit_result <- model_fit(new_deaths ~ weekdays1 + weekdays2 + weekdays3 + weekdays4 + weekdays5 + weekdays6 +
                          f(smoothing_var = t, model = "IWP", order = 3, k = 30), 
                        data = covid_canada, method = "aghq", family = "Poisson")
```

We can take a look at the posterior summary of this model:
```{r}
summary(fit_result)
```

We can also see the inferred function $f$:
```{r warning=FALSE}
plot(fit_result)
```


We can use the `predict` function to obtain the posterior summary of $f$ or its derivative at `new_data`.

For the function $f$:
```{r warning=FALSE}
predict_f <- predict(fit_result, variable = "t", newdata = data.frame(t = seq(from = 605, to = 617, by = 0.1)))
predict_f %>% ggplot(aes(x = x)) + geom_line(aes(y = mean), lty = "solid") +
  geom_line(aes(y = plower), lty = "dashed") +
  geom_line(aes(y = pupper), lty = "dashed") +
  theme_classic()
```

For the first derivative:
```{r warning=FALSE}
predict_f1st <- predict(fit_result, variable = "t", newdata = data.frame(t = seq(from = 605, to = 617, by = 0.1)), degree = 1)
predict_f1st %>% ggplot(aes(x = x)) + geom_line(aes(y = mean), lty = "solid") +
  geom_line(aes(y = plower), lty = "dashed") +
  geom_line(aes(y = pupper), lty = "dashed") +
  theme_classic()
```


For the second derivative:
```{r warning=FALSE}
predict_f2nd <- predict(fit_result, variable = "t", newdata = data.frame(t = seq(from = 605, to = 617, by = 0.1)), degree = 2)
predict_f2nd %>% ggplot(aes(x = x)) + geom_line(aes(y = mean), lty = "solid") +
  geom_line(aes(y = plower), lty = "dashed") +
  geom_line(aes(y = pupper), lty = "dashed") +
  theme_classic()
```
