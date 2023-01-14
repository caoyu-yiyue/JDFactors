# JDFactors

A repository of finance calculating, for research of modeling joint dynamic correlation of 5 factors model. The name of JDFactor is just Joint Dynamic Factors.

It's now a pure R project, but with some early python codes for the old version. By the definition of command line API, it's easy to build the result using GNU Make by one line of code:

```shell
make all
```

For the result with the robust test, run:

```shell
make all_verbose
```

But there's no original data here, for the copyright reason. And there're some Chinese comments in the code, for the early days of my quick commenting. I'm sorry that now I don't have time to translate it into English for more people.

## Target

Frankly speaking, the codes do the following works:

1. An ARMA model for every factor, and finding the best arma order by investigating orders from 1 to 10.
2. ARMA-GARCH model fit for every factor, using the arma order selected.
3. GARCH-Copula model with the norm, t copula, combining DCC and ADCC copula.
4. Rolling fit for investment performance test. Refit in expanding window every 30 days, and filter model by new data. Determine weights of factors by optimizing effective function, and calculating portfolio return.
5. Robost test by different risk tolerant coefficients and different data frequencies.

Of course more complicated works in detail. Maybe some are not efficient codes early days, but I'll be happy if I can help others facing the same problems.

