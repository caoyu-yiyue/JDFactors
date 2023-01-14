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

## Outline of Mathematical Background

Some mathematical Backgrounds for ARMA-GARCH-Copula here, for more clarity.

### ARMA-GARCH

* ARMA model

$$
r_{i, t}=\mu_{i, t}+\epsilon_{i, t}=\delta_{i}+\sum_{j=1}^{p} \psi_{i, j}\left(r_{i, t-j}-\delta_{i}\right)+\sum_{k=1}^{q} \lambda_{i, k} \epsilon_{i, t-k}+\epsilon_{i, t}
$$

* NGARCH model using here

$$
\sigma_{i, t}^{2}=\omega_{i}+\alpha_{i} \sigma_{i, t-1}^{2}\left(z_{i, t-1}-\eta_{i}\right)^{2}+\beta_{i} \sigma_{i, t-1}^{2}
$$

### Copula

#### Copula Funtion

For a collection of random variables $\mathbf{Y}=\left(Y_{1}, \ldots, Y_{n}\right)$, the joint distribution $F$ with margin distributions $F_1, \ldots, F_n$ can be written as:
  
$$
F\left(y_{1}, \ldots, y_{n}\right)=C\left(F_{1}\left(y_{1}\right), \ldots, F_{n}\left(y_{n}\right)\right)
$$

Function  $C:[0,1]^{n} \rightarrow[0,1]$ is the copula function. If the all variables are continuous, there'll be only one copula function. And for every variable, if it is continuous, doing probability integral transforms to it, $i \in \{1, \ldots, n\}, U_i = F_i(Y_i)$, $U_i$ will be standard uniform distribution. So the copula function can be transformed to the form with $\mathcal{U}(0,1)$:

$$
C\left(u_{1}, \ldots, u_{n}\right)=F\left(F_{1}^{-1}\left(u_{1}\right), \ldots, F_{n}^{-1}\left(u_{n}\right)\right)
$$

For $F^{-1}_i$ is the margin quantile function. The joint density function of $\mathbf{Y}$ is:

$$
f\left(y_{1}, \ldots, y_{n}\right)=\prod_{i=1}^{n} f_{i}\left(y_{i}\right)\cdot \frac{\partial C\left(F_{1}\left(y_{1}\right), \ldots, F_{n}\left(y_{n}\right)\right)}{\partial F_{1}\left(y_{1}\right)\ldots \partial F_{n}\left(y_{n}\right)}=\prod_{i=1}^{n} f_{i}\left(y_{i}\right)\cdot c\left(F_{1}\left(y_{1}\right), \ldots, F_{n}\left(y_{n}\right)\right)
$$

So the copula function and the joint distribution can be separated. The form of copula function may differ, for normal copula:

$$
C^{Norm}\left(u_{1}, \ldots, u_{n} ; \mathbf{P}\right)=\Phi_{\mathbf{P}}\left(\Phi^{-1}\left(u_{1}\right), \ldots, \Phi^{-1}\left(u_{n}\right)\right)
$$

And for t copula:

$$
C^{T}\left(u_{1}, \ldots, u_{n} ; v, \mathbf{P}\right)=T_{v, \mathbf{P}}\left(T_{v}^{-1}\left(u_{1}\right), \ldots, T_{v}^{-1}\left(u_{n}\right)\right)
$$

#### Static Copula

For example, a normal copula at time $t$:

$$
F\left(r_{1, t}, \ldots, r_{5, t} \mid \mathcal{F}_{t-1} ; \mathbf{P}\right)=\Phi_{\mathbf{P}}\left(\Phi^{-1}\left(u_{1, t}\right), \ldots, \Phi^{-1}\left(u_{5, t}\right)\right)
$$

For $u_{i, t}=F_{i}\left(r_{i, t} \mid \mathcal{F}_{t-1} ;\boldsymbol{\theta_{m, i}}\right)$, and correlation matrix $\mathbf{P}$ is composed by $\mathbf{z}_{\mathbf{t}}^{*}=\left(\Phi^{-1}\left(u_{1, t}\right), \ldots, \Phi^{-1}\left(u_{5, t}\right)\right)^{\prime}$, defined as copula shocks.

#### Dynamic Copula

For a collections of variables $\mathbf{R}_{\mathbf{t}}=\left(R_{1, t}, \ldots, R_{n, t}\right), \{t=1,2, \ldots, T\}$, with $\boldsymbol{\mu_{t}}=E\left[\mathbf{y}_{\mathrm{t}} \mid \mathcal{F}_{t-1}\right]$ and positive definite conditional variance-covariance matrix $\boldsymbol{\Sigma_{t}}=E\left[\left(\mathbf{r}_{\mathrm{t}}-\boldsymbol{\mu_{t}}\right)\left(\mathbf{r}_{\mathrm{t}}-\boldsymbol{\mu_{t}}\right)^{\prime} \mid \mathcal{F}_{t-1}\right]=E\left[\boldsymbol{\epsilon_{t}} \boldsymbol{\epsilon_{t}}^{\prime} \mid \mathcal{F}_{t-1}\right]$, a MGARCH model:

$$
\mathbf{r_{t}}=\boldsymbol{\mu_{t}}+\boldsymbol{\epsilon_{t}}=\boldsymbol{\mu_{t}}+\boldsymbol{\Sigma_{t}}^{\frac{1}{2}} \mathbf{n}_{\mathrm{t}}
$$

with $\mathbf{n}_{\mathbf{t}} \sim I I D\left(\mathbf{0}, \mathbf{I}_{\mathbf{n}}\right)$. $\boldsymbol{\Sigma_t}^\frac{1}{2}$ is the Cholesky Decomposition of $\boldsymbol{\Sigma_t}$. Discompose the variance-covariance matrix:

$$
\boldsymbol{\Sigma_{t}}=\mathbf{V}_{\mathbf{t}}^{\frac{1}{2}} \mathbf{P}_{\mathbf{t}} \mathbf{V}_{\mathbf{t}}^{\frac{1}{2}}
$$

$\mathbf{V_t}$ is a diagonal matrix of $\sigma^2_{i,t}$,  $\mathbf{P_t}$ is a positive definition correlation matrix.
 
Now define a standard shock, $\mathbf{z_t} = \mathbf{V_t}^{-\frac{1}{2}}\boldsymbol{\epsilon_t}$,  $\mathbf{P_t}$ is the conditional variance-covarince matrix of $\mathbf{z_t}$. Conditional correlation matrix have the following process:

$$
\mathbf{Q}_{\mathrm{t}}=(1-a-b) \bar{\mathbf{Q}}+a \mathbf{z}_{\mathrm{t}-1} \mathbf{z}_{\mathrm{t}-1}^{\prime}+b \mathbf{Q}_{\mathrm{t}-1}
$$

$\mathbf{\bar{Q}}$ is the long perioid (unconditional correlation matrix), a and b are non-negative, and a + b > 1. Transform $\mathbf{Q_t}$ as following to make the diagonal to 1.

$$
\mathbf{P_t}=\operatorname{diag}\left(\mathbf{Q_t}\right)^{-\frac{1}{2}} \mathbf{Q_t} \operatorname{diag}\left(\mathbf{Q_t}\right)^{-\frac{1}{2}}
$$

Rewrite the process as:

$$
\mathbf{Q_{t}}=\mathbf{\bar{Q}}+a\left(\mathbf{z_{t-1}} \mathbf{z_{t-1}}^{\prime}-\mathbf{\bar{Q}}\right)+b\left(\mathbf{Q_{t-1}}-\mathbf{\bar{Q}}\right)
$$

For ADCC, add a asymmetrical parameter:

$$
\mathbf{Q}_{\mathbf{t}}=\left(\mathbf{\bar{Q}}-a \mathbf{\bar{Q}}-b \mathbf{\bar{Q}}-g \mathbf{\bar{Q}}^{-}\right)+a \mathbf{z}_{\mathbf{t}-1} \mathbf{z}_{\mathbf{t}-1}^{\prime}+b \mathbf{Q}_{\mathbf{t}-1}+g \mathbf{z}_{\mathbf{t}-1}^{-} \mathbf{z}_{\mathbf{t}-\mathbf{1}}^{-\prime}
$$

with $z^{-}_{i,t-1} = min(z_{i,t-1}, 0)$, and $\mathbf{\bar{Q}^{-}} = E[\mathbf{z^{-}_{t-1}}\mathbf{z^{-\prime}_{t-1}}]$.
