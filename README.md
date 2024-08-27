# msCCA

## Installation

Install from github:

install_github("LeyingGuan/msCCA2")

Download the source:

install.packages("msCCA2_1.0.tar.gz")

## Tutorial and Example
### R6 msCCAl1 object
The most flexible use of the implemented method is through the R6class object, msCCAl1, which allows easy sequential exploration of the data. We explain the most critical arguments and functions here.
#### Object initialization
##### Important initialization quantities:

X:   be a list of matrix (of the same sample size) to determine the constructin function of the leading multi-block sparse CCA via L1-norm constraints
beta_init: the initial list of projection coefficients for different data matrices
init_method: beta initialization method if beta_init is NULL. Takes value in ("rgcca", "pma", "convex", "soft-thr"). The "convex" choice uses convex relaxation, which provides theoretical guarantees under stringent assumptions and can be very slow for large-scale data set. "soft-thr" is the suggested version that provides no-worse empirical performance compared to "convex" but much faster.
l1norm_max: largest value of L1 norm allowed (at step 0), default is the minimum between sqrt(p) and sqrt(n/4).
l1norm_min: smallest value of L1 norm considered, default is sqrt(2).
eta: as specified in the manuscript, which is a scale-free parameter determining the "relative" step size during proximal gradient descent. Default value is $1/sqrt(n)$.
eta_ratio: the quantity multiplying eta in the bound decaying description in the manuscript. Default value is $sqrt(1/self$n)$.
rho_maxit: maximum number of proximal gradient iterations.
rho_tol: early stop when the ratio between the L2-norm of the beta change and eta is smaller than rho_tol. Default is 1E-3.

##### Initializing msCCAproximal_l1 to be the msCCAl1 object.
In our experiments, all hyper-parameters have been set as the default values.
```ruby
msCCAproximal_l1 = msCCAl1$new(X = xlist, beta_init =NULL, init_method = "soft-thr",  l1norm_max=NULL, l1norm_min = NULL,
                                 eta = sqrt(n), eta_ratio = sqrt(n),rho_tol = 1E-3, rho_maxit = 5E3)

```

