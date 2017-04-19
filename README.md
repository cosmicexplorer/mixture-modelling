mixture-modelling
=================

Implementation of EM to fit (multivariate) mixture models to data using any user-specified metric or parameterized distribution. I was using [mixtools](https://cran.r-project.org/web/packages/mixtools/index.html), but my data was a frame of latitude and longitude, and mixtools only supports the standard Euclidean metric.

Drawing implementation notes from [this post](https://tinyheero.github.io/2016/01/03/gmm-em.html) (thanks!).
