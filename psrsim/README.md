
# psrsim

<!-- badges: start -->
<!-- badges: end -->

## Installation

To be able to compile the C++ code in the package, you need to first install the [GNU Scientific Library (GSL)](https://www.gnu.org/software/gsl/). On *nix, you can do so using the following commands.

```
wget https://mirror.ibcp.fr/pub/gnu/gsl/gsl-latest.tar.gz

tar xvzf gsl-latest.tar.gz

cd gsl-2.8  # or whatever the subdirectory name is

./configure
make
sudo make install
```

You can install the development version of psrsim like so:

``` r
remotes::install_github("dalejbarr/psrsim")
```

