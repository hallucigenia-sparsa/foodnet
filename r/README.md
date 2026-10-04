# foodnet for R

The R companion package of food**net**: it receives the parameters of a consumer-resource model (CRM) from
the food**net** page and shapes them for a simulator, with
[miaSim](https://bioconductor.org/packages/release/bioc/html/miaSim.html) as the example. It assumes no
simulator: everything it returns is a plain matrix, vector or list.

## Install

```r
install.packages("remotes")
remotes::install_github("hallucigenia-sparsa/foodnet", subdir = "r")
```

If that fails with "HTTP error 404", a GitHub token stored on this machine is being used and cannot see the
repository. Installing from a clone needs no GitHub access: `remotes::install_local("<the repository>/r")`.

## Receive the parameters

```r
library(foodnet)
crm <- foodnet_listen()     # then press Send to R (under Get CRM parameters) on the foodnet page
```

or, without opening a port, fetch them from the address the page shows:
`crm <- foodnet_crm("http://127.0.0.1:PORT/crm.json?token=...")`. The `crm.json` in a downloaded parameter
zip works too.

Printing `crm` shows what to read before simulating, every time: the cells never assayed, the links seen
only in another medium (so their size is unknown), the values that pool experiments that disagree, and the
taxa without a growth rate.

## Use them

| function | returns |
|---|---|
| `crm_consumed(crm)`, `crm_produced(crm)` | taxa by resources, the amounts in mM; NA is never zero |
| `crm_rates(crm)` | the growth rates (1/h); `missing =` fills the ones nobody measured, only if you say so |
| `crm_resources(crm)` | the initial concentrations of the medium (mM) |
| `crm_efficiency(crm)` | the efficiency matrix E: consumed shares positive, by-products per unit taken up negative |
| `as_miasim(crm)` | the arguments of `miaSim::simulateConsumerResource` |
| `crm_readme(crm)`, `crm_write(crm, dir)` | the README food**net** wrote, and the parameters as files |

```r
args <- as_miasim(crm)
tse <- do.call(miaSim::simulateConsumerResource, c(args, list(t_end = 48, t_store = 100)))
```

`as_miasim` stops when a taxon has no growth rate or a resource no starting concentration, rather than
passing a number nobody measured. `crm_efficiency` warns, naming them, when it counts links seen only in
another medium as 0. How to scale E is a modeling choice, so nothing does it for you.
