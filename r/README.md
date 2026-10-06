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

Printing `crm` shows what to read before simulating, every time: the cells never assayed, inconclusive or
without a phase, the links seen only in another medium (so their size is unknown), the values whose
experiments disagree, the taxa without a growth rate, and parameters from pooled media or the stationary
phase.

## Use them

| function | returns |
|---|---|
| `crm_consumed(crm)`, `crm_produced(crm)` | taxa by resources, the amounts in mM; NA is never zero |
| `crm_rates(crm)` | the growth rates (1/h); `missing =` fills the ones nobody measured, only if you say so |
| `crm_resources(crm)` | the initial concentrations of the medium (mM) |
| `crm_efficiency(crm)` | miaSim's efficiency matrix E: yields (biomass per mM taken up) positive, by-products per unit of growth negative |
| `crm_backcheck(crm, monod_constant)` | each taxon simulated alone, against its own monoculture |
| `as_miasim(crm, x0, monod_constant)` | the arguments of `miaSim::simulateConsumerResource` |
| `crm_readme(crm)`, `crm_write(crm, dir)` | the README food**net** wrote, and the parameters as files |

```r
crm_backcheck(crm, monod_constant = 1)
args <- as_miasim(crm, x0 = crm$biomass_start, monod_constant = 1)
set.seed(1)
tse <- do.call(miaSim::simulateConsumerResource, c(args, list(t_end = 48, t_store = 480)))
```

miaSim (1.18, `consumerResourceModel`) grows a taxon by its growth rate times the sum over resources of
`E * R / (R + K)`, takes up each resource at `R / (R + K)` per unit of abundance whatever the size of E,
and makes each by-product at `|E|` times its growth. So `crm_efficiency()` sets a positive entry to the
taxon's biomass change divided by its growth rate and its total uptake (a yield), and a negative one to
the amount produced times the growth rate divided by the biomass change. A taxon simulated alone then gains
its measured biomass and makes its measured by-products in proportion to what it takes up; how much it
takes up of each resource follows the Monod constants, which foodnet does not measure, and
`crm_backcheck()` shows how close a choice comes. Abundances are in each taxon's growth curve unit
(`crm$biomass_unit`), so the starting abundances must be too. `crm_efficiency(crm, scale = "shares")` is
foodnet 0.1.0's matrix (shares of uptake), for models that read E that way.

`as_miasim` stops when a taxon has no growth rate, a resource no starting concentration, or when starting
abundances or Monod constants are not given, rather than passing a number nobody measured (miaSim would
draw them at random). It refuses parameters pooled across media or from the stationary phase unless
`allow =` names them. `crm_efficiency` warns how many NA cells of each kind it counts as 0. miaSim keeps
`t_store` of its 0.1 h steps; `t_store = 10 * t_end` keeps them all, so the run reaches `t_end`.
