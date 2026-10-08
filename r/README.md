# foodnet for R

The R companion package of food**net**: it receives the parameters of a consumer-resource model (CRM) from
the food**net** page and shapes them for a simulator, with
[miaSim](https://bioconductor.org/packages/release/bioc/html/miaSim.html) as the example. It assumes no
simulator: everything it returns is a plain matrix, vector or list.

## Install

```r
install.packages(c("remotes", "BiocManager"))
remotes::install_github("hallucigenia-sparsa/foodnet", subdir = "r", ref = "v0.2.0")
BiocManager::install("miaSim")                          # for the simulations
```

Install the version of your food**net** (the page and its help print the line for yours): the package and
the page talk to each other, and a newer package refuses an older page.

If that fails with "HTTP error 404" or "cannot open URL", either the release it names is not published yet, or
a GitHub token stored on this machine is being used and cannot see the repository. Installing from a clone needs
no GitHub access: `remotes::install_local("<the repository>/r")`.

The page's **CRM example** button sets up five gut species for which the steps below run as they are; the
page's help walks through them ("Simulate with miaSim, step by step").

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
| `crm_efficiency(crm)` | miaSim's efficiency matrix E, in each taxon's own unit of abundance: positive for a resource taken up, negative for a by-product, the opposite of foodnet's signed matrix (built from the consumed and produced matrices, never from that one) |
| `crm_scale(crm)`, `crm_unscale(crm, abundance, args = args)` | that unit, in the growth curve's unit, and simulated abundances back in the curve's unit |
| `crm_subset(crm, taxa, resources)` | the same parameters for fewer taxa or resources |
| `crm_chemistry(crm)` | each resource's formula and charge (ChEBI), carbon atoms and degree of reduction |
| `crm_electron_balance(crm)` | each taxon's electrons out over electrons in, with its range over the replicates: a check, which changes no value |
| `crm_phase(crm, "stationary")` | the stationary phase, when the search asked for both phases |
| `crm$interval_start`, `crm$interval_end` | taxa by resources, the hours each value was measured over |
| `crm$consumed_lower`, `crm$consumed_upper`, `crm$produced_lower`, `crm$produced_upper` | each amount's lowest and highest replicate (mM); a 0 from 0 to the detection limit at least |
| `crm_backcheck(crm, monod_constant)` | each taxon simulated alone, against its own monoculture |
| `as_miasim(crm, x0, monod_constant)` | the arguments of `miaSim::simulateConsumerResource` |
| `crm_readme(crm)`, `crm_write(crm, dir)` | the README food**net** wrote, and the parameters as files |

```r
# NA cells are refused unless you say how; here, as 0
crm_backcheck(crm, monod_constant = 1, na = "zero")
args <- as_miasim(crm, x0 = crm$biomass_start, monod_constant = 1, missing_resource = 0, na = "zero")
tse <- do.call(miaSim::simulateConsumerResource, c(args, list(t_end = 48, t_store = 480)))
abundance <- crm_unscale(crm, SummarizedExperiment::assay(tse), args = args)
```

A run is deterministic as `as_miasim()` shapes it. In miaSim 1.18, `simulateConsumerResource` adds a random
immigrant at rate `migration_p` (default 0.01) even with `stochastic = FALSE`, so `as_miasim()` passes
`migration_p = 0`; drift, epochs and external events are off unless `stochastic = TRUE`, and measurement
noise unless `error_variance > 0`. To explore noise, change them in the list before the call, since it already
holds `migration_p` (`args$migration_p <- 0.01; args$stochastic <- TRUE`); to turn it all off again, set
`args$migration_p <- 0; args$stochastic <- FALSE; args$error_variance <- 0`. Leave `norm = FALSE`: relative abundances cannot
be turned back by `crm_unscale()`.

miaSim (1.18, `consumerResourceModel`) grows a taxon by its growth rate times the sum over resources of
`E * R / (R + K)`, takes up each resource at `R / (R + K)` per unit of abundance whatever the size of E,
and makes each by-product at `|E|` times its growth. It has no uptake rate, so the unit of abundance
decides how fast a taxon eats: in cells/mL a culture would empty its medium within minutes, whatever its
growth rate. foodnet gives each taxon a unit of its own, `crm_scale(crm)` = `dx * C / (mu * S)` growth
curve units (biomass change `dx`, growth rate `mu`, total uptake `C`, `S` the sum of the squared uptakes).
In it, E is each consumed resource's share of the uptake (they sum to 1, so a saturated taxon grows at its
measured rate and a trace substrate weighs little) and `-produced * C / S` on each by-product, so a taxon
alone gains its measured biomass and makes its measured by-products in proportion to what it takes up.
Below saturation it grows slower than measured. The same data in another unit give the same simulation. `as_miasim` puts
the starting abundances (in each growth curve's unit, `crm$biomass_unit`) into these units, and
`crm_unscale` turns simulated ones back. How much a taxon takes up of each resource follows the Monod
constants, which foodnet does not measure; `crm_backcheck()` shows how close a choice comes, and how long
each taxon takes to grow against its phase. `crm_efficiency(crm, scale = "shares")` is foodnet 0.1.0's
matrix (shares of uptake), for models that read E that way.

`as_miasim` stops when a taxon has no growth rate, a resource no starting concentration, or when starting
abundances or Monod constants are not given, rather than passing a number nobody measured (miaSim would
draw them at random). It refuses parameters pooled across media or from the stationary phase unless
`allow =` names them. `crm_efficiency` warns how many NA cells of each kind it counts as 0. miaSim keeps
`t_store` of its 0.1 h steps; `t_store = 10 * t_end` keeps them all, so the run reaches `t_end`.
