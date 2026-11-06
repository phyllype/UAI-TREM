# uaitrem 0.5.0

## Speed

* The sparse Cholesky factorization and the selected inverse now use subtrees of the elimination
  tree (Geist and Ng, 1989), and the selected inverse reads the columns of the dense tail by direct
  index. The results are bit for bit the same as before, with any number of threads. In simulation
  (10 iterations of REML, 16 threads, median of 3 runs), a pedigree model of 1 million animals went
  from 36.7 to 11.2 min, and `model_mt()` with 100,000 animals from 6.2 to 3.3 min.
* `pedigree()` no longer clears a full vector for each animal when it calculates the inbreeding. In
  simulation, 1 million animals went from 71 s to 2 s, with the same F bit for bit.
* `model_mt()` uses the ordering and the symbolic factor of the fit again for the residuals, and does
  not calculate them a second time. The results are bit for bit the same. In simulation, 100,000
  animals went from 5.6 to 4.5 min.
* The approximate minimum degree ordering (Amestoy, Davis and Duff, 1996) is available with the
  environment variable `UAITREM_ORDEM=amd`. The default stays the exact minimum degree. This order
  changes the rounding and the point where a fit stops. In our checks the -2logL changed by 9e-11
  (relative) or less, and the components by 5e-5 or less. Near a covariance boundary the likelihood
  is flat, and the components changed by up to 1e-4.

## New

* The new name of the package is `uaitrem`, for UAI-TREM: Unravelling Additive and Indirect
  genetic effects with Time-dependent Random Effect Models (before, BreedingR). In R:
  `library(uaitrem)`; the options are `uaitrem.*`, the environment variables `UAITREM_*` and
  the condition classes `uaitrem_*` (before, `BreedingR.*`, `BREEDINGR_*` and `breedingr_*`).

* `model_mt()` takes the route of `model()` when a covariance group of several terms
  (`group =`) ends near a correlation of +-1: it fits the same terms without the group and then
  the full model from that optimum, and keeps the better by the same rules. The component
  names lead to the terms without the trait (`animal@t1`) and without the coefficient
  (`rn[0]`); a correlation near +-1 inside one term has no reduced model. Gate in
  `tests/testthat/test-reduzido-mt.R` (a direct-maternal model of two traits simulated with no
  genetic variance).
* `gibbs()`: when the mixed model equations lose their Cholesky factor inside a chain of a
  Gaussian trait, the chain no longer stops: the locations are proposed from the equations
  with a small multiple of their diagonal added (the smallest of 1e-12 to 1e-6 of the
  diagonal that factors), and a Metropolis-Hastings step against the true conditional accepts
  the proposal or keeps the state. The ratio is in closed form,
  `log alpha = delta (x_new' D x_new - x_old' D x_old) / (2 var(residual))`, so the step keeps
  the posterior. `message` says how many iterations needed it, how many proposals were
  rejected and the state at the first one. Where the equations are near singular the proposal
  is narrower than the conditional and most proposals are rejected: the chain moves through
  the other components. With `family = "probit"` the same step runs on the locations given the
  liability, and the PX-DA move, which needs the exact solution, is left out in those
  iterations. Gates in `tests/testthat/test-gibbs-piso.R`: with the rescue forced at every
  iteration (delta 0.01 of the diagonal) and the components given, the posterior means of the
  levels fall within 0.1 `sqrt(PEV)` of the BLUP and their posterior SDs within 5% of
  `sqrt(PEV)`; in a probit sire model the posterior means fall within 0.1 posterior SD of the
  exact chain and the SDs within 5% (simulation; with the components estimated and the PX-DA
  move on, the posterior mean of the sire variance was 0.4666 against 0.4646 and 0.4635 of
  two exact chains, Monte Carlo SE 0.001).
* `gibbs()`: the draw of `var(residual)` is truncated at `1e-10 var(y)`, the floor of the
  residual in `model()`, and `message` says how many draws reached it. On a residual that the
  data hardly separate from the other variances (only their sum identified) the chain could
  go down to about 1e-13 of the scale, where the mixed model equations lose their Cholesky
  factor, and stop: in a synthetic case of 100 animals of 12 dams one chain stopped at 12063
  of 19000 kept draws, and it now runs to the end with the same posterior mean. A draw below
  the floor is replaced by a draw of the truncated conditional by its inverse distribution
  function, which gives the truncated conditional exactly and leaves a chain that never
  reaches the floor unchanged. Gates in `tests/testthat/test-gibbs-piso.R`. When a chain
  still stops on the Cholesky factor, the message gives the state at the failure (the
  residual and, by group, the smallest variance and the smallest eigenvalue over the
  largest, on the scale of `var(y)`), and when a group passed 10 `var(y)` it says that a term
  whose levels each sit in one class of a fixed effect is identified between classes only by
  the prior.
* `model()`: a variance whose standard error is 10 times the phenotypic variance or more is
  named in `message` as holding almost no information, beside the SINGULAR verdict, which
  reads the eigenvalues of the information matrix and left it out (in a maternal model,
  var(maternal) of 0.020 with a standard error of 12647). The rule reads the phenotypic scale
  and not the estimate: a small variance with a standard error above it covers zero, and that
  standard error is the information the data hold.
* `model(residual_by =)`: a residual variance per class of a column (an age class, a band of
  days), `var(residual|class)` for each class, estimated with the other components. The
  engine fits the rows scaled by `sqrt(s2_1 / s2_k)`, where the model has one residual, and
  adds the class block: the score of each class from the leverages of its records, its
  vector in the average information, its EM step and the Jacobian of the scaling in the
  -2logL. `h2()` reads the residual at the mean over the records of the class variances (the
  rule of a fit with weights), `h2_curve(bases = list(residual =))` at the class of each
  point, and `residuals()` divides by the variance of the class of the record. Gates in
  `tests/testthat/test-residuo-classe.R`: the -2logL equals the dense V form, which builds R
  directly, to 1e-9; the score equals central finite differences; classes of equal variance
  give the -2logL of `model()` with one residual; a simulated design (planted 0.5, 1 and 2)
  is recovered inside 3 SE, and the SE of the information matrix is within 20% of the one of
  the numerical Hessian. In our own data set (a daily feeding trait, 177336 records of 2988
  animals, a random regression of order 2) six bands of 11 days in the test lowered the
  -2logL by 2638 with 5 more components; the residual variance went from 4584 in the first
  band to 2653 in the fourth and 3365 in the last (SE 23 to 39).
* `model_mt()`: a random term can apply to some of the traits only, with `traits =`.
  The usual case is a trait with repeated records and a trait with one record per animal:
  `cbind(y1, y2) ~ cg + animal(id) + pe(id, traits = "y1")`. Before, `pe()` applied to
  both traits. On the one-record trait the data do not separate it from the residual, so
  the information matrix was singular and every standard error was `NaN`. Now the term has
  no component and no solution on the traits that it does not name. The univariate fits
  of the default start leave it out of those traits, and the other fitters stop on
  `traits =`. The terms of one `group =` can name different traits:
  `rn(id, base = c("phi0", "phi1"), traits = "y1", group = "g") + animal(id, traits =
  "y2", group = "g")` estimates the genetic covariance between the coefficients of a
  random regression on `y1` and the breeding value for a one-record trait `y2`. A fixed
  effect can be on some traits only: `fixed(cg, traits = "y1")` for a class (the new
  marker `fixed(cg)` is the unmarked class `cg`) and `cov(x, traits = "y1")` for a
  covariate. On the other traits the effect has no equation, and `dropped_x` does not
  list it. A column of a term may be missing on a row where none of the traits of the
  term is observed (the contemporary group of a daily trait on the row of a final
  weight): that cell does not enter the likelihood, and the fit puts a value there that
  changes nothing.

  When no row holds both traits (each record in a row of its own), the data do not
  estimate the residual covariance between them. The fit now holds it at 0, gives it no
  standard error and says so in `message`. Before, the zero row of that covariance made
  the information matrix singular, and every standard error of the fit was `NaN`.

  The gates are in `tests/testthat/test-mt-traits.R`: with the covariances between traits
  at 0, the bivariate -2logL equals the sum of the two univariate -2logL to 1e-8, in the
  sparse and in the dense V form; the sparse -2logL equals the dense V form with
  covariances; and the score equals central finite differences in all parameters, with
  the `y2` record in its own row and in the row of a `y1` record. A group of a random
  regression on `y1` and `animal()` on `y2` passes the dense V form and the finite
  differences too, and so does the identity with a fixed class on `y1` only and a
  covariate on `y2` only.

  Recovery, simulation only (`validation/mt_traits_recovery.R`): five generations, 2000
  animals with records, six `y1` records and one `y2` record per animal, 20 replicates.
  All 40 fits converged, 0.8 s each (median).

  | component | true | mean, own row | mean, shared row | empirical SD | mean SE |
  |---|---|---|---|---|---|
  | var(animal@y1) | 4 | 3.958 | 3.960 | 0.41 | 0.45 |
  | cov(animal@y2,animal@y1) | 1.2 | 1.257 | 1.259 | 0.16 | 0.19 |
  | var(animal@y2) | 2 | 2.109 | 2.110 | 0.28 | 0.30 |
  | var(pe@y1) | 3 | 3.127 | 3.127 | 0.24 | 0.28 |
  | var(res@y1) | 5 | 5.005 | 5.005 | 0.08 | 0.07 |
  | var(res@y2) | 4 | 4.030 | 4.030 | 0.19 | 0.23 |

  The coverage of estimate +/- 1.96 SE was 0.90 to 1.00. In the shared-row layout the
  residual covariance (true 0) came out at -0.007, empirical SD 0.11.

  Recovery of a random regression on `y1` with a one-record `y2`, simulation only
  (`validation/mt_traits_rr_recovery.R`): the same five generations, a Legendre regression
  of order 1 on six `y1` records for the genetic and the permanent effect, `rn()` on `y1`
  and `animal()` on `y2` in one group, a fixed slope on `y1` only, 20 replicates. All 20
  fits converged, 2.5 s each (median).

  | component | true | mean | empirical SD | mean SE | coverage |
  |---|---|---|---|---|---|
  | var(rn[0]@y1) | 4 | 4.060 | 0.52 | 0.48 | 0.95 |
  | cov(rn[1]@y1,rn[0]@y1) | 0.5 | 0.516 | 0.23 | 0.18 | 0.95 |
  | cov(animal@y2,rn[0]@y1) | 1.2 | 1.197 | 0.16 | 0.20 | 1.00 |
  | var(rn[1]@y1) | 1 | 0.987 | 0.13 | 0.14 | 0.95 |
  | cov(animal@y2,rn[1]@y1) | -0.4 | -0.438 | 0.10 | 0.11 | 0.90 |
  | var(animal@y2) | 2 | 1.987 | 0.30 | 0.28 | 0.90 |
  | var(pe[0]@y1) | 3 | 3.017 | 0.32 | 0.31 | 0.95 |
  | cov(pe[1]@y1,pe[0]@y1) | 0 | 0.002 | 0.18 | 0.13 | 0.80 |
  | var(pe[1]@y1) | 0.5 | 0.518 | 0.08 | 0.10 | 0.95 |
  | var(res@y1) | 5 | 5.031 | 0.08 | 0.08 | 0.95 |
  | var(res@y2) | 4 | 4.028 | 0.22 | 0.22 | 0.95 |

  For `cov(pe[1]@y1,pe[0]@y1)` the mean standard error is below the empirical SD, and the
  coverage is 0.80 in 20 replicates.

* `integer64` id columns (`data.table::fread()`) are read as their exact digits. The
  latin1 and UTF-8 spellings of the same text are one id: the engine compares bytes, and
  every label now goes to it in UTF-8.
* Fits carry `pedigree_used`, `list(lines, founders, kept)`.
* `parentage_check(pedigree, genotypes, exclusion_threshold = 0.01, sex_markers =
  genotypes$sex_markers)`: for every animal genotyped together with its sire or dam, the
  markers called in both, the opposite homozygotes among them, their rate and `conflict`
  (rate above 1%, the `exclusion_threshold` of preGSf90; Hayes 2011). The denominator is
  the markers called in both, where preGSf90 divides by every SNP; nothing is imputed, so
  a low call rate is not refused here. The maternal-grandsire column of a `sire_mgs()`
  pedigree is not a parent pair. `sex_markers =` leaves the sex-chromosome markers out of
  the test, as preGSf90 does with `sex_chr`: PLINK codes a male X as homozygous, so a sire
  and his son are opposite homozygotes on about `2pq` of the X markers with no error at
  all. A repeated pedigree animal is refused, genotyped animals missing from the pedigree
  are a warning with the count, and none in it an error (zero pairs would read as no
  conflict).
* Every single-step fitter (`model()`, `model_mt()`, `gibbs()`, `h_inverse()` and through
  it `snp_blup()`, `estimate_gamma()`) runs the same parentage test on the genotypes it
  is given, without the columns named in `genotypes$sex_markers`, and warns once, naming
  the first pairs, when one exceeds 1%; the fit goes on with the data as given. The
  warning has class `br_parentage`. `options(uaitrem.verify_parentage = FALSE)` (or 0,
  preGSf90's `verify_parentage 0`) turns it off, for genotypes simulated without the
  pedigree.
* `genotypes$sex_markers`, an optional element (column indices, names or one logical per
  column of `m`): the sex-chromosome markers, left out of the parentage test of every
  fitter, of `parentage_check()` and of `qc_genotypes()`. They stay in G.
  `qc_genotypes(sex_markers = )` returns it as one logical per kept column, so its result
  goes to `genotypes =` without the fit warning again about the pairs that the X alone
  makes conflict. A logical of the wrong length (columns cut after the QC) is refused.
* `qc_genotypes()` filters animals too. `min_call_rate_animal = 0.90` (preGSf90's
  `callrateAnim`) drops the animals whose calls over the kept markers that vary fall below
  it, and `pedigree =` drops the offspring in conflict with a genotyped parent (the
  `verify_parentage 3` default of preGSf90). Call rates are filtered worst first, in
  passes: each pass looks at the call rate of every kept animal and marker, and the side
  holding the lowest rate loses its entries below the threshold and below the other side's
  lowest failing rate. A failed sample leaves without taking the markers, and a bad batch
  of markers without taking the animals (in the gate where 120 of 400 markers lost half of
  their calls every animal is below 0.90 over the 400, so a single PLINK-order pass,
  `--mind` before `--geno`, would drop all 60; the passes keep the 60 and drop the 120).
  MAF, HWE and the parentage test follow, and the passes repeat until nothing changes, so
  what `qc_genotypes()` keeps always passes the fitters' call-rate check. The matrix is
  read once: the counts are kept and updated as animals and markers leave, so
  heterogeneous missingness, which takes hundreds of passes, costs little more than one
  pass (see Performance). The input can be `list(ids, m)`; the result then carries `ids`
  and goes straight to `genotypes =`. New elements: `ids`, `rows`, `columns`, `n_animal`,
  `n_parentage`, `passes`, `removed` (one row per removed animal, with the reason, the
  rate and the pass), with a pedigree `parentage`, and with `sex_markers =` `sex_markers`.
  Removed animals are a message, a warning when more than half of the animals or markers
  went or when rows leave a bare matrix without rownames.
* `read_plink()` reads the .fam, the .bim and the id columns of a .raw as text:
  `id = "fid_iid"` names each sample `FID_IID` for files that repeat an IID across
  families, and a .bed comes back with `map` (the .bim: `chr`, `snp`, `cm`, `pos`, `a1`,
  `a2`) so the sex-chromosome markers can be found. A line with the wrong number of
  fields, and a .raw dosage that is not a number, are refused with the file, line and
  column.
* `options(uaitrem.min_call_rate_animal)` sets the call-rate line of every function that
  imputes a missing genotype by the mean (default 0.90; 0 turns the check off), and the
  default of `qc_genotypes(min_call_rate_animal =)` follows it.
* Missing-value codes and absurd values. Before the fit, every fitter that reads a
  numeric trait warns about two patterns, naming the trait, the count, the first row and
  the value:
  - a common missing-value code (-9, -99, -999, -9999, -99999, 99, 999, 9999, 99999 or
    -123456789) held by two records or more, outside the range of the others
    and at least 3 robust SDs from their median;

  On simulated data, ten records at -999, or three 1e4 above the rest, took h2 from 0.35
  to 1e-10 without a word. The same patterns in a `cov()` covariate give a warning that
  names it. `profile_theta()` gives these warnings once, from its outer fit, and none
  when `fit =` is given; the profile values are the same.
* `h2(fit, weight = , scale = )`. `weight` is the
  weight at which `var(residual)` of a fit with `weights` is read; by default, the
  harmonic mean of the weights of the records that entered the fit (see Breaking changes).
  `weight = 1` gives the ratio of the components as estimated, the number airemlf90
  reports, and `weight =` on a fit without weights is refused. `scale = TRUE` multiplies
  the pedigree components of a fit with metafounders by
  `var_scale = 1 + mean(diag(Gamma)) / 2 - mean(Gamma)` (Legarra et al. 2015), the usual
  base; the default `FALSE` reads them as estimated, on the base of Gamma. A value of
  `scale` other than `TRUE` or `FALSE` is refused.
* One scale rule in `h2()`, `h2_curve()`, the share column of
  `print()`/`summary()` and the h2 row of `summary()` of a `gibbs()` fit. Two warnings
  give the other reading:
  - metafounders with `var_scale` more than 10% from 1 when `scale` is not given, with the
    value on the usual base.

* The extractors read the state of the fit. `h2()`, `accuracy()`, `solutions()`
  and `rg()` warn on a fit that did not converge (the value is read at the last iterate).
  When the fit's SINGULAR verdict puts the numerator in a set of components the data do
  not separate, they return `NA` with a warning naming the components. A set that touches
  only the denominator, the residual side or another trait leaves the value, with a
  warning that says why it holds. A fit with no estimates gives `NA` and says so.
  `snp_blup()` names the conjugate-gradient stop. In a multi-trait fit `h2()` gives one
  warning for all traits.
* `accuracy()` warns, with the count and the first level, when a PEV exceeds the prior
  variance of its level (accuracy set to 0), as for a non-genotyped animal of a single
  step whose prior under H is not 1 + F.
* `residuals()` and `fitted()` for `model()`, `model_mt()` and `gibbs()` fits. The fit
  stores the conditional residual `y - X b - Z u`, the fitted values `X b + Z u` and `X
  b` and, in `model()`, the leverage `h_ii` of each record in the mixed model equations
  (Hodges 1998), read from the factor and selective inverse of the PEV step at no extra
  factorization; `model_mt()` compute theirs (and `resid_var`, the variance of each
  residual) with one more assembly and factorization at the final components, about one
  more REML iteration, and keep the fit with `NA` residuals if that pass fails. The
  variance of a residual is not `s2e`: `s2e (1 - h_ii) / w_i` averaged 0.35 s2e in
  contemporary groups of two and 0.69 s2e in groups of forty on a simulated animal model
  with h2 = 0.4, and moves with h2 and the layout, so `type = "standardized"` divides by
  it (a record with `1 - h_ii < 1e-8`, a group of one record, is `NA`). A `gibbs()` chain
  returns the posterior means (`"conditional"`, `"marginal"`; `fitted()` on the liability
  scale for the probit); with several chains they pool the chains that ran to the end,
  the same chains as `b` and `ebv`. Cox-Snell, martingale, deviance and randomized
  quantile), descriptive and not calibrated.
* `mcmc_diagnostics()`: convergence diagnostics for the chains of `gibbs()`, or for any
  draws: a matrix iterations x chains, an array iterations x chains x parameters, a list
  of chains, a draws object of the posterior package, or an `mcmc` / `mcmc.list` of coda.
  - One row per component, plus the heritability and the correlations computed draw by
    draw. Each row has the rank-normalized split R-hat, bulk and folded; bulk- and
    tail-ESS, and the Monte Carlo standard errors of the mean, sd and quantiles (Vehtari et
    al. 2021); the Pareto k-hat of the tails (Vehtari et al. 2024), with the MCSE of the
    mean left blank above 0.5; autocorrelations at lags 1, 10 and 50; the Geweke (1992) z
    and the Heidelberger-Welch (1983) test, combined over the chains, with the spectral
    density at zero estimated in the final window only (on log or ranks when the draws are
    heavy-tailed, see `start_scale`); and how far the chains disagree in location and in
    scale.
  - Across components, the multivariate ESS of Vats, Flegal & Jones (2019), by lugsail
    batch means over the chains, and the multivariate PSRF of Brooks & Gelman (1998) with
    their factor (1 + 1/M).
  - With two or more chains, the R-hat of every breeding value and fixed effect, flagged
    against a limit that grows with the number of levels and with the autocorrelation time
    of the slowest component. Levels of a group with two effects are named with the effect,
    "F016 (maternal)".
  - A verdict, `ok`, `run longer`, `burn-in too short`, `not mixed` or `chain stopped`,
    with every reason listed.
  - `gibbs()` runs it and keeps it in `$diagnostics`. `print()` and `summary()` show R-hat,
    bulk- and tail-ESS, MCSE and status per row, and end with the verdict and the notes.

  Gates (`tests/testthat/test-gibbs-diagnostics.R`):
  - R-hat, the bulk, tail, mean and basic ESS, and the MCSE of the mean, sd and quantiles
    equal posterior 1.7.0 to 1e-10 on 18 sets of draws, among them odd numbers of draws
    and 6 to 11 draws per chain, both for the functions and through the table. The Pareto
    k-hat matches to 1e-8 on the 10 sets with 21 or more draws per chain. With 29 or fewer
    draws in all, or with both tails constant, `khat` is `NA` where posterior returns
    `-Inf` (no fit in either); neither exceeds 0.5, so the verdict and the MCSE columns do
    not change.
  - The mPSRF equals `coda::gelman.diag()` when M = p, and the Brooks-Gelman formula when
    not.
  - The multivariate ESS of a VAR(1) is within 3% of its analytic value with the lugsail
    correction, and 5% high without it.

  Calibration (`validation/mcmc_diagnostics_calibration.R`):
  - 10 to 26 bulk draws per chain: never `not mixed`. Without the ESS floor, the
    disagreement rules would say so 12.5% to 61% of the time, and 37% to 97% counting the
    isolated chain.
  - One of 4 chains shifted 0.5 and 1 posterior sd: `not mixed` 83% and 100%.
  - A start 50 sd away with no burn-in: `burn-in too short` 100% with 4 chains and 86%
    with one, where R-hat passed 1.01 in 2% of the 4-chain runs.
  - 50 000 levels at an autocorrelation time of 5: the largest location R-hat passed 1.01
    in 95% of converged runs. The limit flagged none in 40 when given the true largest
    time + 2 as ceiling, and found a level with one chain shifted 0.5 posterior sd 90% of
    the time.
  - The ceiling a fit computes, the slowest row + 2, was checked on 860 simulated animals
    in three models. Largest level time over the slowest row's: 0.17 for the animal model,
    0.33 for direct-maternal (chains far from mixed), 0.20 for a two-coefficient reaction
    norm.

  Scale (`validation/mcmc_diagnostics_scale.R`, 1 thread): on a `gibbs()` of 4 chains x
  5000 kept draws over 20 000 levels of an iid effect (40 000 records), the diagnostics
  take 0.22 to 0.23 s (3 runs). Ten seeds gave `ok` with no level flagged. The levels
  there mix almost as independent draws, so this measures cost and false alarms, not the
  location rule nor an animal model.
* `posterior::as_draws_array()`, `posterior::as_draws()` and `coda::as.mcmc.list()` take a
  `gibbs()` fit, using the chains that ran to the end. This opens the rest of both
  packages, such as `posterior::rstar()` or `coda::raftery.diag()`. Both packages are
  suggested, not required.
* The fixed-effect design is sparse. `model()`, `model_mt()` and `gibbs()` built X dense:
  n x p on the whole table, again on the used rows, once more inside the Gram-Schmidt
  that found its rank (each orthogonalised column kept as an n-vector), and again in the
  final X. That came to about 5 n p doubles and O(n p^2) time in one thread: 4.95 GB for
  one `model(maxiter = 0)` on 209 646 records with 675 fixed columns. X is now kept by
  rows, one entry per record and fixed term. Its rank comes from a Cholesky of X'X in
  column order. A column that this step cannot decide reliably goes to an explicit
  residual (corrected semi-normal equations), and every decision there is certified: a
  residual at or below the threshold drops the column, and one with at most 1e-3 of
  itself left in the span keeps it. If a decision cannot be certified, the whole rank
  falls back to the previous Gram-Schmidt, run column by column without a dense X. That
  takes a column that nearly collinear columns before it span only with very large
  coefficients, such as x1, a near copy x2 = x1 + 1e-8 z, and then x1 + 1e8 (x2 - x1).
  Unscaled polynomial covariates and near copies on their own stay on the sparse route
  (tested to degree 8 on covariates with means from 1 000 to 2.46 million, and copies
  within 1e-4 to 1e-9). The criterion is unchanged: a column goes when its distance to
  the columns kept before it is below 1e-9 of its norm, so the last column of a dependent
  set goes. The rank never exceeds the number of records.

  `dropped_x` is the same as before on every fit of the test files other than this
  change's own (compared before the choice of the `model_mt()` reference level listed
  under Fixed, which changes some drops of a multi-trait fit). In that file the previous
  build kept more fixed columns than records in a few small designs (6 of 40 seeds of a
  10-record `model()` loop, and a `model_mt()` trait on 9 records); the new rank does
  not. On 1 324 adversarial random designs `dropped_x` is also the same, except in 23
  where the previous Gram-Schmidt kept a dependent column: more columns than records, or
  a column within 1e-9 of the span that its rounding let through. The new rank drops that
  column. Where the fallback runs, the result is the previous Gram-Schmidt's, rounding
  included, except that it never keeps more columns than records: it keeps the first
  columns the previous one kept, up to the number of records.

  Measured on simulated data, 3 isolated processes each, with `dropped_x` checked against
  the levels counted in the table (`validation/fixed_design_memory.R`):
  - 200 000 records, 701 candidate fixed columns (a 680-level contemporary group, two
    classes, two covariates and a nested one): `model(maxiter = 0)` peaks at 259-267 MB of
    working set (362-363 MB committed) in 2.6-4.2 s, depending on machine load. The fixed
    part is 145-153 MB of working set and 232 MB committed (three sets of 3 replicates).
  - 628 857 records: 617 MB of working set, 1.06 GB committed, 8.2 s. The fixed part is
    286 MB of working set and 675 MB committed.
  - 20 000 records: the previous build took 496 MB and 6.5-7.1 s for what now takes 33 MB
    and 0.3 s.

  The rank keeps a dense triangular factor of p(p + 1)/2 doubles for p kept columns, with
  O(p^3) work in one thread: about 0.4 GB and 4.5 minutes with 10 000 fixed columns (one
  run). The factor is sized before that work for the candidate columns with an entry on
  the used records (or the records, if fewer), dependent ones included. The fit stops
  with a message when that size passes 4 GB, or when the dense fallback would need more
  than 4 GB.

  Their rank uses the same engine route with the relative criterion of `qr()` (1e-7).
  Their products sum in the order of R's reference BLAS, so with that BLAS the fits are
  the same bit for bit. Weighted fits no longer rescan every random-effect incidence once
  per weighted record.

  The identity with the previous build was checked on Windows with R's reference BLAS,
  using `validation/fixed_design_identity.R`, before the changes to `model()`,
  `model_mt()` listed under Fixed (the -2logL from the residuals, the per-trait walk, the
  new start). Across the test files (all but the file of this change), the 3 085 values
  returned by 24 fitting and accessor functions are identical bit for bit. These include
  -2logL, theta, b, EBV, PEV, `dropped_x`, residuals and fitted values.
* Genotypes are read where they are. The three types give the same result bit for bit. A
  single step with an APY core of 2000 on 12 000 genotyped animals and 20 000 markers,
  median of 3: 52.2 s and 6.56 GB at the peak before, 29.8 s and 3.37 GB with the same
  double matrix, 27.8 s and 1.83 GB with a raw one.
* `pegs()` reads the genotypes where they are too: it held three double copies (the matrix
  imputed in R, the copy into the engine and a column-major one inside it), and now reads each
  marker from the R object at each update, straight from memory when the matrix is double and
  the marker complete. Same results for the three types; with 3000 animals, 5000 markers and 3
  traits, 100 passes take 5.8 s with a double matrix as before and 7.5 s with a raw one (median
  of 3), the price of converting each marker at every pass for an eighth of the memory. A
  marker with no observed genotype is now an error instead of a column of NaN.
* `read_blupf90_snp()`: the SNP_FILE of the BLUPF90 programs into a raw matrix, in two
  passes over the file so nothing else is held; `ids =` keeps a subset. `read_plink()` takes
  `storage = "integer"` or `"raw"`.
* `apy_core_select(method = "lanczos")`: the count of eigenvalues that explain 98% of G
  by stochastic Lanczos quadrature (Ubaru, Chen & Saad, 2017), from products with the
  genotype matrix only, never the Gram matrix nor its eigendecomposition. Each probe
  contributes the average of the Gauss rules of its last `ceiling(steps / 4)` leading
  tridiagonals: a single Gauss rule per probe is biased at the threshold, and more
  probes do not remove that bias, because the Ritz values near it fall in almost the
  same place for every probe. The count comes with a standard error over the probes,
  from the linearized count. On 1000 animals, eight independent draws of 20 sets of 30
  probes put the spread of the counts at 0.82 to 1.31 times the reported standard error
  at the gated levels (90% and 98%, mean about 1.0, `validation/apy_core_lanczos_se.R`).
  On 10 000 animals and 10 000 markers (`validation/apy_core_lanczos.R`, ten seeds) the
  exact count at 98% was 5953 and the estimate 5934 to 5975 (mean 5957.3, standard error
  about 12); with the same probes on the exact eigenvectors the quadrature error was
  +2.3 counts (0.04%), against +21.4 (0.36%) for a single Gauss rule. `"auto"` keeps the
  exact route up to 4000 animals or markers and takes Lanczos above.
* `simulate_breeding()` preallocates: it grew the haplotype matrices by `rbind` once per
  animal, O(n^2 m) copying, and 10 000 animals with 10 000 markers did not finish in hours;
  they take 26 s now, with the same population for the same seed.
* `sire_mgs(ped, dam = "dam")`: the MIXED pedigree, dams where they are recorded
  and maternal grandsires where they are not. A row with a known dam takes the
  sire-dam rules and only a row without one the grandsire path; the rules are per
  row, so both kinds share one A^-1. A grandsire that disagrees with the sire of
  the dam in the same row is an error, and a dam whose sire is unknown receives it
  from the grandsire column of her offspring, with a message. `pedigree()` reports
  the type `"mixed"` and which path each row took (`via_mgs`).
* `validation/sire_mgs_recovery.R`: the sire and maternal-grandsire model at scale. Twenty replicates of 1300 bulls over five generations and 36 000 daughter records: var(sire) 0.0762 (SE 0.0011) against the true 0.075, and the same files read as sire / dam give 0.0964, 29% too high, the error the declaration exists to prevent.
* `snp_blup()` takes any relationship structure `model()` takes. A group of dimension q
  with covariance `K0` gets q sets of marker effects, and since the joint covariance is
  `K0` times the scalar one, the precision is `K0^-1` times the scalar precision. `g`
  comes back as a marker x component matrix when there is more than one component, with
  the columns named as in theta. The applications of `A22^-1` in each iteration run in
  parallel.
* `solutions(fit, pedigree)`: one call, one data.frame with `id`, `ebv`, `se` and
  `acc`, sorted by breeding value. It is what an evaluation is for, and until now
  the caller assembled it by hand, joining the vector from `ebv()` to the vector
  from `accuracy()` by name.
* `apy_core_select()` and `apy_core = "auto"`: the APY core size is the number of
  eigenvalues of G that explain 98% of its trace (Pocrnic et al., 2016), the
  animals are drawn at random with a recorded seed, `size =` and `include =` (for
  proven sires) override the draw, and the object reports the cost of one
  factorization with that core against the exact one. The fit keeps the core it
  used in `fit$apy`.
* `sire_mgs()` and `pedigree(type = "sire_mgs")`: a pedigree of sires and
  maternal grandsires, declared, never guessed. The grandsire path weighs 1/4 and
  the Mendelian variance is `11/16 - F_s/4 - F_k/16` (Henderson, 1975, 1976);
  the declaration reaches every fitter and the single step. A third column named
  like a grandsire column (`mgs`, `mgsire`, `maternal_grandsire`) is refused
  unless declared.
* `setTimeLimit()` now stops the following between columns or panels, with the error
  "interrupted (Ctrl-C or a time limit)" and without leaving the factor allocated:
  - the sparse Cholesky factorization;
  - the minimum-degree ordering;
  - the tiled dense factorization of every fitter;
  The interruption comes back as a C++ exception, not as a jump over the engine's
  destructors. Ctrl-C goes through the same `R_CheckUserInterrupt()` check but was not
  tested interactively.
* `summary()` for `model_mt()`, `gibbs()` (posterior quantiles and the heritability
  sampled draw by draw), `snp_blup()`.
* A singular average-information matrix is reported by `model()` and `model_mt()`, naming
 the components the data do not separate. Their standard errors are NaN and the reported
 point depends on `start =`; only combinations of them are estimable.
* `fit$dense_block = c(dense = k, columns = n)` in `model()`, `model_mt()` and `gibbs()`:
  the size of the dense corner of the factor, the k of the k^3 each factorization pays.
  `model()` prints it with `verbose = TRUE`.
* In a single step, `accuracy()` divides a genotyped animal by its diagonal of G* also in
  `model_mt()`, which carried `1 + F` until now. `accuracy()` also takes on the
  log-hazard scale.
* The component table of `print()` and `summary()` has a `correlation` column for
  each covariance, and a `share` column that divides each variance or covariance
  by the phenotypic variance of its trait.

* Metafounders with genotypes: `model()`, `model_mt()`, `gibbs()` and `h_inverse()` build
  H(Gamma), with G centred at allele frequencies 0.5 and scaled by m/2 (G05), A22 taken
  from A(Gamma), and G* = (1 - w) G05 + w A22 without the affine adjustment to A22, which
  is the base correction Gamma already makes (Legarra et al., 2015; Garcia-Baccino et
  al., 2017). Every unknown parent must be a metafounder, or it is an error. `snp_blup()`
  centres its markers at 0.5 with scale m/2 and solves the same system (equal to
  `model()` at the same theta). The pair was refused before. Measured on 20 simulated
  replicates (1400 animals, two base populations with different allele frequencies and
  means, truncation selection, the last two generations genotyped, the last one predicted
  without phenotypes): `estimate_gamma()` within 0.011 of the true Gamma; against the
  plain single step, the bias of the genetic trend went from -0.141 (SE 0.020) to -0.081
  (0.018) genetic standard deviations and the accuracy from 0.282 to 0.294, the same with
  the estimated and the true Gamma.

* `estimate_gamma(method = "ml")`: the maximum likelihood of Gamma for a single
  metafounder (Legarra et al., 2024b). With one metafounder
  `A^gamma = (1 - gamma/2) A + gamma 11'`, so the Gaussian likelihood of the
  markers depends on gamma through three summaries of A22 and G, and is maximized
  on `[0, 2)` directly; the result carries the log-likelihood and a standard
  error from the curvature (optimistic, markers taken as independent). Gates:
  the identity against A(gamma) built by the pedigree code, the estimate equal to
  the maximum of the dense likelihood, and recovery of the base gamma.

* `sire(sire, mgs = "mgs")`: the sire and maternal-grandsire model, the record
  carrying 1 on its sire and 1/2 on its maternal grandsire in the same effect; an
  unknown grandsire leaves the sire only. Gate: -2logL and BLUP equal to the
  dense GLS with the incidence built by hand, over the A of the sire and
  maternal-grandsire pedigree.

* `pegs()`: the multivariate SNP-BLUP of Xavier and Habier (2022), ported from
  this project's validated Julia engine and giving the same numbers on the same
  data: the effects of one marker for every trait solved together by randomized
  Gauss-Seidel, variances by pseudo-expectation, bending (Hayes and Hill, 1981),
  and the `hcs` and `xfa` structures for many environments. Wide data or long
  data by `environment =`. It assumes the residuals of different traits are
  uncorrelated (the same trait in different environments); on traits of the same
  record the residual correlation leaks into `r_g`, and `model_mt()` is the route.

* `gibbs(chains =, cores =)`: independent chains, each from a seed drawn from R's
  generator, so `set.seed()` governs all of them and a PSOCK cluster (`cores =`,
  which also works on Windows) gives the same draws as running them in series.
  The result pools the chains (locations by the total variance over all draws)
  and reports `rhat`, the rank-normalized split R-hat of Vehtari et al. (2021),
  also exported as `rhat()`.
  Validated against Harville (1974): with flat priors the marginal posterior of
  the components is the REML likelihood, so a 2-D quadrature of the -2logL of
  `eval_internal()` gives the exact posterior. On 200 animals, four chains of
  60 000 iterations gave posterior means 0.6977 and 0.7757 against 0.6972 and
  0.7755 from the quadrature (z = 0.15 and 0.08 against the Monte Carlo error),
  the posterior standard deviations within 1 percent, R-hat 1.0007
  (`validation/gibbs_harville.R`).

* `gibbs(prior = "half_t")` and `half_t(df = 2, scale = 1)` put the hierarchical half-t
  prior of Huang and Wand (2013) on every random term; a scalar term is the case d = 1.
  Each standard deviation is half-t with 2 degrees of freedom and scale A_k, with density
  flat at zero, and each correlation of a covariance group is uniform on (-1, 1). The prior
  is proper, and both full conditionals are conjugate (C given U and a is inverse Wishart,
  a_k given C is inverse gamma with shape (nu + d)/2), so the chain stays exact Gibbs. The
  scale is read from the design: A_k = sd(y) / r_k (1 / r_k on the liability with
  `family = "probit"`), where r_k^2 is the mean squared row of the term's incidence times
  the mean diagonal of its covariance structure over the levels with records. The fit
  returns it in `prior_scale`, and the chain starts each variance at its share of
  var(y) / 2 divided by r_k^2. `half_t(df, scale)` changes the degrees of freedom and
  multiplies the scale, for every variance or by variance name. With `prior` unset, the
  half-t is now the prior of the scalar terms of a Gaussian trait, of a covariance group
  with 2d + 2 levels or fewer, and, with df = 4, of every term of a probit (see Breaking
  changes). Gates (`tests/testthat/test-gibbs-meia-t.R`,
  `validation/gibbs_half_t_exactness.R`):
  - with no data, the sampler's own conditionals recover the half-t marginal of each
    standard deviation and the (1 - r^2)^(nu/2 - 1) density of each correlation for d = 1,
    2, 3 and nu = 2, 4 (KS sqrt(n) D at most 1.72 in 400 tests of 1e5 draws, threshold 3),
    while the shape (nu + nl)/2 has no stationary distribution there (sqrt(n) D at least
    14.9 in 100 seeds);
  - with the same seed, y, 10 y and 0.1 y, a random-regression covariate multiplied by 10
    or 0.1, and K multiplied by 10 or 0.1 give the same draws rescaled (largest relative
    error 2.3e-13 over 200 and 1000 iterations);
  - the chain is the same with 1 and 4 threads, and in series or in parallel (pedigree, full
    single step and APY).
* `gibbs(prior_residual = )` sets the prior of the residual apart from `prior`:
  `"jeffreys"`, `"flat"`, `"uniform_sd"`, `c(df =, scale =)` or a half-t with A = sd(y)
  times the multiplier of `half_t(scale =)`. Left unset it is `"flat"` when `prior` is unset
  (see Breaking changes), `"uniform_sd"` when `prior` is a half-t, and `prior` itself when
  that names one of the other priors.
* `validation/gibbs_scalar_prior.R` measures the prior of scalar random terms and of the
  residual in `gibbs()`.
  - Scalar terms. One term alone, `random(id)` or `animal(id)`, on 2000 records with 5, 20
    or 200 levels and a planted variance of 0 or 0.05 (var(y) about 1), 30 data sets per
    cell, 2 chains of 50000, residual on `"jeffreys"`. A chain counts as stuck when it stays
    below 1e-6 var(y) for 1000 iterations or more in a row.
    - `prior = "jeffreys"`: every chain stuck at a planted 0 (60 of 60 in each of the six
      cells, draws down to 10^-260 var(y)), and 2 to 14 of 60 per cell at 0.05 in five of
      the six cells (35 of 360 chains, in 22 of 180 fits; one of them stops early, see
      Fixed, and before that fix it stopped the whole `gibbs()` call). The 95% intervals
      covered 0.05 in 0.77 to 0.97 of the data sets per cell, 160 of 180 in all.
    - `"half_t"`: no chain stuck at 0.05. At a planted 0, 8 of 360 chains met the rule, all
      with 200 levels (longest run 1814 iterations, smallest draw 10^-9.1 var(y)).
      Coverage at 0.05 was 0.93 to 1.00.
    - `"uniform_sd"` matched `"half_t"` on 20 and 200 levels. On 5 levels its posterior
      mean was larger (0.206 against 0.133 for `random(id)`, REML 0.050).
    - `"flat"`: no chain stuck, but on 5 levels the posterior mean was 1.34
      (`random(id)`) and 1.40 (`animal(id)`) against a REML of 0.050 and 0.047.
    - At a planted 0 with 200 levels both proper priors mix slowly: the median bulk ESS is
      116 to 128 from 2 x 49000 draws.
    - `prior_residual = "jeffreys"`: one data set kept 9.3% of a chain's draws with
      var(residual) below 0.1 (REML 0.597), and its bulk ESS was 255 after 200000
      iterations.
    - `"uniform_sd"`: the same data set kept 2.5%, and every bulk ESS of var(residual) was
      427 or more.
    - `"half_t"`: 4.5%, and every bulk ESS was 422 or more.
    - Coverage of var(residual) was 0.93 in all three.
  - The half-t at the size of a herd (mode padrao, with `"uniform_sd"` on the residual and
    no ASIS step, both named in the mode): the population of `simulate_breeding()` with 100
    founders and 4 generations of 500 (2100 animals), 30 data sets per design, 2 chains
    doubled from 20000 until the smallest bulk ESS reached 400.
    - One record per animal, `y ~ cg + animal(id)`: posterior means 0.293 and 0.702 against
      REML 0.291 and 0.701 (planted 0.3 and 0.7), coverage 0.90 and 0.93, a bulk ESS of 400
      and an R-hat under 1.01 in 30 of 30 at 2 chains of 20000.
    - Two records per animal, `y ~ cg + animal(id) + pe(id)` (4200 records): 0.303 / 0.142 /
      0.549 against 0.300 / 0.144 / 0.548 (planted 0.3 / 0.15 / 0.55), coverage 1.00 / 1.00 /
      0.97, the bulk ESS in 30 of 30 and the R-hat in 28.
    - Against `prior = "jeffreys"` on the same chains the posterior means differed by at most
      0.0013 (SE 0.0004).
  - Mode portao gives the margins of the default arm of the REML gate of
    `tests/testthat/test-gibbs.R`: over 100 chain seeds and 20 data seeds the posterior mean
    under the default (the ASIS step on the scalar term, `"flat"` on the residual) sat at
    most 0.021 from REML in `var(animal)` and 0.032 in `var(residual)` (threshold 0.2), the
    same at 1 and 4 threads.
  - Mode asis decides the ASIS step on scalar terms: the twelve cells of the scalar
    measurement above, `half_t()` on the term and `"uniform_sd"` on the residual, with and
    without `asis = TRUE` from the same seeds, 2 chains of 50000. Rule written before the
    run, per iteration: the lower 95% bound of the
    paired ratio of the smallest bulk ESS above 1 at a planted 0.05 (the six cells pooled),
    and at 0 no upper bound below 1, pooled or in any cell. It passed: 2.12 (1.90 to 2.36)
    at 0.05 and 3.47 (2.97 to 4.04) at 0, from 1.09 (1.06 to 1.12; `random()`, 5 levels,
    0.05) to 13.61 (11.67 to 15.86; `animal()`, 200 levels, 0), and the ESS of
    `var(residual)` did not fall in any cell (0.99 to 2.62). With 200 levels and nothing
    planted, a bulk ESS of 400 and an R-hat under 1.01 on both components came in 0 and 1
    of 30 data sets without the step and in 29 and 30 with it. Chains held at zero: 8 of
    360 without it and 3 with it at a planted 0, none at 0.05; coverage of 0.05 171 of 180
    in both; the paired difference of the posterior means of the term -0.018 (SE 0.013) to
    +0.012 (0.012) posterior SD. The step costs O(records), as much as the rest of an
    iteration of this one-term model (1.36 to 1.96 times the CPU, 3 interleaved runs per
    cell), so per second the same rule fails: 1.44 (1.28 to 1.61) and 2.26 (1.93 to 2.64)
    pooled, but 0.82 (0.78 to 0.86) and 0.93 (0.88 to 0.98) with 5 levels at a planted 0,
    0.77 and 0.69 at 0.05, and 0.94 (0.91 to 0.97) on `random()` with 20 levels at 0.05.
  - Mode misto runs the default with `"uniform_sd"` on the residual (`"flat"` and the ASIS
    step on the group, the half-t and the step on `cg`) on
    `random(cg) + animal(id, group = "g") + maternal(dam, group = "g")` (2400 records, 2650
    animals), with and without the step (`asis = FALSE`): the chains are `identical()` at 1
    and 4 threads, in series and in parallel, in both arms (2 data sets), and the step on the
    three terms cost 1.02 and 1.10 times the CPU in two runs of 3 interleaved repetitions on
    data set 1, on a machine shared with other jobs.
* `gibbs(prior_fixed = )` puts a normal prior with mean 0 and that standard deviation on
  every fixed-effect coefficient of the liability as the engine codes it, the intercept and
  the contrasts against the reference class, the last in the order of the data
  (`family = "probit"`; `Inf` is flat). It is now the default with a half-t prior and a
  component estimated (see Breaking changes). With flat fixed effects and a half-t prior,
  `gibbs()` warns that the posterior can be improper when the model has at least as many
  fixed-effect columns as the df (it is when the effects can reproduce the category of
  every record, as with one record per animal; the warning counts columns only). On the
  design of Odegard et al. (2010) (2000 animals with one binary record each, 80 fixed
  classes; `validation/gibbs_probit_odegard.R`, 10 data sets, 4 chains of 25000), every
  flat-b prior, `half_t(df = 4)` included, drifted to `var(animal)` of about 1500 in 8 of 10
  data sets (drift +830 to +915 between chain halves, SE about 150). `prior_fixed = 1` cut
  the drift to +11 (SE 4) under `half_t(df = 4)` and to +2 (SE 5) with `px` as well, the
  default, but under the default 4 of 10 data sets still sat at medians of 15 to 62 after
  25000 iterations, and the 95% interval covered the planted 0.25 in 3 of 10 (without `px`,
  3 of 10 at 11 to 49, coverage 3 of 10; under `half_t()`, 4 of 10 at 16 to 70). With and
  without `px` the target is the same, yet 3 data sets sat above 10 under one arm and not
  the other: 25000 iterations describe these chains, not that posterior. `prior_fixed = 3`
  (without `px`) let it back up to medians of 262 to 408. The script now prints the median
  per data set and arm. Why 4 data sets still drift under the default is open: their full
  sibs agree in category as expected at `h2 = 0.2`, but the two that agree most sat high
  with and without `px`, which fits slow mixing and a tail of the posterior alike.
* `gibbs(px = TRUE)` adds the Haar PX-DA move of Wang and Roy (2018, Algorithm 2) after
  the liability draw. Gate: chains with and without it agree within Monte Carlo error on
  two variance components, the fixed effects and 70 levels (16 chains of 3000 each way,
  largest |z| 3.44 over 20 data sets; thresholds 4.5 and 5). Both reproduce a
  fixed-effect posterior computed by quadrature (|z| up to 2.35, standard deviation
  within 1.05%). The chain is identical with 1 and 4 threads. The lower bound stopped at
  1.00 on one of the three, short of the rule written before the runs; it is on by
  default with `family = "probit"` by a decision that overrides that rule (see Breaking
  changes).
* `gibbs(start = )` sets where each chain starts: one vector, or a list with one per chain.
* `gibbs(asis = TRUE)` adds the ancillarity-sufficiency interweaving step of Yu and Meng
  (2011) to every covariance group, a scalar random term being a group of one, for
  `family = "gaussian"` and `"probit"`, under `"flat"`, `"uniform_sd"`, `c(df =, scale
  =)` and `"half_t"`. It is refused on any estimated group or scalar term under
  `"jeffreys"`: there the target of the step, the conditional of `chol(C)` given the
  non-centred effects, grows like `prod_i L_ii^-i` near a zero diagonal entry, where the
  likelihood stays positive, and cannot be integrated. Run anyway
  (`validation/gibbs_asis.R`, mode jeffreys: direct-maternal group plus a scalar term on
  960 records, `"jeffreys"` on both, 16 chains of 3000, 5 data sets), the step held 12 to
  16 of the 16 chains at a singular correlation in every data set, against none without
  it, and accepted 0.7% to 1.0% of its proposals on the group; under `"flat"` no chain
  was held either way (acceptance 0.79 to 0.86). After the usual draw of a group's
  covariance `C`, the effects are rewritten as `U L^-T` with `L = chol(C)`, and `C` is
  drawn again by Metropolis-Hastings. The proposal is the Gaussian posterior of the lower
  triangle of `L` in the regression of the partial residual on the columns `Z_j U~_m`.
  The acceptance is the ratio of the prior density of `L L'`, given the half-t
  auxiliaries, times the Jacobian `2^d prod L_ii^(d - i + 1)`. The fit reports the
  acceptance in `asis_acceptance` and `chain_asis_acceptance`.
  - Gate: `validation/gibbs_asis.R` and `tests/testthat/test-gibbs-asis.R`. With and without
    the step, 16 chains each way, the chains agree within Monte Carlo error on every
    component, fixed effect and level, over 20 data sets of each of four small models. The
    largest |z| was 3.10 on components, 3.58 on fixed effects and 5.10 on levels; the test
    thresholds are 4.5 and 5.5. On the 2 x 2 groups, the step without the Jacobian fails on
    the components (|z| at least 5.84), and the step without the Metropolis-Hastings
    correction fails on the components or the levels (|z| at least 4.25 or 5.57). The same
    gate on the Gaussian default (`"flat"` on the direct-maternal group, the half-t on a
    scalar term, `"flat"` on the residual; mode exato, model mat_padrao, 20 data sets) gave
    at most 2.06, 3.26 and 3.91, and rejected the step without the Jacobian (|z| at least
    10.40 on the components) and without the correction (at least 7.32). The chain is
    identical with 1 and 4 threads and stays equivariant under `"half_t"`.

* `br_threads()`: OpenMP where a genomic evaluation spends its k^3. The dense
  tail of every sparse Cholesky (the genotyped block after the ordering) is
  factored in tiles of 64 x 64 inside its own column storage (Buttari et al.,
  2009); the inverse of that tail inside the selected inverse, the dense inverses
  of order 256 or more (G* and A22 in the single step) and the columns of A22
  run in parallel too. Every number has one owner thread and a fixed summation
  order, so any thread count gives the same bits: a test fits a genomic model
  with 1 and 3 threads and requires identical -2logL and solutions, and the
  suite passes with 1 and 4. The default is 1 thread; `options(uaitrem.threads
  = )` or `UAITREM_THREADS` set it at load, and `OMP_THREAD_LIMIT` caps it.
  `br_threads(lapack = TRUE)` hands the dense inverses and products to R's
  BLAS and LAPACK, for an optimized multithreaded BLAS.

## Breaking changes

* The defaults of `gibbs()` (`prior`, `prior_residual`, `prior_fixed`, `px` and `asis`
  left unset) changed: each part was decided against the measurements of the New section.
  It was not run on the direct-maternal and `rn()` groups, on the animal models of the
  size of a herd or on scalar terms with 200 levels, which ran with `"uniform_sd"` on the
  residual (the herd-size models also without the ASIS step), and the probit default was
  not run on the many small classes and the low incidence of
  `validation/gibbs_probit_default.R`, which ran without `px`. For the residual, the
  probit and the two sampler steps a rule was written before the runs (for the residual,
  after a pilot; `px` is on against its rule); for scalar terms and covariance groups the
  measurements were asked for, and the reading below was set when deciding.
    - direct-maternal at -0.3 and -0.8 and `rn()` with d = 2 and 3
      (`validation/gibbs_group_prior.R`, mode padrao, 30 data sets each, 2 chains of 16000
      and of 8000), with `"uniform_sd"` on the residual, not run with `"flat"` there: no
      chain of 240 held, coverage 0.900 to 1.000, the correlation of -0.8 +0.038 posterior
      SD (SE 0.009) from REML toward 0, and a bulk ESS of 400 and an R-hat under 1.01 on
      every component in every data set.

    A group of d components with 2d + 2 levels or fewer takes `half_t()` instead, with a
    message, and the fit names it in `prior_fallback`: under `"flat"` its conditional, an
    inverse Wishart with nl - d - 1 degrees of freedom, has no finite mean there (and is
    improper with 2d levels or fewer), so the posterior mean the fit reports would not
    exist. The message gives the group's d and number of levels. With 2d + 3 or 2d + 4
    levels the group keeps `"flat"`: its posterior mean exists, but its posterior SD does
    not (the conditional has a finite variance only above d + 3 degrees of freedom), so
    the `sd` of the fit and the Monte Carlo SE of those components estimate nothing
    there. `prior = "flat"` named keeps its old behaviour (an error with 2d levels or
    fewer).
  - The residual of a Gaussian trait gets `"flat"`, the prior the covariance-group
    measurements had on it; it had `"jeffreys"`. Three priors were measured there, each
    against a rule written before its runs (the one for `"flat"` after a pilot, below):
    - `"jeffreys"` leaves mass near a zero residual with one record per animal: with
      `"half_t"` on the group, one data set kept 9.3% of a chain's draws below 0.1 against
      a REML of 0.597 (`validation/gibbs_scalar_prior.R`, mode residuo, 30 data sets);

    That last clause did not say whose coverage; it was read after the run as that of the
    scalar term (113 against 113 of 120). The residual's was 109 against 110 of 120 at a
    planted 0.05: read that way the rule fails by one data set, and `"uniform_sd"` would
    have stayed. `prior_residual = "uniform_sd"` gives that residual with the rest of the
    default; with `prior = "half_t"` named the residual still gets `"uniform_sd"`.
  - `asis` and `px` default to `NULL`. With `prior` unset the ASIS step runs on every
    estimated term of a Gaussian trait, covariance groups and scalar terms, and on every
    estimated covariance group of a probit; with `family = "probit"` and a component
    estimated the PX-DA move runs too. Rule written for both before the runs: keep the step
    if the lower 95% bound of its paired ratio of the smallest bulk ESS passes 1. The ASIS
    step passed it; PX-DA did not on one design and is on by a decision that overrides it.
    - ASIS on a scalar term (`validation/gibbs_scalar_prior.R`, mode asis: `random()` or
      `animal()` with 5, 20 and 200 levels on 2000 records, planted 0 and 0.05 of var(y),
      30 data sets per case): 2.12 (1.90 to 2.36) times per iteration at 0.05 and 3.47
      (2.97 to 4.04) at 0, above 1 in every one of the twelve cases (1.09 to 13.6), which
      passed the rule; with 200 levels and nothing planted, the data sets with a bulk ESS of
      400 and an R-hat under 1.01 went from 0 and 1 of 30 to 29 and 30. Per second the same
      rule fails in that one-term model, where the step costs as much as the rest of an
      iteration (1.36 to 1.96 times the CPU): 0.82 and 0.93 times with 5 levels at a
      planted 0, 0.77 and 0.69 at 0.05, and 0.94 on `random()` with 20 levels at 0.05. In
      `random(cg) + animal() + maternal()` on 2400 records the step on the group and on `cg`
      together cost 1.02 and 1.10 times the CPU in two runs on a shared machine (mode
      misto). It was not measured on probit scalar terms, which keep it off.

  So every default Gaussian chain with an estimated component changes (through the
  residual, each scalar term and each covariance group, and the ASIS step), and every
  default probit chain with an estimated component. With `family = "probit"` those three
  priors now also need `allow_improper = TRUE`. With `theta_fixed`, or with no component
  estimated, the probit fixed effects stay flat and `px` and `asis` stay off. Tests and
  validations of the sampler's mechanics name the prior they were calibrated with;
  `tests/testthat/test-gibbs-padrao.R` checks each part of the default against the same
  configuration named. The gates that read a default chain
  (`tests/testthat/test-gibbs.R` and `test-gibbs-probit.R`) take their thresholds
  from the default itself, measured over 60 to 120 chain seeds and 15 to 60 data seeds
  with the same result at 1 and 4 threads (the modes portao and fronteira of the
  validation scripts, and margem of `validation/gibbs_probit_default.R`).
* `fit$prior` of `gibbs()` is the prior of the covariance groups as the chain used it, never
  `NULL`; `fit$prior_scalar` is new (the prior of the scalar terms); `fit$prior_residual` and
  `fit$prior_fixed` hold the values used (`NULL` for the probit residual; `Inf` for flat
  fixed effects) instead of the arguments as given. `prior_fixed = Inf` asks for flat fixed
  effects, and a Gaussian fit takes only that. `fit$prior_fallback` (new) names the
  covariance groups the default moved from `"flat"` to `half_t()`; `fit$px` is the PX-DA
  move as run; `fit$asis` is `NULL` when the step did not run, and under the default of a
  probit its rows for scalar terms show no attempt (`asis_acceptance` `NA`).

* The ids of the data are matched to the pedigree by one rule in every fitter (`model()`,
  `model_mt()`, `gibbs()`, `snp_blup()`, the `eval_internal` family). `accuracy()`,
  `solutions()`, `snp_effects()`, `predict()` redo it from the pedigree given. Fits of
  the same data change:
    - the engine fitters dropped the records of an absent animal, sire or dam in silence,
      and only `n_used` fell (110 of 120 records with 10 animals absent; theta 0.249/0.614
      against 0.331/0.637 with the complete pedigree on the survey case);
  - The pedigree is pruned to the animals of the data (every row, the trait observed or
    not), the genotyped and their ancestors. This is the informative set of renumf90
    (`PED_DEPTH` deep enough to reach every ancestor), and a message counts the lines left
    out. The likelihood and every breeding value kept are unchanged. But `ebv()`,
    `solutions()`, `accuracy()` and `predict()` no longer have an animal with no record,
    no genotype and no recorded descendant; its EBV was the parent average. Give such an
    animal a row with the trait missing to keep it. A pedigree of 50 220 lines for 220
    records, 50 000 of them random matings unrelated to the data, ran past 120 s at
    12aca18 and now fits in 0.02 s.
  - Data that used to fit now stop, each with a message that names the column, the count
    and the first culprit:
    - an id written differently from its pedigree line (spaces, tabs, carriage returns or
      no-break spaces at the ends, upper against lower case, leading zeros, a trailing .0,
      decomposed accents, latin1 bytes read as UTF-8). It stops with both spellings and
      the difference. Before, those records left in silence (51 of 120 entered with
      leading zeros);
    - in `pegs()`, a data id written differently from its `genotypes$ids` line. It stops
      with both spellings; before, its record left as one without genotype (3 records with
      a trailing space were dropped and only counted);
    - a parent written differently from its own line, and a 0 written another way in the
      pedigree (0 followed by a carriage return, "00" with no line). Both stop with the
      two spellings;
    - a pedigree that is not a data.frame, or has fewer than three columns;
  - `accuracy()` accepts the pedigree the fit was given with lines that are no ancestor of
    any level; they are pruned as in the fit. A pedigree with ancestors added or levels
    missing is still refused.
* Genotypes that used to run now stop, with an error that names the rows and ids:
  - in every function that takes `genotypes =` (`model()`, `model_mt()`,
    `gibbs()`, `h_inverse()`, `snp_blup()`,
    `pegs()`, `g_matrix()`, `apy_core_select()`, `eval_internal()`,
    `estimate_gamma()`, `snp_effects()`,
    `parentage_check()`, `qc_genotypes()` with a list, through
    `model()`): an id repeated in `genotypes$ids` (the single step stopped with "the
    matrix is not positive-definite"; `snp_blup()`, `pegs()`, `g_matrix()` and
    `apy_core_select()` ran with one animal too many); a missing id (`NA`) in
    `genotypes$ids` (`g_matrix()` ran with an animal named "NA"); rownames of
    `genotypes$m` that disagree with `genotypes$ids` (the genotypes went to the animals at
    those positions in `ids`, without a word); genotypes with no rows (a
    `list(ids = character(0), m = ...)` gave a fit on the pedigree alone, silently) or no
    columns;
  - an animal with a call rate below 0.90 over the markers that vary, wherever a missing
    genotype is imputed by the marker mean: `model()`, `model_mt()`,
    `gibbs()`, `h_inverse()`, `snp_blup()`,
    `pegs()`, `g_matrix()`, `apy_core_select()`, `eval_internal()`,
    `snp_effects()`, and
    `apy_inverse()` with a double or integer matrix.
    `parentage_check()` imputes nothing and does not refuse it,
    `estimate_gamma()` refuses any missing genotype with its own message, and
    `qc_genotypes()` drops those animals instead.
    `options(uaitrem.min_call_rate_animal = 0)` restores the old behaviour.
* The single-step fitters warn on genotyped parent-offspring pairs with more than 1%
  opposite homozygotes (see New). A test or a script that simulates genotypes
  independently of the pedigree needs `options(uaitrem.verify_parentage = FALSE)`; any
  value of that option other than TRUE, FALSE, 1 or 0 is an error. Sex-chromosome markers
  make a sire and son conflict without any error: name them in `genotypes$sex_markers` or
  leave them out.
* `qc_genotypes()` now drops animals by default (`min_call_rate_animal = 0.90`); pass
  `min_call_rate_animal = NULL` for the markers alone. The marker call rate, MAF and HWE
  are counted over the kept animals, so `n_call`, `n_maf` and `n_hwe` can differ from
  before on the same data. A filter that leaves no animal or no marker is an error that
  says what removed them (it returned an empty matrix: `qc_genotypes(matrix(NA, 3, 4))`
  gave `n_call = 4`). The HWE chi-square is the same test, now computed from the genotype
  counts of the kept animals.
* `read_plink()` returns ids as written in the file: an IID `0001` was `1`, `T`/`F` were
  `1`/`0`. An IID repeated in the file is refused (pass `id = "fid_iid"` when the families
  tell the samples apart). A .bed result has the extra element `map`.
* A fit whose time axis cannot represent the serial correlation the data carry now stops,
  or warns, instead of returning a wrong fit with `converged = TRUE`. The data are fitted
  again with the records off the grid moved to its nearest free point (or with the order
  of the records within subject as time), and with no serial correlation (`model()`, or
  `model_mt()` under `cbind()`) on the records of the fit. A negative correlation the
  axis cannot represent, with a likelihood-ratio gain above 10.83 (0.1%, 1 df), stops the
  fit with an error that says what to do and names the gaps off the grid (subject, both
  times, gap length); from 3.84 to 10.83 a warning says the same. A numeric time is
  printed there with up to 15 significant digits and no scientific notation, so a time in
  epoch milliseconds (as `data.table::fread()` reads a timestamp, double or `integer64`)
  names the records exactly: "t = 1704153600000 to 1704153600500 (gap 500)". A
  correlation the axis can represent restarts the walk from the refit, and a gain still
  above 10.83 stops the fit ("the spacing of the records contradicts their serial
  correlation"). Data that fitted at 12aca18 and now stop (cases):
* With `missing_code`, a covariate cell equal to the code is refused on every row, like
  NA, when it lies more than 10 robust SDs from the other cells, and warned about
  otherwise (a 0/1 covariate under `missing_code = 0`). A -999 age had been fitted as an
  age: the regression fell from 0.0116 to 0.0017.
* A factor or text where a number is read is refused with the way to convert it.
* In the same fitters, a factor whose labels are numbers out of numeric order is refused,
  and so is an `ordered()` factor whose levels are in text order. This is what `factor()`
  and `ordered()` make of a numeric text column without `levels =` ("1" < "10" < "2"). A
  deliberate order such as `factor(x, levels = 3:1)` must now be declared with
  `ordered()`. Labels that write one number twice ("1" and "1.0", " 2" and "2") are
  refused, as factor levels and as text.
* An unordered factor trait whose levels are in alphabetical order runs in `gibbs(family
  = "probit")` with a warning. This is what `factor()` and `read.csv(stringsAsFactors =
  TRUE)` give a text column. The warning shows the order used, the upper category of a
  binary trait and the record count of each category. `ordered()` declares the order and
  silences it.
* A time stored as a factor or as text is refused, naming the column and the labels that
  do not read as numbers. Before, a factor entered by its codes 1, 2, ... without a word,
  and a factor entry gave a misleading "entry >= stop" error. Numeric text ("40"), which
  used to be read, is refused too, as `model()` refuses a text trait.
* `h2()` reads the marker of the term:
  - in a `sire()` model it returns `4 var(sire) / Vp` (it returned `var(sire) / Vp`, a
    quarter of the heritability);
  - in `sire(sire, mgs = )` the phenotypic variance also holds `f var(sire) / 4`, with `f`
    the fraction of the records with a known maternal grandsire (per trait in `model_mt()`
   );
  - on a `pe()` or `random()` group it stops: `var / Vp` there is a variance ratio, which
    the share column gives.

  The h2 row of `summary()` of a `gibbs()` fit follows the same rule (it was the first
  `var()` over the sum).
* The default group of `ebv()`, `solutions()`, `accuracy()`, `h2()`, `h2_curve()`,
  `plot()` and of the `gibbs()` h2 row is the group of the direct genetic term. It was
  the first group, so `random(lot)` written first made `ebv()` return the lots and `h2()`
  return `var(random) / Vp`. Inside a named group `h2()` reads the direct term in any
  formula order.
* On a fit with `weights`, `h2()`, `h2_curve()` and the share column read
  `var(residual)` at the harmonic mean weight of the records that entered the fit, the
  mean residual variance of those records. Weights w and 10 w now give the same h2, as
  they give the same breeding values. `weight = 1` gives the old number, which is also the
  one airemlf90 reports.
* Designs that ran now stop, before any iteration, in every fitter that estimates
  components. Before, the component walked to 1e8 or more, or to zero, with `converged =
  FALSE`; on a simulated one-sire-per-group design, var(sire) ended at 4.3e8, 1.2e13,
  1.1e14 or 1e-10 depending on the seed and the start, see
  `validation/structure_cases.R`. Or the chain of `gibbs()` walked its prior (var(animal)
  0.73 and var(residual) 0.011 against 0.4 and 0.6 with `cg = id`). The designs:
  - a random term whose levels are a function of a fixed class effect: one sire per
    contemporary group; a fixed litter with `maternal(dam)` or `pe(dam)` whenever each
    litter has one dam, whatever the number of litters per dam; `cg = id`. In `model_mt()`
    this is checked on the rows of each trait too;
  - a random term with a single level;
  - a term whose records all carry the unknown parent code 0;

  `maxiter = 0` with `n_em > 0` estimates (the EM moves the start) and stops on these
  designs too.
* `maxiter = 0` without `start =` is no longer read as given components: the extractors
  warn that the fit did not converge.
* `accuracy()` on a fit with no estimates returns `NA` with a warning (it stopped with "no
  component 'var(animal)'").
* Data that used to fit now stop. Every fitter refuses a missing label in a column that
  gives a level, and a missing value in a covariate. A missing label is NA, NaN, Inf,
  -Inf, blank text or the text "NaN"; a missing value is NA, NaN, Inf or -Inf. Rows whose
  trait is missing are checked too, and three cases that fitted at commit f2c0d69 (n_used
  47 of 48) now stop:
  - a blank label on a row with a missing trait stops every fitter, `model()` included.
    Before, the every-row rule held only in the engine fitters, and only for numeric `NA`,
    `NaN` and `Inf`;
* `weights` given as text or a factor are refused.
* Error texts changed. Code that matches them needs the new forms:
  - "Inf or -Inf in the column 'g'" (was "Inf");
  - "weights: NA or Inf in the weight column 'w' in 2 row(s), the first at row 4".
* `snp_effects(blend =)` now defaults to `NULL`, and a value different from the fit's is
  refused. A fit made before this version has no `fit$genomic` and must be refitted.

## Fixed

* The advice on repeated records ("the model has no pe(id)") did not see a `pe(id, base =)`,
  the permanent effect of a random regression. So a model with `animal()` and `pe(base =)`
  got the advice. A `pe()` with a base now counts as the permanent effect.

* `model()` and `model_mt()`: when a covariance group stops at a wall (a correlation at
  +/-1), the standard errors are conditional on that wall. The fit holds the group there and
  optimizes the other components conditional on it. So `vcov` is now
  `2 J_F (J_F' AI J_F)^-1 J_F'`, with `J` the derivative of `theta` in the coordinates of the
  step and `F` the coordinates off the wall. It has one rank less. Before, it was the full
  `2 AI^-1`. A fit without a wall keeps its numbers to the last bit. The simulation was a
  random regression of rank 1, with r = 1 in the fit. There the standard error of the
  permanent-environment variance was 0.0518. Now it is 0.0486, against 0.0489 from the
  curvature of the profile -2logL (`profile_theta()`). Gates in
  `tests/testthat/test-ep-parede-bloco.R`.

* `model()` stops when the step makes no progress for 25 iterations, as `model_mt()` and
  `model_ar1()` already did. No progress means no gain above 1e-8 of the -2logL and no EM
  step that moves. Before, the walk went on to `maxiter` at the same point. The verdict
  follows the Newton decrement, as at `maxiter`. Inside the tolerance, it is the
  convergence of a flat ridge. Above it, the message says that the point is a stall.

* With `br_threads()` above 1, a small fit no longer opens a parallel region at each
  evaluation. The trace tr(K^-1 Z) of a block now opens one from 1024 levels. The dense
  block at the end of the selected inversion now opens one from 193 columns. The fit is the
  same to the last bit with any number of threads, on both sides of these thresholds. Gate in
  `tests/testthat/test-threads-limiar.R`.

* `model_mt()` refuses two traits that are exactly collinear on their common records (the
  fit gave only boundary messages). Gate in `tests/testthat/test-estrutura-media5.R`.

* A base of `rn()` or `pe(base =)` whose columns differ in scale by more than 100 times (t
  and t^2 without scaling) gets a warning with the cause and the way out (`legendre()`, or
  scaled columns): the fit stopped at `maxiter` with "Raise maxiter=". Gate in
  `tests/testthat/test-estrutura-media4.R`.

* A singular `gamma` (`gamma = 0`, or two metafounders with identical rows) is refused by
  the fitters with the cause and the way out (a small diagonal for the limit of unknown
  parent groups, one label for one population). The help said that the fitters took it;
  they stopped with "the starting theta is INADMISSIBLE" or "the pedigree A^-1 is not
  positive-definite". `a_inverse()` still takes it, through the pseudo-inverse. Also: a
  record with a metafounder label is refused (it entered as an animal); several parents
  cited without a line are counted and listed (the message named the first only);
  `qc_phenotypes()` says why for a text trait or a single record (it stopped with generic R
  errors). Gates in `tests/testthat/test-estrutura-media3.R`.

* More cases of the structure of the data stop with the cause and what to do, in place of
  a generic error or a silent pass. The notation of other packages (`cg:sex`, `cg*sex`,
  `factor()`, `I()`, `poly()`, `(1 | g)`, a transformed trait such as `log(y)`) gets a
  message that says what to write here. `rg()` with a trait that is not in the fit, or the
  same trait twice, says so (it stopped with "subscript out of bounds"). `apy_core =` with a
  repeated id (it stopped with "not positive-definite"), or with positions in place of ids,
  is refused. A single step with fewer than 3 genotyped animals is refused: the allele
  frequencies of G come from them. A trait of `cbind()` in `model_mt()` with no record is
  refused by name. Gates in `tests/testthat/test-estrutura-media2.R`.

* Six cases of the structure of the data that passed without a word now stop with the
  cause, or warn where the data can be legitimate. A trait name with a character that the
  component names use (`@ | , ( ) [ ]`) is refused in `model_mt()` (`h2()` gave `Inf` for
  `p@2`). `gibbs()` warns when `burnin >= n_iter` (the chain keeps no draw). A base of `rn()` or
  `pe(base =)` with linearly dependent columns is refused when the fit estimates. Two
  genotyped animals with the same genotype at every marker get a warning with the pairs (a
  duplicated sample, twins or clones). A marker with the same genotype in every animal (all
  1, a failed call) leaves G like a monomorphic one, and a panel of only such markers is
  refused. `gamma =` without `metafounders =` is refused (it was ignored). Gates in
  `tests/testthat/test-estrutura-media.R`.

* `model()` and `model_mt()`: a design whose fixed effects keep as many columns as the
  records of a trait (no residual degree of freedom for REML) is refused when the
  components are estimated, with the cause. Before, the fit stopped with "did NOT converge"
  and a Newton decrement of 1e14, and the -2logL was not defined. With the components given
  (`maxiter = 0`) the solutions are computed, and `message` says it. In `model_mt()` the
  count is per trait. Gates in `tests/testthat/test-gl-residuais.R`.

* `model()`: a `start=` with the residual variance below its floor (1e-10 `var(y)`) stopped
  with "the starting theta is INADMISSIBLE", which was false (the value is positive; the mixed
  model equations did not factor there). It now starts the residual from the automatic start,
  and `message` says so. A `start=` that ends with the residual at its floor also gets a walk
  from the automatic start, and the better fit is kept: with one record per animal the corner
  `var(residual) -> 0` is a local optimum of the boundary, and in simulation (300 animals) a
  start at the floor certified it 91 units of -2logL above the interior optimum. Gate in
  `tests/testthat/test-partida-piso.R`.
* `model_mt()`: a scalar component held at the boundary (a variance at its floor in the
  certificate) has no standard error, its row and column of `vcov` are 0, and the standard
  errors of the other components come from the information conditional on it, as in
  `model()`; `message` names it. Gate in `tests/testthat/test-ep-condicional-mt.R`: with
  `pe(id, traits =)` at the boundary, the standard errors of the other components equal those
  of the model without it at the same point.
* The walk from the reduced model, when a covariance group ends near a correlation of +-1,
  looked for the groups near +-1 in a list that lost every group but the last: a new pair of
  components with no group in common emptied the list (`x[-integer(0)]` is empty in R). With
  two disjoint groups near +-1 only one was seen. Gate in `tests/testthat/test-lista-correcoes.R`.

* `qc_phenotypes(fence_by =)`: Tukey fences within each class of a column (an age class), so
  that a trait that grows with age is not cleaned by one fence over all ages. A class with
  fewer than 8 observed records keeps the fences of all records; the fences of each class are
  in `fences_by`.
* The warning about covariate values more than 10 robust SDs from the median also says to
  check them for a typing error, which is corrected in the data, and not only to treat them as
  a missing-value code.
* `model()`: an NA in the level column of a term without relationship (`random(pen)`, a
  fixed class) also says the term and how many rows hold an NA, as blank text already did.
* `print()` of `qc_phenotypes()` shows the summary: the result has the class `br_qc`, and
  the cleaned records stay in `$data`. Before, typing the result printed every record before
  the summary, into the console and any session log.
* `print()` of `model_mt()` gives the Newton decrement next to relDelta, as `model()` does,
  and says "record(s) with at least one trait" (it counted rows with any trait as
  "complete"). The coefficient of a covariate prints as `cov(x)` and not as `x=x`; `fit$b`
  keeps its names.
* `h2_curve(bases =)`: the base at each point of the gradient of a scalar term declared with
  `base =` (a `random()` on an indicator of a class of age, for example). The value of a data
  column at a point is not known from the fit; before, the variance of such a term entered at
  every point with base 1, without notice, and now that case warns.
* The check of suspect values before a fit no longer flags legitimate
  values of a concentrated variable. With more than half of the values at the median the
  robust scale was 0, and its reserve marked as "more than 10 robust SDs from the median" the
  real age of a weight at 60 days when almost every animal was weighed at 60, and the 1 of a
  0/1 trait of low incidence. The scale now comes from the values away from the median, and a
  variable of two values is an indicator, which the rule does not read.
* `gibbs()`: when every chain stops early, the message says that the results come from the
  longest chain and that it stopped too. Before, it named only the other chains, and the one
  used read as complete.
* `summary()` of a `gibbs()` fit gives the posterior of the ratio of each variance over the
  phenotypic variance of `h2()`, not only the ratio of the genetic term (the maternal and
  common-environment ratios had to be computed by hand from the samples). The location note
  of `mcmc_diagnostics()` comes out one time for the terms that share its limit, and says that
  more chains or a longer run lower the limit.
* `mcmc_diagnostics()` of a probit `gibbs()` chain: the location note says that a level with
  every record in one category has an improper posterior only for the fixed effects under a
  flat prior. For a random term it says that the prior of the term keeps such a level finite
  and that it mixes slowly.
* `model()`: the boundary notes of `message` name the components at the zero boundary and
  the groups at the singularity boundary, and say whether a group has a correlation at +-1 or
  a combination of its effects with no variance and no extreme correlation (in a group of 3
  or more). The note no longer says "the boundary is real", which the fit does not check.
* `model()`: `var(residual)` has a floor, 1e-10 of the variance of the records, the floor of
  the variances of the covariance groups, and the same rules: a residual that goes to zero
  is held there, the others are optimized conditional on it, it has no standard error, and
  `message` names it. Before, the residual walked in log toward minus infinity: with two
  equal records per animal the fit stopped at a residual of 9e-14 with a Newton decrement
  of 1e11, and the message did not say that the residual was at the boundary.
* `model()`: a scalar component resting at the zero boundary (the fit already said "fixed
  there while the others were optimized conditional on that") has no standard error, and the
  standard errors of the others come from the information conditional on it,
  `2 [AI_ll^-1]` on the free block. Its row and column of `vcov` are 0, so the delta method
  of `h2()` and `se_function()` reads it as a constant. The SINGULAR verdict reads the same
  free block. Before, all the standard errors came from the whole matrix.
* `model()` and `model_mt()`: near a correlation of +1 or -1 in a covariance
  group, the -2logL comes from a transformed system. With the group covariance C = F F', the
  effects u = (F x I) u* have Var(u*) = I x K, and the equations of u* have no C^-1. Before,
  log|G| + log|C| of the mixed model equations cancelled two large numbers, and the error
  grew as 1/(1 - |r|) and with the number of levels. With 1 - |r| = 1e-12 the error against
  the dense V form was 0.002 to 0.03 units in simulated data of 280 to 840 records, and the
  fitter could accept a step because the wrong value was lower. Now the two agree to 3e-8.
  The transformed system is used when the smallest eigenvalue of the correlation matrix of a
  group is below 1e-6; the score and the information matrix still come from the usual
  equations. Gate in `tests/testthat/test-lista-correcoes.R`.
* `model()`: when a covariance group of two or more terms (`group =`) ends near a
  correlation of +1 or -1 (the smallest eigenvalue of its correlation matrix below 1e-4),
  certified or not, the fit also takes the route that the message gave before: it fits the
  same terms without the group, then walks the full model from that optimum with the
  covariance at 0. The better of the two walks is returned, by the rules of the second walk
  of `model_mt()`, and `message` says which one. Near +-1 the score and the information
  matrix of the walk lose digits, and the walk creeps to the boundary or stops in another
  basin. In our own data set (direct-maternal, 3507 animals) the walk from the var(y) start
  stopped without a certificate at 300 and at 600 iterations, with r = -0.999998, 0.19
  -2logL units above the point of the reduced model; now the fit certifies that point. In
  simulated direct-maternal data with no genetic variance (12 data sets of 300 records), 11
  walks ended near +-1: one certified a point 0.0017 lower, and the other ten kept the first
  walk at about twice the iterations. A group of one term (`rn()`) has no reduced model and
  stays as before. Gate in `tests/testthat/test-lista-correcoes.R`.
* `model_mt()`: at the limit of re-entries into the cone, the message no
  longer advises the univariate fits with the covariances at zero when the fit already
  started from them (the default start). It says that, and advises a start with that
  covariance nearer to zero, or a likelihood-ratio test against the fit that holds it at
  zero.
* `model()`: the Newton decrement of the convergence certificate inverts the information
  matrix by its eigenvalues, each floored at 1e-9 of the largest. On a flat ridge the matrix
  is singular to the rounding, and the plain inverse gave 1e-6 at one point and 1.7e7 at a
  point with the same -2logL to 1e-8. A large gradient along a direction of no curvature
  still refuses the certificate.
* `model()`, `model_mt()` and the other Gaussian fitters stop when the trait
  has the same value in all records: there is no variance to partition. Before, the walk ran
  to `maxiter` with the residual going to 1e-17.
* The EM warm-up of `model()` puts the rank of X into the update of `var(residual)`:
  `(e'e + s2e (p + q - tr)) / n`, the EM of REML with the fixed effects as missing data under
  a flat prior. Before, the p was missing, and the fixed point of the EM was the residual of
  ML, not the zero of the REML score. With 60 fixed levels for 240 records, 4000 EM rounds
  went to a residual of 0.660 against 0.936 of REML; now they go to the REML optimum.
* A fit that stops at `maxiter` with relDelta 0 says that the walk froze, so more
  iterations give the same point, and no longer advises to raise `maxiter` (`model()` and
  `model_mt()`). No gate reaches such a stop yet.
* The help of `maxiter` says when a fit at the ceiling has `converged = TRUE`: the Newton
  decrement at the last point certifies a flat direction of the likelihood, and `message`
  says that the components were still moving. Before, the help said that a fit at the
  ceiling always had `converged = FALSE`.
* The formula reader takes an explicit intercept, `y ~ 1 + cg + animal(id)`, and a unary
  plus, `y ~ +cg + animal(id)`, as the same formula without them. Before, the first stopped
  with "did not understand the term: 1" and the second with "subscript out of bounds".
  `0 +` and `- 1` stop with the reason: the fitters always fit an intercept.
* `legendre()` takes any order. Before it stopped above order 6.
* A term declared twice (two `rn(id, ...)` on the same column) says how to name the second
  one, with `nome =`.
* A term without relationship whose column holds "0" among animal ids of the pedigree
  (`pe(dam)` without `maternal(dam)`) now warns: the term reads "0" as one level, so all the
  records with an unknown dam share one effect. A column of pens or lots has no warning.
* An inadmissible start of `model_mt()` names the group whose covariance matrix is not
  positive-definite, and the fit returns `neg2logl = NA` instead of 0.
* `sparse_chol()`, `selected_inverse()` and `sparse_solve()` stop on triplets of both
  triangles of the matrix. Before, each triplet off the diagonal counted for both of its
  positions, so those entries counted twice, and the result came from another matrix.
* The help of `weights =` in `model()` says that the -2logL leaves out `-sum(log(w))`.
* `model_mt()`: the message about components resting at a covariance
  boundary now names them (up to 8, and how many more). Before it gave only their number.
* `profile_theta()` gives `eval_internal()` only the arguments that it takes: a
  `verbose = FALSE` for `model()` made every evaluation fail, and every point came back
  with the penalty 1e12 and `converged = TRUE`. A grid point outside the admissible space
  (a variance of 0, a group outside the cone) is now `NA` with `converged = FALSE`, and a
  Nelder-Mead stop at `maxit` is restarted once from its point before it counts as not
  converged.
* `model_mt()` keeps `vcov`, `2 AI^-1` with the names of `theta`, so `se_function()`
  gives the delta-method standard error of r_g and h2. Before, `se_function()` refused with
  a reason that did not hold. A component with no standard error has its row and column
  `NA`. Gate in `tests/testthat/test-mt-traits.R`: the square root of the diagonal equals
  `se`.
* The SINGULAR verdict gives the advice of the components that it names: dams and
  granddams with records for `maternal()`. Before, a maternal model got the advice about
  pens.
* `gibbs(prior = "jeffreys")`: a chain whose scalar variance had fallen toward zero could
  stop the whole fit with "the matrix is not positive-definite", every chain lost
  (`validation/gibbs_scalar_prior.R`: `animal(id)` on 5 levels with 0.05 planted, 2 chains
  of 50000). The scale of the inverse Wishart draw, `u'K^-1 u`, had reached the subnormal
  range. That chain now stops with "the scale matrix of group 'animal' is not
  positive-definite inside the chain (its variance drawn near zero)", the verdict says
  `"chain stopped"`, and the chains that ran to the end are pooled. Gate:
  `tests/testthat/test-gibbs-padrao.R`, three of four chains started at a variance of 1e-307
  on a variance whose true value is zero (at least one stopped in each of 60 seeds, the
  fourth ran to the end in 60 of 60).
* A chain stopped by a failed Cholesky of the mixed-model equations said "C is not
  positive-definite inside the chain", which reads as the covariance matrix of a group. It
  now says "the coefficient matrix of the mixed-model equations is not positive-definite
  inside the chain (a variance drawn near zero)".
* The extreme-category refusal of `gibbs(family = "probit")` explains the case of the
  normal prior of the fixed effects: there the posterior of the level is proper but set by
  that prior alone, and the level is still refused.

* "no row enters the analysis" no longer covers different causes with one sentence.
  `model_mt()` said "no estimable fixed effect for any trait" when no row entered. In the
  joint mode, when no row has both traits observed, the message says so ("the joint mode
  takes only the rows with both traits ('q' and 'b') observed, and none of the 120 row(s)
  has both"); 12aca18 said "the SECOND trait of the joint analysis must be binary (it has
  0 value(s))". The joint mode's own specification errors still come first.
* A missing dam or sire (NA, blank text) in a relationship column is still refused, and
  the message now says that an unknown parent is written 0 and keeps the record in the
  analysis, instead of only "drop those rows or fill them".
* A genotyped id written differently from its pedigree line (a trailing carriage return,
  say) stops `h_inverse()` and every single-step fitter with both spellings. 12aca18 said
  "genotyped animal 'G1_001' is not in the pedigree", printing the id without the
  character that made it differ.
* Genotypes against the pedigree and against themselves (survey of data-structure cases,
  group E2, at commit f2c0d69): a repeated genotype id, rownames that contradict the ids,
  an animal imputed mostly by the marker means (a failed sample or a low-density panel,
  which shrinks its genomic relationships and makes it look average and more accurate than
  it is), parent-offspring pairs that Mendel forbids and `read_plink()` ids changed by
  `read.table()` all ran without a word. They are now refused or handled as described
  under New and Breaking changes; the gates are in `tests/testthat/test-estrutura-e2.R`
  (27 tests, every one fails on the build before this change).
* A bare matrix given to `genotypes =` stops with "genotypes must be a list with 'ids' and
  'm'" (it was R's "$ operator is invalid for atomic vectors").
* Covariates enter the fixed design standardized (mean and SD of the used rows) in every
  fitter, and the rank is decided on that scale. `b`, `se_b` and the `gibbs()` posterior mean and SD of `b` come back in the unit of the data, and
  -2logL keeps its value. Before:
  - age + 1e9 left theta, EBV and b all NA;
  - age + 3e8 gave b = 5.90;
  - age + 1e11 and age x 1e-12 went to `dropped_x`;
  - `gibbs()` with age + 1.7e9 gave b = 0.

  Dates and timestamps are now fitted at full precision. For the nested covariate the
  warning says that beyond about 1e6 SDs the fit may not converge, and beyond about 1e7
  SDs the equations may fail to factor (reported as an inadmissible or singular starting
  theta, or as C not positive-definite in `gibbs()`), then gives the advice: add the
  class as a fixed effect, or subtract a reference value. A nested covariate whose class
  is a fixed effect aliased with another one is also scaled only, with no warning: `herd
  + cg + cov(x, nested = "cg")` with cg nested in herd, where the rank test drops more
  than one level of cg that has used records. A very large offset can still make such a
  fit fail as at 12aca18; fit cg before herd, or leave herd out.
* The help of `gibbs()` said that a design that does not identify a component stops the
  chain. The structure checks stop it only in the structural cases they list. The help
  now says so. Its `prior` entry no longer says that the spurious mass at zero is
  negligible for a well identified component: that holds for a scalar variance, not for a
  covariance group.
* The messages of these checks quote the labels they list, so a blank or a leading space
  shows ("'2', ' 2'"). A list cut at its first items ends in "(and N more)".
* `?predict.breeding_fit_thr`: the joint fit returns the probability of the upper binary
  category, `fit$binary_levels[2]` (the second level of a factor), not "the larger binary
  value".
* Extractors and separability (survey of data-structure cases, group E6). Reviewed three
  times, with the gates in `tests/testthat/test-estrutura-e6.R` (56 tests) and the
  reproductions in `validation/structure_cases.R`.
* `rg()` documents its reading of the fit state in Details (the text had landed inside the
  `t2` argument).
* The test suite uses `expect_warning(x <- ...)` instead of the value of
  `expect_warning()`, which is the condition under testthat 3.2.1 (DESCRIPTION
  accepts >= 3.1.7), and a gate refuses that idiom.
* Blank text and the text "NaN" in a level column were levels of their own, in silence.
  - Blank text is what `read.csv()` and `data.table::fread()` give for an empty cell of a
    text column; spaces, tabs and line breaks alone count as blank. "NaN" is what
    `factor()` makes of a numeric `NaN`.
  - What it did at f2c0d69: with "" in two rows of `random(lote)`, -2logL was
    54.430335165401253, bit for bit the value of a level named "BRANCO" in those rows. A
    blank `cg` gave the column "cg=" (in the reference class). A blank
    id of `animal()` left the fit like an id outside the pedigree (n_used 46 of 48).
  - Every fitter now stops with the term, the column, the number of such rows and the
    first of them, on every row.
  - The texts "NA", "Inf", "-Inf", "nan" and "0" remain level names, since they may be
    real ones.
  - In the maternal-grandsire column of `sire(mgs =)`, blank text is the unknown
    grandsire, like "0". A single space " " used to stop with "maternal grandsire ' ' is
    not in the level set". The text "NaN" there is refused.
* A covariate with a missing value left the model without a word.
  - What it did at f2c0d69: `cov(x)` with one `NA` in 48 rows became NaN and left through
    the rank check (`dropped_x` "x=x"). The fit had the -2logL of the model without
    `cov(x)`, 50.430239486919447, with `n_used` 48. A nested covariate lost only that
    row's class ("x=c1"). `random(g, base = "phi1")` with one `NA` gave -2logL 0. `Inf` did
    the same as `NA`, and `model_mt()` and `gibbs()` dropped the covariate too.
* `pegs()` made a blank environment one more trait (`n_records` 2 39 39). An `NA`
  environment dropped its records uncounted (`dropped` 0), and a blank or `NA` id was
  counted as a record without genotypes. The id and the environment now follow the
  level-column rule.
* `predict()` in the joint mode without a class factor stopped with "no column '' in
  newdata". It now uses the intercept.
* The rule in R compares bytes, as the engine does. It used `trimws()`, which stopped on
  a label that `read.csv()` read from a latin1 file without `fileEncoding`: "input string
  k is invalid UTF-8", with no term, column or row.
* `qc_phenotypes(classes =)` returns `missing_labels`: per class column, the number of rows
  with a missing label that every fitter refuses. Those rows no longer count as a class in
  `small_classes`.
* Weights: a text or factor column is refused instead of being read as codes. `NA`, `Inf`
  or a weight <= 0 is reported with the column, the count and the first row.
* `?model` states the rule once ("only the observation may be missing"), with the
  conventions of the field:
  - ASReml stops on a missing explanatory value unless `!MVINCLUDE` or `!MVREMOVE` is
    declared; asreml-R defaults to `na.method(x = "fail")`.
  - A BLUPF90 class code 0 becomes one level "0" here.
  - `missing_code`, like `OPTION missing` in BLUPF90, marks the observation only: a
    covariate holding -999 is fitted as the number -999.

  Gates: `tests/testthat/test-rotulo-ausente.R`. Every case has an arm that runs the same
  data with a real label in place of the missing one and checks the number the missing one
  gave at f2c0d69. One battery of values runs through both implementations of the rule,
  the engine's and R's.
* `model_mt()` no longer certify a covariance boundary that the step merely stopped at.
  Their step freezes a boundary block whose curvature has collapsed, whatever the sign of
  the score, and the certificate left that block out of the Newton decrement. On real
  data two bivariate fits came out `converged = TRUE` 46.8 and 16.3 -2logL units above
  the interior optimum. A certificate with a block at a boundary now has to pass four
  checks on the likelihood itself:
  - the Karush-Kuhn-Tucker condition: points inside the cone along the negative
    eigenvectors of the block's gradient matrix, with the step length of the quadratic
    model of the average information, and as a fallback each Cholesky diagonal raised to
    0.1%, 1% and 10% of its scale;
  - the boundary face: the coordinates of the block that move along it, by central
    differences of the -2logL (step 1e-2 in log-Cholesky), in a joint Newton step with the
    free components;
  - a second walk from inside the cone, within the same `maxiter`, for every fit that ends
    with a block at a boundary: each such block made diagonal, its variances at least 10%
    of the `var(y)` start;
  - after a certificate passes, a polish at the boundary: the face step again, and for a
    trait whose variance sits at the floor the same point with its covariances at zero,
    accepted while each lowers the -2logL by more than 1e-5 and the new point still
    certifies.

  A point that lowers the -2logL by more than 1e-4 sends the walk there, or refuses the
  certificate when the walk cannot move (at `maxiter`, after 10 re-entries into the cone
  or 20 moves along the face). The polish never turns a certified fit into an uncertified
  one and never raises the -2logL.
  - Simulated starts at r = +-(1 - 7.5e-7), which the previous build certified 5.2 to 5.5
    units above the optimum in `model_mt()` and 15.0 above in now reach it;
    with `maxiter = 20` the same start ends `converged = FALSE`. Of 48 starts with one
    litter variance at 1e-12 (seeds 1 to 12), the previous build certified 17 points above
    the optimum, the KKT check alone 6 (0.0056 to 4.73 units above), and with the second
    walk none.
  - These are searches, not a proof. A better point that none of the checks reaches still
    passes.
  - `model()` keeps its check by the sign of the score: its step never freezes a component
    whose score points inward.
* `model_mt()` walk on each trait divided by the power of two just above its standard
  deviation, and evaluate the result in the original units.
  - The walk no longer reads the unit of a trait. The relative floor of the step, the
    curvature-collapse test, `relDelta` and the slacks of the EM rescue all used to compare
    numbers in different units.
  - On simulated data the previous build certified a point 4.70 -2logL units above the
    optimum, with theta 19% off, for y1 times 1000. For y1 times 1e-4 it stopped without a
    certificate 133.7 units above.
  - Now the estimates agree up to the scale factors (1e-5) and the -2logL up to the
    Jacobian `2 (n - p) log c` (1e-9); at 24 000 records, 3e-8 and 1.8e-9 for y1 times
    0.01. A power-of-two factor gives the same fit bit for bit.
  - Components resting at a boundary are the exception. They sit on a flat face and can
    move along it when the data round differently: 0.4% and 0.6% with ten components at a
    covariance boundary, the -2logL equal within 2.5e-5 after the Jacobian.
  - On real data, a boundary estimate had moved 0.66% and the -2logL 0.011 under a
    rescaling of one trait. That case was not rerun with this version.
* The SINGULAR verdict and the standard errors no longer depend on the unit of `y` in
  `model_mt()`.
  - A component with no information at all is still flagged. A correlation form of the
    matrix, tried first, normalized it away (diagonal 3.1e-31, standard error 9.6e15, no
    warning).
  - `model()` is unchanged here, since all its components share the unit of `y^2`.
* `model()`: the variance floor of the step is now 1e-10 of the variance of `y`, fixed for
  the whole walk.
  - It was 1e-10 of `max(1, s2e)`. In the user's unit that held the floor at 1e-10
    whenever the residual variance was below 1.
  - With `y` times 1e-5 a permanent-environment variance stuck on that floor, and the fit
    was declared converged with the additive variance 14% low. With 1e-6 it stopped 1.685
    -2logL units above the optimum.
  - Now `y` times 1e-6, 1e-8 and 1e6 give the same estimates to 1e-11.
  - A variance at the zero boundary is now reported at 1e-10 var(y) instead of 1e-10.
* The REML -2logL of `model()` and `model_mt()` is formed from the residuals, `e'R^-1 e +
  u'G^-1 u`. It used to be `y'R^-1 y - b'r`, which subtracts two numbers of order n (mean
  y)^2 and lost digits when a constant was added to `y`.
  - All of them now agree to 3e-14, and to 6e-15 at 21 000 and 24 000 records
    (`validation/unit_shift_start_scale.R`).
  - On real data the collapse had been off by 9e-6.
* `maxiter = 0` evaluates the start and returns `converged = FALSE` in `model()`,
  `model_mt()`: an evaluation is not a convergence. In `model()` this is a change: from a
  start at the optimum, with `n_em = 0`, it used to return `converged = TRUE`, with the
  message "converged by the Newton decrement while the components were still moving
  (relDelta 0)". `model_mt()` returned `TRUE` as well.
  - `model()` says when the evaluated point was moved by `n_em` EM rounds.
  - A start beyond the floors of the step, such as a correlation within about 1e-6 of
    +-1, is raised to them before the walk in `model_mt()`. There the
    sparse -2logL was 0.0028 to 51 units below the dense V form, and the fit stopped at its
    first iteration.
* `model_mt()` chooses the reference level of a class among the levels with records of
  every trait.
  - A drop for one trait (`term=level|trait`) is therefore a level with no record of that
    trait.
  - `dropped_x_reason`, also returned by gives the reason of each drop: "no
    record", "aliased", "no record of y2" or "aliased within y2".
  - The previous build named the empty level `CG=x`, without the trait, and dropped a level
    that had records as `CG=y|trait`.
  - Rank, likelihood and breeding values do not change.
* `geweke_z()` divided the difference of the window means by their iid variances.
  - On the four chains of a converged `gibbs()` of 860 animals it gave z = -4.12, 4.93,
    2.56 and -6.01 for `var(animal)`.
  - It now uses the spectral density at zero (Geyer's initial monotone sequence),
    estimated in the final window only, and combines the chains. On the same run the
    combined z is -0.09.
  - `gibbs()$geweke` looked at the first chain only; it is now the combined z.
  - The `|z| < 4` check in test-gibbs.R had hidden this.
* With several chains, `gibbs()$ess` was the sum of the per-chain ESS, which cannot see
  chains that disagree. It is now the bulk-ESS of Vehtari et al. (2021).
* `ess()` is now Geyer's initial monotone sequence on FFT autocovariances, as
  `posterior::ess_basic(split = FALSE)`, and takes a matrix of chains. The previous
  estimator (positive pairs of `acf()`, no monotone step) agreed in the median. On single
  chains of 2000 draws, 10% differed by more than 5% to 8%, and the worst by 24% to 39%
  (200 chains per autocorrelation).
* `rhat()` folded `|x - median|` after splitting the chains, with the median of the
  halves; posterior folds first. With an odd number of draws they differed in 18% of the
  runs of 4 chains of 21 draws, by up to 2.4%.
* A chain of `gibbs(chains =)` that stopped early (an inadmissible theta or a failed
  Cholesky inside the chain) kept fewer draws. `cbind()` recycled it silently in `rhat`,
  and the pooled locations averaged it in. Now:
  - only the chains that ran to the end are pooled and diagnosed (the longest alone if
    none did), and the residuals and fitted values of the fit come from the same chains;
  - the verdict says `chain stopped`;
  - `chain_n`, `chains_used` and `chain_messages` say what each chain did.
* `gibbs()` with several chains: the pooled `ebv_sd` and `b_sd` are now the sd of all draws
  stacked (denominator N - 1). They were the root of the mean within-chain variance plus
  the population variance of the chain means, larger by about (M - 1) / (2 M n) for M
  chains of n draws (0.07% for 3 chains of 500). `ebv` and `b` do not change.
* `gibbs()` with one chain returned no `rhat`; it is now the split R-hat of the chain's two
  halves.
* `gibbs()` accepted `thin = 0`, which crashed R, and non-integer `n_iter`, `burnin` and
  `thin`; these are refused now.
* `validation/gibbs_harville.R` takes the Monte Carlo error of the means from
  `mcse_mean`; the z against the quadrature is unchanged (0.15 and 0.08).
* DESCRIPTION declared `testthat (>= 3.0.0)`, but the tests use `expect_no_warning()`
  (testthat 3.1.5) and `local_mocked_bindings()` (3.1.7); the bound is now 3.1.7.
  `posterior` and `coda` join Suggests.
* `estimate_gamma(method = "pseudo_em")` failed with "not positive-definite" when a
  metafounder had no genotyped descendant. As documented, its row and column are now zero,
  and the rest of Gamma equals the estimate without that metafounder (gate: 1e-10, for a
  metafounder that is only declared and for one whose descendants are all ungenotyped).
* `snp_effects()` with a marker that has every genotype missing stopped with "missing
  value where TRUE/FALSE needed". The marker is now left out, as it is in the fit.
* A genotype list whose `ids` and rows of `m` differ in number is now refused with the
  same message in every function that takes genotypes ("genotypes$ids has N identifier(s)
  for M row(s) of genotypes$m"). Before:
  - `apy_core_select()` chose a core from the mismatched list.
  - `h_inverse()` and `model()` already refused it, with a different message.
* A covariance group of terms without a relationship matrix, `random(a, group = "g") +
  random(b, group = "g")` (or `pe()`), paired its levels by position. Each term took the
  levels of its own column in order of appearance, so level i of `a` covaried with level
  i of `b`, which is another label, and the fit depended on the row order. On 300 records
  with 20 + 20 levels and fixed components, -2logL was 340.573 in the original order and
  340.078 with the rows permuted, against 338.860 from dense mixed-model equations paired
  by name. With 8 + 16 levels the Kronecker block had the size of the first term (the MME
  route gave 340.868, the V form 328.760), and the breeding values of the second term
  came out under the first term's labels. The terms of such a group now index the union
  of the level names of their columns, sorted (integer labels in numeric order, then the
  rest), and pair by name. A level present in one column only is still an effect of the
  other term, with no record there. Results no longer depend on the row order in
  `model()`, `model_mt()`, `gibbs()` and `snp_blup()` (a row-permutation gate for each).
  `model()` also equals the dense equations paired by name to 1e-10, and a bivariate
  `model_mt()` with zero between-trait covariances equals the sum of the two univariate
  fits to 1e-10. In `snp_blup()`, permuting the rows moved the group's breeding values by
  up to 1.28 (sd 0.31). Refit every fit with a multi-term iid group made before this
  version: `accuracy()` refuses such an old fit only when its blocks carry different
  level names (8 + 16); one whose two columns had the same number of levels names both
  blocks with the first term's levels and passes without an error, still on the old
  pairing.
* A group of iid terms whose columns share no level with records in both (sire and dam
  ids that never repeat between the sexes) has a covariance that does not enter the
  likelihood. REML used to run its 300 iterations and end with `converged = FALSE` and
  NaN standard errors for every component. `model()`, `model_mt()` and `gibbs()` now stop
  before fitting and name the group, the two terms and their columns. With given
  components (`start =` with `maxiter = 0` and `n_em = 0`, or `theta_fixed =` in
  `gibbs()`) the group is accepted.
* `accuracy()` divides each level by its own prior variance, matched to the PEV by level
  name. These groups need no pedigree, and `solutions()` gives their `acc` column without
  one (in a multi-trait fit, with `trait =`).
* `accuracy()` no longer matches the pedigree by position. The pedigree was rebuilt and
  read row by row, so the same pedigree with its rows in another order divided animals by
  the F of other animals (up to 0.065 on 210 simulated animals). The single-step prior
  `h_prior` comes named by genotyped animal from `model()`, `model_mt()`, `h_inverse()`
  and is placed by name. The pedigree must still hold exactly the animals of the fit.
* `eval_internal()` and `eval_internal_mt()` report a misnamed column as "no column(s) in
  the data: x", as the fitters do, instead of "undefined columns selected".
* Numeric ids are written the same way on both sides. The engine labels a numeric data
  column by the full integer ("100000"); the R side wrote the numeric ids of the
  pedigree, the genotypes and the keys with `as.character()`, which gives "1e+05" for the
  double 100000, and the records of those animals left the fit, visible only as a smaller
  `n_used`. On 150 animals with ids 1e6 and 2e6 among 100002 to 100149, two records each,
  `n_used` was 296 of 300 and -2logL 434.12, against 300 and 441.69 with the same ids as
  text. Every numeric id now goes through the engine's own formatter, so double, integer
  and text ids, and mixed types across the pedigree columns, give the same `n_used`,
  -2logL, components, breeding values and accuracy (gated to 1e-10), in `model()`, the
  single step, `snp_blup()` (double ids in `subjects` and integer ids in `changes` used
  to stop it). Fits made before with such ids named their levels "1e+05": refit before
  calling `predict()` on them.
* A non-integer numeric id was labelled with 6 significant digits: 123456.7 became
  "123457", the label of the integer 123457, and in the pedigree the animal 1.1234567
  and the sire 1.1234568 became one animal. Labels are now the shortest writing with 15
  to 17 significant digits that reads back to the same double, so distinct numbers never
  share a label.
* `qc_genotypes()` accept a raw genotype matrix (1 byte per genotype, 5 for missing), as
  `read_blupf90_snp()` and `read_plink(storage = "raw")` return it, with the same
  filtering, counts and values as a double or integer matrix; a raw matrix used to be
  refused with the number of entries of the whole matrix reported as out-of-code
  genotypes. `qc_genotypes()` returns the filtered matrix in the storage it came in.
* `apy_core_select(method = "lanczos")`: the earlier estimate sat 0.3 to 0.5% above the
  exact count because of the Gauss quadrature at the threshold, not because of the
  probes, and its `count_se` was the spread of the per-probe counts at a fixed
  threshold, which jumps by a whole Ritz weight: at 10 000 x 10 000 it was 69 against an
  actual spread of 12.5 over ten seeds, so "within one standard error" tested little.
  The averaged rule and the new standard error fix both (paired quadrature error at 98%
  from +21.4 to +2.3 counts, standard error 12.2 against a spread of 11.9). The
  "conservative side" claim is withdrawn.
* `h2_curve()` summed the components of the other groups raw. With the permanent
  environment also a random regression on the same basis, the curve of Example
  10.2 of Mrode and Pocrnic was off by a ratio of 0.63 to 1.04 along the
  lactation, so its shape was wrong. Every random regression is now evaluated at
  the point of the gradient.
* `model_mt()` measured the rank of X before dropping the rows whose levels have no line
  in the pedigree, as `model()` does not: a contemporary group made only of such animals
  reached the equations as an empty column, and the fit said "did NOT converge" or
  "SINGULAR" about the wrong thing.
* On a flat ridge with the Newton decrement already under tolerance, `model()`,
  `model_mt()` returned `converged = FALSE` at `maxiter`. The verdict at the end of the
  budget now follows the certificate, and the message says the likelihood is flat along
  some direction. In the two mirrors the "stopped at maxiter" diagnostic was also hidden
  by any informational note.
* The rank of X is measured PER TRAIT in a multi-trait fit, so a design where a
  fixed level only occurs for one trait fits instead of dying. The unit that
  drops is the pair (column, trait), reported in `dropped_x` as `CG=5|y1`; after
  the empty pairs go, what remains for a trait is checked for collinearity too.
  Gated three ways: the sparse and the dense V routes agree to 7.5e-12 on the
  design that could not be fitted before, the components land where the separate
  single-trait fits land, and a design with nothing to drop keeps its layout. The
  multi-trait evaluation also says whether it failed on the covariance or on the
  design, which used to share one message.
* `summary()` existed for one of the seven fit classes; on the others it fell
  through to `summary.default`, which returned a table with the shape of a result
  and no meaning. `summary()` and `print()` also divided by different
  denominators on the same fit; both now come from one table.
* `gibbs()` accepted metafounders together with genotypes, which every other
  fitter refuses because the genomic side does not know Gamma; now it refuses too.
* The C++ engine printed some messages in Portuguese (the starting theta, the EM
  rescue, the factor pattern); they are in English now.

## Version

* The version carries a development suffix. `v0.4.0` is a tag, and what comes
  after it is not `0.4.0`.

# uaitrem 0.4.0

## New

* `start = ` in `model_mt()`, which its own stopping message had been recommending with
  no way to act on it.
* Metafounders take a full Gamma matrix and not only its diagonal, so the
  ancestral relationship BETWEEN two base populations is a parameter rather than
  an assumption. A singular Gamma goes through the pseudo-inverse instead of
  being refused: `gamma = 0` is the unknown-parent-group limit.
* Every iterative fitter prints where the estimates are while it runs, grouped
  by term, so a covariance matrix reads as a matrix.

## Fixed

* The convergence certificate in `model_mt()` was computed in theta while the step had
  already moved to log-Cholesky coordinates, with an active set of its own. A clamped
  direction was frozen in the step and charged in full in the certificate: a fit sat 980
  iterations with the decrement stuck at 3.25e+03 at a point a coordinate descent
  improved by 0.200 units.
* `accuracy()` divided by the prior the animal actually has, in a two-term
  covariance group.
* The mirrors measured the rank of X over the whole table instead of over the
  rows that entered the equations.
* A marker argument with a misspelled name passed silently. `animal(grupo = )`
  now says the argument does not exist, suggests `group`, and lists what the
  marker takes.
* Metafounders combined with genotypes are refused rather than silently
  correcting the base twice.
* `br_version()` reported 0.1.0 while the package was at 0.3.0: the string is a
  constant in `src/entrada.cpp` and nothing tied it to DESCRIPTION. It is now
  compared against `packageVersion()` by a gate, so the two cannot drift again.

## License

* GPL-3, replacing MIT, decided before the repository was ever made public and so
  with nothing already granted to anyone. MIT let anyone take the engine, extend
  it and ship the result closed; GPL-3 asks that a distributed derivative carry
  the same licence and ship its source. Research use and publication are
  unaffected, which is where citation comes from.
* `CONTRIBUTING.md` asks a contributor to grant the right to license their
  contribution under a licence the project may adopt later. Without it the
  copyright spreads with the first merged pull request and the licence can no
  longer be changed without finding every author. Nothing about it is
  retroactive: a released version keeps the licence it was released under.
* The MIT two-line stub in `LICENSE` is gone, because it exists only for template
  licences. `License: GPL-3` in DESCRIPTION is a standard R abbreviation and needs
  no file: `LICENSE.md` carries the full text for anyone reading the repository,
  and stays out of the tarball as R asks, since R already ships a copy.

## Gates

* `R CMD check` runs on every push through GitHub Actions, on macOS, on Windows
  and on three versions of R under Ubuntu.

## Docs

* Every help page has a `\value` section, and a reference to another function is
  a link rather than the literal text `[model()]`: roxygen markdown is on.
* `CONTRIBUTING.md`, plus `URL` and `BugReports` in DESCRIPTION, so the package
  says from the inside where it lives and where a problem should go.

# uaitrem 0.3.0

## License

* MIT, consistently declared in DESCRIPTION, LICENSE (the R two-line
  convention), LICENSE.md, CITATION.cff and .zenodo.json (was other-closed).
  FAPESP-funded work ships open source.

## Docs

* Counts counted: 54 exported functions, 48 exercised by the hands-on vignette,
  every number re-derived by execution. The reference list cites like a
  reference list again - the when-this-would-matter notes moved out of it.
  Suite: 860 expectations, 0 failures; R CMD check on the
  tarball: Status OK.

# uaitrem 0.2.0

The package was reviewed chapter by chapter against Mrode & Pocrnic (2023), every
published example in scope became a test, and what the book covers that the package
lacked was built. The suite stands at 828 expectations, 0 failures; `R CMD check` on
the tarball reports Status OK with 0 errors, 0 warnings, 0 notes.

## New

* Fixed-effect solutions. Every fitter returns `fit$b`, a named vector
  (`term=level`, plus `|trait` in the multivariate fitters), sliced from the
  solution the solver already had. `coef(fit, effects = "fixed")` reads it;
  `print()` and `summary()` show the block. The parametrization is
  implicit-intercept (a dropped column is a zero solution), so comparisons with
  textbook solutions are by contrast; this is documented where `b` is.

## Fixed

* `converged` is now certified. Besides the relative step, the Newton decrement of the
  free components (`g' AI^-1 g`, boundary components excluded) must pass tolerance.
  Previously a fit could report `converged = TRUE` with the score far from zero when the
  AI and EM steps stalled together on the boundary of a 2x2 covariance group: measured in
  the field at 6.9 units of -2logL short of the optimum.
* `profile_theta()` profiles. At each grid value the other components are
  re-optimized; it used to evaluate a slice through the estimate, which made
  intervals systematically too narrow.
* `accuracy()` runs on metafounder fits and uses the metafounder F; with
  gamma = 0.7 the difference against the unrelated-base F was 0.38 of accuracy.
* `maxiter` default raised to 300 in the iterative fitters, and the
  iteration-cap message now reports the step, the decrement, and how to ask for
  more.

## Gates

* New test files carry the published numbers of Mrode & Pocrnic (2023), chapters 3, 4, 5,
  8, 10, 15 and 17-19, page-cited constant by constant.

## Docs

* Every named method now carries author and year, in prose, roxygen and code comments
  alike, and the reference list holds 62 verified entries, identical in README.md and
  REFERENCES.md. Two real credit errors were fixed: the matrix-free A22^-1 identity is
  Masuda et al. (2017), not Vandenplas, and the third author of Takahashi et al. (1973)
  is Chin, not Chen.
* The full maternal model taught by the README and the skill is Eqn 8.1: ONE
  permanent environment, the dam's. The old recipe with a second `pe()` on the
  animal is non-identifiable on single records and is documented as such.

# uaitrem 0.1.0

Initial version: AI-REML with covariance groups, single-step genomics (H, APY,
ssSNPBLUP), reaction norms, direct-maternal effects, multi-trait, a Gibbs sampler,
and the supporting pedigree and quality-control tools.
