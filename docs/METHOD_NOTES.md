# Method notes

How food**net** turns mGrowthDB monocultures into produced and consumed arcs, the decisions behind each
step, and the questions still open. A change to any of these is a scientific decision: it gets an issue
where a human settles it before it is built (see [AGENTS.md](../AGENTS.md)).

## The specification (Karoline, 2026-10-04)

food**net** is the sister tool of grow**net**: the same structure (page and command line), the same
inputs, the export to graph formats, Cytoscape and R, the option to collect growth rates, and merging to
genus and across studies as advanced options. Its purpose differs: "it builds bipartite taxa-metabolite
networks that show which taxon produces and consumes which metabolites. It works with batch monocultures
with metabolite data."

## Decisions

1. **Growth phases, not a fixed window.** Karoline: "metabolite production and consumption during
   exponential phase can differ from stationary phase. So let's do this differently: estimate when the
   exponential phase ends and take that as the threshold", with a radio choice "Exponential phase" (the
   default), "Stationary phase", "Both", and "This will not treat diauxic shifts well, but we are also not
   able to identify them always clearly, so OK for now." A time window in the advanced settings overrides
   the phases, with a start and an end.
2. **Where exponential growth ends** (the agent's choice, to confirm): the first sampled time at which the
   culture reaches 90% of its maximal abundance on the linear scale, counted from its start. The share is
   an advanced setting. The linear scale was chosen on data: on R. intestinalis in study SMGDB00000007, a
   90% rule on log abundance ended the phase at 12 h while glucose was still taken up until 16 h, when the
   cells peaked. The growth curve is per strain when the replicate has one, else the culture's own, with
   cell counts before optical density. A replicate without a curve borrows the median boundary of its
   experiment, and says so.
3. **The value is the net change.** Karoline: "net change is OK": the concentration at the end of the phase
   minus at its start, interpolated linearly between metabolite samples (they are often sampled at other
   times than the cells), never extrapolated past the last sample (`window_beyond_data`). Averaged over
   replicates. The split into phases replaces the need for a two-phase rule: "it does not produce/consume in
   parallel but sequentially" (Karoline, on R. intestinalis).
4. **The detection limit is 0.2 mM**, an advanced setting, applied to the mean change (where the data of the
   hand-checked reference matrices break).
5. **Short records are used and flagged.** Karoline: "If metabolites are recorded for less than 24h, they
   are still being used but there should be a warning about it."
6. **The value medium.** Karoline: "It should be the medium that supplies data for most of the taxa,
   whichever medium that is. However, if the 2nd input field is filled, then it's the medium (or media) of
   that field." Other media give presence or absence only: "the majority medium was used for values and
   the other media were only used for presence/absence. This is a good default, but an advanced option can
   allow ignoring media differences and reporting values regardless. Another advanced option should allow
   reporting everything as booleans." A tie for most taxa goes to the medium with more replicates and is
   reported. Medium names are compared without case, punctuation or a parenthesized abbreviation, since
   studies spell Wilkins-Chalgren with and without "(WC)".
7. **Pooling across studies, with conflicts reported.** Karoline: "can be pooled with the same medium, but
   when there's a conflict, this needs to be reported." A conflict is experiments that disagree on the
   direction, or on whether the compound changed beyond the limit.
8. **Duplicate deposits are counted once** (the agent's choice). Found while building: study 2's RI_WC is
   study 7's ri2 rounded to two decimals, and study 7's bt2 is study 2's BT_WC. Two experiments of one taxon
   in one medium are one deposit when every replicate series of one matches a series of the other on every
   shared compound.
9. **Acid and base forms are one metabolite** (the agent's choice). mGrowthDB records acetic acid
   (CHEBI:15366) in some studies and acetate (CHEBI:30089) in others; without joining them, the same
   compound in two studies would be two nodes. The pairs are listed in `src/foodnet/compounds.py`.
10. **The matrices.** Karoline: "Instead of an adjacency matrix, there will be 2 matrix formats: one matrix
    with taxa as rows and metabolites as columns and another with 2 matrices: 1 for consumption and the
    other for production." NA is never zero. With Both, a column per metabolite and phase.
11. **CRM parameters.** Karoline: "raw plus helper is fine, and yes, send initial medium concentrations.
    Growth rates can come in that order: all replicates for which metabolite measurements are taken;
    another experiment in the same medium." The R package builds the efficiency matrix; scaling is the
    user's choice. With Both, the CRM uses the exponential phase.
12. **No metabolite is left out by default.** Karoline: "we'll handle this specific problem in mGrowthDB,
    we don't want to skip any metabolites by default" (on ethanol, whose HPLC channel in study
    SMGDB00000011 reads the ethanol used for cleaning). "Leave out these metabolites" stays as an advanced
    setting, empty by default.
13. **A filled second box is a limit.** Karoline, 2026-10-05, after "when I gave a list of studies, the
    results also included studies that were not in my list": "by default, when something is entered in the
    2nd field, only data matching what was entered are shown (I think this comes closer to what users want),
    but in advanced settings, we can switch on showing supporting evidence from other studies. By default,
    when the 2nd field is left empty, always all data are considered as discussed." The growth-rate fallback
    stays inside the limit too. Refined the same day, when four study ids pooled mMCB into the Wilkins-Chalgren
    values while reproducing the hand-checked reference: ids set the scope (which data), a medium name the value medium within
    it, and with ids alone the majority rule runs inside the scope, its other media giving presence.
14. **A second time window for named metabolites.** Karoline, 2026-10-05: "to include the different time
    window for trehalose ... I'd like to avoid a table where each metabolite gets 1 time window. perhaps 2
    time windows, with the option to add metabolites by name". The named metabolites take the second window
    instead of the phases or the main window; an empty end is each culture's last sample. Names match the
    metabolite or the name its study recorded it under. Columns, the image, the README and the report say
    which interval each value covers. With 0 to 48 h and trehalose to the end, every trehalose value of
    the hand-checked reference is reproduced.
15. **Values are used as mGrowthDB serves them, zeros included.** Karoline, 2026-10-05, on the isolated
    zeros in study SMGDB00000009 (B. fragilis formate at 4 h and 48 h): there is no reason to treat a zero
    as a missing value. Formerly open decision 3.
16. **Entries seen only in another medium are NA by default, TRUE or their value on request.** Karoline,
    2026-10-06: NA "is a cautious default. I'd like to have an option to set them to TRUE and another option
    to show the value but with a different background color in the image". The value is the mean over the
    other media's replicates; the image draws it outside the gray scale, since it is not comparable with the
    value medium's amounts. A CRM takes TRUE entries as missing.
17. **Studies are found through the species list**, not mGrowthDB's search, which matches per-strain
    measurements only and so misses a monoculture measured at the culture level (study SMGDB00000009).

18. **A medium is told apart by its alterations and atmosphere.** Karoline, 2026-10-06, on "Wilkins-Chalgren"
    also matching Wilkins-Chalgren with mucin: users "will not put by themselves '-mucin' unless they have extra
    knowledge ... we're basically hitting a limitation in mGrowthDB, which does not systematically report
    medium composition. To fix the issue ...: instead of just matching, check descriptions that suggest
    something altered the medium or atmosphere and treat it as another medium if you do. This default can be
    switched off in the advanced settings." `foodnet.media` reads what a description or name says was added or
    taken away, and the recorded gas composition (an unrecorded one does not split a medium). A medium name in
    the second box gives values from the medium it matches for the most taxa; the others it matches give
    presence only and are named. Experiments can also be excluded by id.
19. **Settled 2026-10-06 (Karoline), formerly open:** the end of exponential growth at 90% of the maximal
    abundance on the linear scale ("OK"); CRM parameters with Both use the exponential phase ("OK, but
    document": the help, the CRM README and the README say so); single-replicate values are shown with a
    caution ("OK"; since decision 20 a single replicate decides nothing by default, and is shown so only when
    judging by the mean alone); the growth-rate fallback takes another monoculture of the taxon in the same
    medium ("keep as is, but with the stringent medium matching"), which decision 18 now makes strict.

20. **The experiment is the unit, and a change must be pinned down.** Karoline, 2026-10-06, after a review
    of 0.1.0 showed a study with ten replicates outvoting one with two, and +0.84 +/- 1.04 mM on a 27 mM
    background drawn as production: "experiment level + inconclusive". Each experiment is judged on its
    replicates, by how many there are:
      * three or more: a one-sided 90% t-interval on the mean ("Confidence interval"): a change needs the
        near bound beyond the limit, no change the whole interval inside it;
      * two (13 of the 32 experiments with metabolites in mGrowthDB): both replicates beyond the limit on the
        same side, or both inside it ("Pairs agree"), since an interval on one degree of freedom left clear
        pairs such as -3.2 and -1.5 mM inconclusive;
      * one: inconclusive ("1 is inconclusive"); none exists today, and identical non-zero replicates (one
        series deposited twice) count as one.
    Anything else is `inconclusive`: NA with its own evidence, never an arc and never a 0. Across experiments,
    those that decide must say the same (else `conflict`); an inconclusive experiment whose mean lies beyond
    the limit on the other side is a conflict too, as are experiments none of which decides but whose means
    lie beyond the limit on both sides. One that does not contradict is named (`experiment_left_out`), and
    its mean still counts in the value, since the experiments left out are those with the smaller effects.
    The value is the mean of the experiments' means; `n`, `n_experiments`, the spread and the test are over
    them; amounts more than twofold apart are the caution `amounts_differ`. Three rules came first and were
    dropped over four rounds of review: the mean plus and minus one standard deviation (a spread, not an
    inference: it called changes more readily with two replicates than with ten), every replicate beyond
    the limit at any n (one failed sample erased E. coli LF82's -8 mM pyruvate uptake, seen in four of five
    replicates), and the interval at every n (too strict for pairs, and a single replicate decided alone
    gave an arc for half of true zeros). The 0.2 mM limit stays; compounds at another scale can have their
    own (thiamine=0.01, which also finds mGrowthDB's "thiamine(1+)"). Advanced settings can judge by the
    mean alone, as 0.1.0 did. The multiple-testing family is the tests that make the arcs; the p-value is
    reported, the rule above decides.
21. **Merging to genus applies the detection limit**, and taxa that disagree make the genus cell
    inconclusive (`conflict`), instead of a median that can sit below the limit or hide a consumer.
22. **One value medium per search, made stable and visible.** Karoline, 2026-10-06: "keep one, make it
    stable". Only taxa whose cultures give values vote; ties go to replicates with values, then the name;
    the page names the taxa with more data in another medium (their values depend on what was searched with
    them), a narrow win, and the taxa with no phase; every matrix CSV names the value medium in its first
    header cell. Choosing per taxon would make a matrix mix media, which a CRM cannot use.
23. **Cells say why they have no value.** New evidence states: `inconclusive` (decision 20), `no_phase`
    (assayed, but the cultures gave no end of exponential growth or no stationary phase), and
    `seen_elsewhere` (0 in the value medium, a change that way in another). 0.1.0 wrote the first two as
    `not_assayed` or 0, and the third as `below_limit` while its README promised NA.
24. **Exponential growth ends at the 90% rule, or earlier where the growth rate drops.** Karoline,
    2026-10-06 ("Earlier of 90% and rate"). The 90% rule ends growth late on two kinds of real curve: a
    plateau that keeps creeping up (E. coli LF82, study SMGDB00000009: growth ends near 10 h, but the count
    rises 20% more by 168 h, so 90% of the final maximum comes at 84 to 108 h) and a late rise that is not
    growth (B. hydrogenotrophica, study SMGDB00000004: qPCR rises again at 30 and 48 h while OD falls, DNA
    from lysing cells). The growth-rate rule ends growth at the first sample after the fastest growth from
    which the specific rate stays below a tenth of its maximum over two consecutive intervals; the boundary
    is the earlier of the two. On the curves checked it ends E. coli at 8 h in all five replicates and B.
    hydrogenotrophica at 15 to 17 h. A fourth round of the review counted all of mGrowthDB: the rate rule
    places 19 of 83 replicate boundaries, most of them one sample (4 to 8 h) before the 90% rule (B.
    thetaiotaomicron, R. intestinalis, B. fragilis, F. prausnitzii), and E. coli 76 to 100 h before it; on
    coarse curves (24 h apart) it finds nothing and the 90% rule decides. Values whose boundary the rate rule
    moved by more than one sample carry `growth_rate_boundary`. A first attempt read the 90% rule on a running
    median of three; the second round of the review showed it moved E. coli later still (84 to 120 h),
    did not remove a two-point late rise, and clipped a real one-sample peak (R. intestinalis ri2_B), so it
    was dropped. The first growth curve that gives a boundary is used; relative 16S is no growth curve;
    boundaries on fewer than three samples (`coarse_sampling`) or further apart across replicates than a
    sampling interval (`boundaries_differ`) are flagged, and so is a boundary the rate rule placed
    (`growth_rate_boundary`). A third round of the review showed the rate rule ending E. coli's growth at 8 h
    while it fermented its glucose from 8 to 12 h (after pyruvate ran out, cells +15%); Karoline chose to keep
    the boundary and flag it: where the rate rule moved the boundary, a compound whose replicates kept
    changing beyond the limit in the first interval after it carries `still_changing`, on its exponential
    value (which lacks that change, and which the CRM takes) and on its stationary value when that change is
    at least a quarter of it. On E. coli LF82 this marks the glucose taken up from 8 to 12 h, so the
    exponential phase's -3.1 mM, against -8.2 mM over 0 to 12 h, does not pass as its glucose uptake. Not chosen: a whole-run fallback for cultures
    without a boundary (the unblanked OD of study SMGDB00000010, which starts near 0.7, still reads five
    taxa that made 4 to 6 mM butyrate as not grown; they are `no_phase`, and a time window gives them
    values), and flags for failed samples and transient peaks.
25. **Media differ by amount too, and by more phrasings.** An alteration keeps its amount (0.1% and 0.75%
    linoleic acid, study SMGDB00000014; 0 to 5 mg/L pantothenate, study SMGDB00000019, formerly one medium).
26. **CRM parameters follow miaSim's equations.** Karoline, 2026-10-06: "rebuild from miaSim". Read from
    miaSim 1.18's `consumerResourceModel`: a positive entry of E is a yield (uptake is `R/(R+K)` per unit
    abundance whatever E is), a negative one the by-product per unit of growth. miaSim has no uptake rate,
    so the unit of abundance sets how fast a taxon eats (a second round of the review found a first version,
    with E as a yield in cells/mL, emptying the medium within 0.1 h, with the growth rates cancelling out).
    So each taxon gets a unit of its own, `crm_scale()` = `dx * C / (mu * S)` growth curve units (`S` the
    sum of the squared uptakes), in which E is each consumed resource's share of the uptake and
    `-produced * C / S` on each by-product: alone, a taxon grows at its measured rate while saturated, gains
    its measured biomass and makes its measured by-products in proportion to its uptake, and the simulation
    does not depend on the unit (a third round checked this, and the community's mass balance). Shares,
    not `1/n`, so a trace substrate does not slow a taxon down. The payload (now `foodnet.crm/v1`) carries each
    taxon's biomass change over the phase, keyed by technique and unit; how uptake splits between resources
    follows the Monod constants, which foodnet does not measure, and `crm_backcheck()` simulates each taxon
    alone to test a choice, including how long it takes to grow.
    `as_miasim()` requires starting abundances and Monod constants (miaSim would draw them at random) and
    refuses pooled media and the stationary phase unless allowed. 0.1.0's uptake shares remain as
    `scale = "shares"`.
    A final round found E. coli LF82's fitted rate (easylinear, 0.375 per hour) below the mean rate its own
    curves show over the phase (0.65): the fitting window took in the plateau of a fast, sparsely sampled
    curve, and a simulation grew it at half speed. The payload carries each taxon's phase mean rate, the page
    and R warn when it exceeds the fitted rate, and `growth = "phase_floor"` uses it there; the default stays
    the fitted rate (grownet's method, Karoline's choice of 2026-10-04).
28. **A culture without an end of exponential growth gives its change over the whole run.** Karoline,
    2026-10-06 ("the physiologist proposal sounds like a good answer"; "In the phase column, marked"). The
    growth curve only places the boundary between the phases; it is no claim that cultures metabolize only
    while growing. A culture with no boundary (one read as not grown, which an optical density without its
    blank does to five Db-MM taxa of study SMGDB00000010 that made 4 to 6 mM butyrate, or one without any
    growth curve) now gives the change from its first to its last metabolite sample, in the column of the
    phase asked for (with Both, the exponential column; the stationary one says `no_phase`), with the
    evidence and caution `whole_run`, arcs in phase `whole_run`, and a dashed frame in the image. Where any
    culture of a cell has a phase value, the whole-run changes of the others are left out of it and named.
    Whole-run values do not vote for the value medium. The CRM takes them, and its caveats name these taxa,
    since stationary uptake is in their values; no growth rate is taken from their curves (give one).
    A tenth review round showed the rule would also give values to cultures that truly did not grow (a dead
    inoculum, evaporation or abiotic drift as a "consumed" arc), since the 1.5-fold rule cannot tell them
    from unblanked OD. Karoline, 2026-10-06 ("Absolute OD rise"): an OD curve that rose by 0.1 or more still
    counts as grown and gives its whole-run change; a culture that grew by neither rule gives no value
    (`not_grown`, no arc, named in a warning). An eleventh round (the physiologist and the critic, separately)
    showed that this threw away A. soehngenii, whose unblanked OD did not rise but whose replicates turned
    glucose and lactate into butyrate. Karoline, 2026-10-06 ("Activity counts"; "Yes, 0.05 OD"): a culture
    whose growth curve shows no growth still gives its whole-run change, marked `growth_unclear`, when a
    compound was used up beyond the limit and another made beyond it; `not_grown` is kept for cultures whose
    compounds did not move so. An OD curve that grew by the fold rule must also rise by 0.05, since near the
    blank a fold is noise; a culture without any growth curve is marked `growth_unknown`. The whole-run
    taxa are named in the first header cell of every matrix, and their names stay the same in every file.
    A twelfth round showed that test too lenient (noise passed it in up to 85% of simulated non-growing
    cultures, and evaporation read as production and consumption) and the 0.05 floor too strict (Variovorax's
    low inoculum, 0.002 to 0.05 OD). Karoline, 2026-10-06: "raise the threshold on change in metabolite
    concentrations for non-growing organisms so it exceeds what would be expected based on evaporation", and
    "Smoothed fold, no floor". So a culture that did not grow counts as active only when a compound was used
    up (to below the limit, or falling as a trend) and another made (rising as a trend: Spearman's rank
    correlation with time of 0.8 or more), each by more than the larger of the limit and 10% of its level
    (what evaporation could account for); its whole-run values must clear the same threshold
    (`within_evaporation`), and any whole-run value must exceed twice its series' own scatter
    (`within_scatter`). For OD, the 1.5-fold test reads the maximum of a running median of three.
    A thirteenth round refined both: volatile compounds (ethanol, methanol, acetone, gases) never count as
    used up, since they can leave as vapor; the evaporation share is a setting (10% by default, more for open
    plates) and explains rises only; the scatter test applies to every value, phase or whole run, measured
    within the window the change spans (the spread of its sample-to-sample steps, so a sharp depletion is no
    scatter) and only to changes beyond the limit, so a flat series stays a measured 0; and a value flagged
    by either test keeps its vote (a fourteenth round showed that dropping it selects replicates by their
    values), so it can stop a call, but cannot make one: a change needs at least two unflagged replicates
    that show it themselves (a fifteenth round found a Bacteroides pair, +0.20 and a flagged +0.21 mM, read as
    butyrate production).
27. **Hardening.** The R listener takes parameters only with a one-time secret it writes to a file only
    this user can read, since base R cannot bind a port to 127.0.0.1; a result whose records could not all
    be read says `incomplete` in its files and the command line exits 3; the page never serves another
    search's data for a job that is gone.

## Known limits

- **The validation is not independent.** The reproduction of a hand-checked reference matrix (106 of 108
  cells) used ids, a time window and a second window chosen to reproduce it, and the comparison with a
  published figure compared directions of dominant fermentation products, which agree often by chance. Neither
  tests the default exponential-phase rule on held-out studies.
- **False arcs.** Each experiment is decided at one-sided 90% confidence, so a compound at the limit is
  called beyond it about one time in ten; in simulation a true zero gave an arc in at most 4% of cells with
  two or more replicates. Over the 80 or so values of a search of all of mGrowthDB a few arcs are expected to
  be false; the q-values are reported, not used to decide.
- **Pairs err more readily than triplicates.** A pair decides when both replicates lie beyond the limit, with
  no estimate of noise: on noisy data (a standard deviation of 1 mM) about a third of pair-decided arcs on a
  true zero are false, against about one in seven for three replicates, and a real change of 0.3 mM is read
  as no change about one time in sixteen (one in seventy for three). Such values, changes and zeros alike,
  carry `pair_decided`; a matrix 0 with that caution in cautions.csv rests on two replicates.
  Identical replicates count as one only within twice the limit, so 0.39 mM three times is inconclusive and
  0.41 mM three times is a change.
- **The exponential phase can miss a slow-down.** Where the growth rate ends growth early (E. coli LF82),
  part of a compound's use falls in the stationary phase; the values say `still_changing`, the page, the
  zips, crm.json and the R package warn, and a time window over both phases gives the whole change.
- **The scatter test needs four samples in a window.** A window of fewer (Db-MM's 0 to 24 h exponential phase
  holds two) gets no scatter check; such values carry `coarse_sampling` or `pair_decided`. Volatile compounds
  are recognized by name or ChEBI id from a fixed list.
- **A net change hides what was made and used again** within a phase (formate in E. coli), and a single
  failed sample at a phase end becomes the value; neither is flagged.

## Open decisions

None at the moment.
