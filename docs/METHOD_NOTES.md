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
    caution ("OK"); the growth-rate fallback takes another monoculture of the taxon in the same medium ("keep
    as is, but with the stringent medium matching"), which decision 18 now makes strict.

20. **The experiment is the unit, and a change must clear the limit across its spread.** Karoline,
    2026-10-06, after a review of 0.1.0 showed a study with ten replicates outvoting one with two, and
    +0.84 +/- 1.04 mM on a 27 mM background drawn as production: "experiment level + inconclusive". A cell's
    mean, standard deviation and test are over its experiments' means (over the replicates when one
    experiment gives it); `n` stays the replicates and `n_experiments` is new. A change needs the mean plus
    and minus one standard deviation beyond the limit on one side; no change needs it inside the limit;
    anything else, and experiments that disagree, is `inconclusive`: NA with its own evidence, never an arc
    and never a 0. The 0.2 mM limit stays; compounds at another scale can have their own (thiamine=0.01).
    Advanced settings can judge by the mean alone, as 0.1.0 did. Identical replicates give no test
    (`no_variance`), and the multiple-testing family is the tests that make the arcs.
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
24. **A more robust phase boundary.** Karoline, 2026-10-06: "robust boundary". The 90% rule is read on the
    growth curve smoothed by a running median of three (Tukey's rule at the end, the start as measured); the
    first growth curve that gives a boundary is used; relative 16S is no growth curve; boundaries on fewer
    than three samples (`coarse_sampling`) or further apart across replicates than a sampling interval
    (`boundaries_differ`) are flagged. Not chosen: falling back to the whole run for cultures without a
    boundary, flags for failed samples and transient peaks, and a boundary from the growth rate. Checked on
    the real curves: smoothing moves single-point noise and a lone late jump, but not a plateau that keeps
    creeping up (E. coli LF82, study SMGDB00000009, +20% from 12 to 168 h, whose boundary lands at 84 to
    120 h); those values carry `boundaries_differ`. See the open decision below.
25. **Media differ by amount too, and by more phrasings.** An alteration keeps its amount (0.1% and 0.75%
    linoleic acid, study SMGDB00000014; 0 to 5 mg/L pantothenate, study SMGDB00000019, formerly one medium).
26. **CRM parameters follow miaSim's equations.** Karoline, 2026-10-06: "rebuild from miaSim". Read from
    miaSim 1.18's `consumerResourceModel`: a positive entry of E is a yield (uptake is `R/(R+K)` per unit
    abundance whatever E is), a negative one the by-product per unit of growth. miaSim has no uptake rate,
    so the unit of abundance sets how fast a taxon eats (a second round of the review found a first version,
    with E as a yield in cells/mL, emptying the medium within 0.1 h, with the growth rates cancelling out).
    So each taxon gets a unit of its own, `crm_scale()` = `n * dx / (mu * C)` growth curve units, in which E
    is `1/n` on each consumed resource and `-n * produced / C` on each by-product: alone, a taxon grows at
    its measured rate, gains its measured biomass and makes its measured by-products in proportion to its
    uptake, and the simulation does not depend on the unit. The payload (now `foodnet.crm/v1`) carries each
    taxon's biomass change over the phase, keyed by technique and unit; how uptake splits between resources
    follows the Monod constants, which foodnet does not measure, and `crm_backcheck()` simulates each taxon
    alone to test a choice, including how long it takes to grow.
    `as_miasim()` requires starting abundances and Monod constants (miaSim would draw them at random) and
    refuses pooled media and the stationary phase unless allowed. 0.1.0's uptake shares remain as
    `scale = "shares"`.
27. **Hardening.** The R listener takes parameters only with a one-time secret it writes to a file only
    this user can read, since base R cannot bind a port to 127.0.0.1; a result whose records could not all
    be read says `incomplete` in its files and the command line exits 3; the page never serves another
    search's data for a job that is gone.

## Known limits

- **The validation is not independent.** The reproduction of a hand-checked reference matrix (106 of 108
  cells) used ids, a time window and a second window chosen to reproduce it, and the comparison with a
  published figure compared directions of dominant fermentation products, which agree often by chance. Neither
  tests the default exponential-phase rule on held-out studies.
- **A net change hides what was made and used again** within a phase (formate in E. coli), and a single
  failed sample at a phase end becomes the value; neither is flagged.

## Open decisions

1. **The end of exponential growth on a creeping plateau.** Smoothing (decision 24) does not place it for
   E. coli LF82; a boundary from the growth rate (the first interval whose rate falls below a tenth of the
   maximum) places it at 8 h in all five replicates, but on noisy qPCR curves (study SMGDB00000004) it ends
   growth early (3 to 8 h). Karoline to decide whether to switch, combine, or keep the 90% rule with its
   caution.
