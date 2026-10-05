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
4. **The detection limit is 0.2 mM**, an advanced setting, applied to the mean change (Figure 3c of the
   community control paper, where the data break).
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
    values while reproducing Figure 3c: ids set the scope (which data), a medium name the value medium within
    it, and with ids alone the majority rule runs inside the scope, its other media giving presence.
14. **A second time window for named metabolites.** Karoline, 2026-10-05: "to include the different time
    window for trehalose ... I'd like to avoid a table where each metabolite gets 1 time window. perhaps 2
    time windows, with the option to add metabolites by name". The named metabolites take the second window
    instead of the phases or the main window; an empty end is each culture's last sample. Names match the
    metabolite or the name its study recorded it under. Columns, the image, the README and the report say
    which interval each value covers. With 0 to 48 h and trehalose to the end, every trehalose value of
    Figure 3c is reproduced.
15. **Studies are found through the species list**, not mGrowthDB's search, which matches per-strain
    measurements only and so misses a monoculture measured at the culture level (study SMGDB00000009).

## Open decisions

1. **The 90% boundary** (decision 2): the share, and linear against log scale.
2. **"Both" for the CRM**: the exponential phase is used. Alternatives: the stationary phase, or both
   phases as separate parameter sets.
3. **Zeros that look like missing values.** In study SMGDB00000009, B. fragilis formate reads 0.0 at 4 h in
   one replicate and at 48 h in another, between values near 10 mM. Read as measurements, the 48 h one
   turns a production of about 2 mM into a mean change near zero. food**net** uses the values as recorded;
   whether to treat an isolated zero as missing is open, and probably a question for mGrowthDB.
4. **Single replicates** are shown, flagged. Whether to hide them by default, as grow**net** once did.
5. **The growth-rate fallback** takes any monoculture of the taxon in the value medium. Whether it should
   also require the same study.
