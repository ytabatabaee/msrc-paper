# msrc-sim mechanistic audit

The pilot uses `ytabatabaee/msrc-sim` commit
`f5157d2c64068eb32e0eabfe9bec629c28a404f9`, package version 0.8.8. The
repository is used in place; no simulator source is copied into `msrc-paper`.

The mechanistic configuration supplies an ultrametric Newick species tree,
branch lengths in generations, a root extension, and default or branch-specific
effective population sizes. Internal node names such as `A`, `B`, and `ROOT`
are population-branch identifiers. In the balanced rooted Newick tree used by
the pilot, each cherry has a daughter ancestral branch of length `g`. The
effective unrooted quartet internal length is therefore
`t_unrooted = 2g/(2 Ne) = g/Ne` under the simulator's diploid-copy scaling.
The ordinary MSC control uses this unrooted length in
`P(Q_S)=1-(2/3)exp(-t_unrooted)` and
`P(Q_A)=P(Q_O)=(1/3)exp(-t_unrooted)`.

Forward rearrangement origins are supplied by `rearrangement.origin_branch`
and integer `origin_time_from_branch_start`. `initial_copy_count` gives the
number of rearranged chromosomes at the origin, out of `2 Ne` copies. The
current pilot uses origins on the ancestral `ROOT` branch at different depths,
so the rearrangement can be inherited by both daughter clades without
conditioning on a desired terminal pattern. Descendant-branch origins are
supported by the input class but are treated as a separate diagnostic in later
work because they do not represent an ancestral origin shared by both clades.

The Wright–Fisher process updates a segregating count each generation with a
binomial draw over `2 Ne` chromosome copies. Neutral trajectories can be lost,
fixed, or remain polymorphic. At speciation, daughter populations are sampled
from the parent terminal frequency. The default split rule is unbiased; the
optional `speciation_split` section can apply supported daughter-frequency
rules, but it is not used in this neutral pilot. Thus persistence through
speciation is an output of the simulated history, not an imposed terminal
state.

For each realized history, the backward genealogy simulator starts one sampled
lineage per taxon with the sampled terminal arrangement. Lineages switch
arrangement states according to the effective cross-arrangement rate and the
realized population frequency. Same-state lineages in the same population
branch coalesce with rates determined by the state frequency and effective
population size. The public replicate API records quartet counts, terminal
patterns, root-end status, and per-replicate diagnostics. The public CLI writes
these records without requiring any analysis-side population-genetic model.

The recombination parameter is `baseline_rate` together with
`effective_cross_arrangement_fraction`; the backward implementation uses their
product `rho`. A fraction of 1.0 is the no-suppression control and 0.1 is the
strong-suppression pilot condition. This is the implemented parameterization;
the pilot does not introduce a `suppression = 1-m` parameter.

The package also exposes `simulate_msc_genealogy`, which is used for matched
ordinary MSC controls on the same species tree. For a four-taxon balanced tree
with internal length `t=g/(2 Ne)`, the control expectation is
`P(Q_S)=1-(2/3)exp(-t)` and `P(Q_A)=P(Q_O)=(1/3)exp(-t)`. The control runner
compares observed frequencies with this expectation and Monte Carlo standard
errors.

The mechanistic-to-conditional comparison is only local. If the frequency is
approximately constant at `p` and time is expressed in the same units, the
mechanistic switching rates are approximately `m01=rho*p` and
`m10=rho*(1-p)`, while state-specific coalescence rates are approximately
`1/[2 Ne (1-p)]` and `1/[2 Ne p]` before finite-copy safeguards. In the
mechanistic simulation `p` varies by branch, generation, and replicate, so no
single conditional `(m01,m10,lambda0,lambda1)` represents an unconditional
mechanistic cell. Experiment 1 therefore remains a conditional benchmark,
not a mechanistic reduction.
