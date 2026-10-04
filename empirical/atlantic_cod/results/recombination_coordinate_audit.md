# Atlantic cod linkage coordinate audit

    Source Data Fig. 1 positions are interpreted as gadMor2 linkage-group base-pair coordinates. The frozen Atlantic cod analysis also uses the gadMor2 cod-phylogenomics coordinate set recorded in `data/atlantic_cod/metadata/region_manifest.tsv` and `data/atlantic_cod/processed/cod_window_trees.tsv`.

    ## Frozen inversion coordinates

    - LG01: 9114741-26192386
- LG02: 18489307-24048607
- LG07: 13610591-23019113
- LG12: 638100-14327837

    ## Source coordinate ranges

    - LG01: 1..26242489 max observed SNP position; 298 linkage SNPs
- LG02: 1..24100282 max observed SNP position; 230 linkage SNPs
- LG07: 1..23066726 max observed SNP position; 249 linkage SNPs
- LG12: 1..13681347 max observed SNP position; 197 linkage SNPs

    ## Window inventory

    - LG01: 110 frozen 250-kb windows
- LG02: 92 frozen 250-kb windows
- LG07: 119 frozen 250-kb windows
- LG12: 105 frozen 250-kb windows

    Coordinate compatibility: compatible. Source SNP positions, frozen inversion boundaries, and frozen SNAPP-window starts/ends are all linkage-group bp coordinates in the same gadMor2 coordinate frame. Stage 9 joins linkage, topology, and time tracks only by LG/start/end physical coordinates. Alternate inversion limits from other workflows were not substituted.

    Frozen-manifest hash check before/after Stage 9: {"after": {"stage4a": "23e21f056f94fa5e492156fbe9033423eb536053c86def00f4814cd46e246b72", "stage4b": "6d8eb64e16947b3207a48390c615b7e4c78e45bc73a6fa940a560be0ab77884b", "stage5": "5ae6f3274f0accd9c6aedc1a8b512901a5c67209bb4024987d5192acabc445db", "stage6": "be05abce694a782fc685fc1c8aa9004e106299474b814f86cb9f4fac8bd44e28", "stage7": "8385afddf18cb7570317269f89ae4e5c4f5ae4b2f7fa944bef256208f9de325b", "stage8": "dc0529677777fcc4c4063fcdd35250b55cb72e89cffd690fd9aeaf2e7bfb0970"}, "before": {"stage4a": "23e21f056f94fa5e492156fbe9033423eb536053c86def00f4814cd46e246b72", "stage4b": "6d8eb64e16947b3207a48390c615b7e4c78e45bc73a6fa940a560be0ab77884b", "stage5": "5ae6f3274f0accd9c6aedc1a8b512901a5c67209bb4024987d5192acabc445db", "stage6": "be05abce694a782fc685fc1c8aa9004e106299474b814f86cb9f4fac8bd44e28", "stage7": "8385afddf18cb7570317269f89ae4e5c4f5ae4b2f7fa944bef256208f9de325b", "stage8": "dc0529677777fcc4c4063fcdd35250b55cb72e89cffd690fd9aeaf2e7bfb0970"}}.
