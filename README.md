## Comparison of equally valid trees

An analysis of whether multiple equally valid phylogenetic trees produced by
distance-based methods (eg UPGMA) are actually close with negligible
differences, or actually very different. Intended as a continuation of
the [2022 paper by Segura-Alabart et al. in *Brief Bioinf*](https://academic.oup.com/bib/article/23/5/bbac312/6652780).

#### General idea

The overall approach is to count the pairs of taxa that do not change the
distance, *ie* a number of edges/nodes on the path between them. Intuitively,
the less such pairs of taxa, the more similar tree topologies are. This can be
calculated either in a pairwise manner (to get a distance matrix for the tree
set), or for a dataset as a whole (to get an overall feeling of whether it
changes significantly or not).

#### Test data
Small test dataset for debug is available in `./test_data` and consists of
three files:

`test_trees.nwk` contains trees from fig. 1A of the Segura-Alabart paper (lines
1 and 3) and a rerooting of first tree (line 2).

`Thai_2007.nwk` and `Karatas_2019.nwk` are two small-ish tree sets from
the same paper. (106 and 108 trees, respectively)

#### License, citation and stuff

All code is available under MIT license, although test datasets may be subject
to other licenses (see original papers).

There is no academic publication for now, but there will be one at some point,
so please keep citation in mind.