SunStar
===============================================================================

Welcome to `SunStar`, an implementation of [generalized
STAR](http://www.dms.uaf.edu/~jrhodes/papers/STARandGeneralizations.pdf). There
is a detailed writeup of the algorithms an implementation in the `doc` folder.
There is a short version contained in this document as well

Purpose
-------------------------------------------------------------------------------

The purpose of `SunStar` is to infer the errors associated with inferring
species trees using the STAR method[1]. It does this by varying weight
schedules while running STAR.

Usage
-------------------------------------------------------------------------------

`SunStar` currently has only a few options. Here they are:

-   `-f` `--filename`: Filename for a set of Newick trees. These trees must all
together be rooted or unrooted. If they are all unrooted trees, then the `-o`
is required.
-   `-r` `--require-ratio`: Specifies the threshold ratio the inferred trees
must meat to be reported. Specifically the output of trees that occur less than
the specified ratio percent of the time will be suppressed.
-   `-t` `--trials`: The number of trials using randomly generated schedule of
weights. If this flag is not present, then a default schedule will be used.
-   `-l` `--logfile`: Filename to log the sequences that are generated by the
`-t` option
-   `-o` `--outgroup`: The name of the outgroup taxa. This taxa must be present
on all of the trees.
-   `-s` `--silent`: Silence progress bar output. Only output the final trees.

[1]: Estimating Species Phylogenies Using Coalescence Times among Sequences (Liu et. al. 2009).

