---
author: Ben Bettisworth, University of Alaska Fairbanks
title: "SunStar: An Implementation of the Generalized STAR Method \\ DRAFT 1"
subtitle: DRAFT 1
date: 2017-03-13
documentclass: article
bibliography: bib.yaml
numbersections: true
geometry: margin=1in
toc: true
abstract: STAR[@liu09] is a method of computing species trees from gene trees. Later,
    STAR was proven to be statistically consistent given a few conditions
    [@rhodes_star]. Using these conditions, it is possible to detect errors in
    the species tree inference process, which will produce instabilities in the
    tree resulting from STAR. We have developed a piece of software that does
    this called \texttt{SunStar}.
...

\newcommand{\OO}{\mathcal{O}}
\newcommand{\SunStar}{\texttt{SunStar}}

Introduction
===============================================================================

\SunStar[^SunStar_name] is a program designed to infer the support for
a species tree that has been inferred from gene trees via the STAR method.

[^SunStar_name]: The name is a slight pun on G-STAR, the original name for this
project. The sun that earth orbits is a g-class star, therefore SunStar.

Background
-------------------------------------------------------------------------------

###Problem Description

Phylogenetics is primarily concerned with the inference of species trees.
Unfortunately, we cannot directly infer species trees. Instead, phylogenetic
methods are limited to traits we can observe. Historically, these traits were
features like morphology of bones and limbs. With the advent of molecular
sequencing, observation of traits using DNA sequences became common, and with
those observations came the early mathematical methods of phylogenetics.

These methods were primarily using gene trees as proxies for species trees. Due
to incomplete lineage sorting (ILS) [@pamilo88] this is not sufficient for
inferring species trees. Informally, this means that genes might diverge
differently than the species diverge. Therefore, just using a single gene tree
as a proxy for species doesn't work. More advanced methods of inferring
species trees get around the ILS method by using information from multiple
genes.

###Prior Work

Previous software that does gene tree summary includes ASTRAL [@astral], NJst
[@njst], STAR[@liu09], BEST[@best] and ASTRID[@astrid]. All of these programs
will compute species trees from existing gene trees. ASTRAL is the exception in
this group, in that it does not compute a distance table from the set of gene
trees passed to it. For the rest of the software, they all compute a version of
a distance table.

In particular, STAR and ASTRID are a very similar method to the method used in
\SunStar, with one important difference, they do not attempt to infer the
support for the tree reported. \SunStar does, using the results from
generalized STAR [@rhodes_star].

Notation, Conventions and Definitions
-------------------------------------------------------------------------------

A _tree_ is a set of vertices(or nodes) and a set of edges that connect them.
A vertex with no children is called a _leaf_, plural _leaves_. A tree can have
a special node that is designated the _root_. The root of a tree is often
labeled as $\rho$. A node is called an _interior_ node if it is not a leaf.

Edges can have have tags which are numbers. These are typically called
_weights_ or _lengths_, and are often a measure of a notion of distance.

A tree is _trivalent_  when all vertices of the tree have degree 3 or 1. A
rooted trivalent tree is a tree that is allowed to have a single vertex with
degree 2. All trees discussed here will be trivalent unless otherwise noted.

A _gene tree_ is a rooted or unrooted tree that has been inferred from gene
sequences, and relates the divergence of those genes. A _species tree_ is a tree
that relates the divergence of species.

An _ultrametric_ tree is a tree where all the distance to the root vertex are
the same.

We will use the notation $F(n) = \OO(f(n))$ to denote the order of a function
$F$. Normally, $F$ is not reported, but instead referred to by this notation.
By convention, this will roughly the smallest $f$ that satisfies the
relationship.

Newick Notation
-------------------------------------------------------------------------------

Newick Notation is a format of specifying trees. An example of a Newick tree is

```
((a:.1,b:.3):1.0, (c:.43,d:.5):1.0);
```

The '`:`' signifies a weight on an edge leading up from the node towards the
root, and is optional.

One advantage of Newick notation is that the grammar[^newick_ref] for this
language is quite simple, and requires only a few productions. This makes it
easy to write a parser to recognize a Newick string. Unfortunately there are
quite a few disadvantages, one of which is the non-uniqueness of a tree and
a string. Specifically, there are many strings which represent the same tree.
The problem is due, largely to the ordering of subtrees. Note that these two
trees are the same

```
(a,b); == (b,a);
```

This causes problems when attempting to compare trees, which is already a hard
task. Fortunately, we can side step these issues by

-   Enforcing an outgroup (to insure a consistent root), and
-   Enforcing an order on the taxa and subtrees.

By doing this, we can compare trees by comparing their Newick strings.

One additional problem is the lack of standardized grammar. Some
implementations allow for taxa labels to start with a number. Some
implementations require the semicolon at the end, others do not. Its difficult
to deal with this issue, so instead I documented the grammar I used in section
\ref{newick-imp}.

[^newick_ref]: See section \ref{newick-imp}

Neighbor Joining
-------------------------------------------------------------------------------

STAR(estimating Species Trees using Average Ranks)
-------------------------------------------------------------------------------

STAR infers a species tree by calculating a distance table of ranks[^1], and
then a tree building algorithm to generate a new tree [@liu09]. The algorithm
is as follows: Given a set of gene trees, do the following for each tree.

-   Set all the edge weights
    -   If the edge is not connected to a leaf, set the weight to 1.
    -   If the edge is connected to a root, set the weight to .5
    -   Otherwise set it to the height of the tree minus the distance to the
    root.
-   Calculate pairwise distance, put them into a table.

Once all the tables have been calculated, do an entry wise sum, and divide each
entry by the number of gene trees. Informally, think of this as taking the
average of the distance table. Once we have the "average" distance table, run
a tree building algorithm on it, such as Neighbor Joining.

[^1]: Here the term rank is used to mean edge length or weight, and will be
  referred to as edge length in the future.

Generalized STAR
-------------------------------------------------------------------------------

Generalized star [@rhodes_star] is version of STAR where the edge weight
schedule is allowed to vary. In this version, edge weights of an equal distance
from the root are given the same weight. Weights need only be non-negative
except for one. Additionally, the edges connected to leaves are adjusted to
make the tree ultrametric, as in STAR. Other than weight assignment, algorithm
proceeds the same.

Additionally, we can attempt to infer the stability of the species tree
produced by STAR by varying the weights and reporting the different trees that
come out of the tree creation step. This is the primary goal of \SunStar.

\texttt{SunStar}
-------------------------------------------------------------------------------

\SunStar currently has only a few options. Here they are:

-   `-f` `--filename`: Filename for a set of Newick trees. These trees must all
together be rooted or unrooted. If they are all unrooted trees, then the `-o`
is required.
-   `t` `--trials`: The number of trials using randomly genereated schedule of
weights(Section \ref{randomized-schedule}). If this flag is not present, then
a default schedule (Section \ref{default-schedule})will be used.
-   `-l` `--logfile`: Filename to log the sequences that are generated by the
`-t` option
-   `-o` `--outgroup`: The name of the outgroup taxa. This taxa must be present
on all of the trees.



Algorithms
===============================================================================

Finding the distance between two taxa
-------------------------------------------------------------------------------

###Specification

In order to calculate the distance table, we have to calculate the distance
between any two taxa. Fortunately we are working with trees, so the path to any
two nodes is unique. Therefore, if we find a path, we have found the only path.
So, the strategy to find a path between to nodes is to calculate the list of
parents for each node. Once we have this list, we walk the list, starting from
the parent side, until we find a difference. Then, the remaining nodes in the
list, plus the node right before the divergence, is the path between the two
nodes.

For example, suppose we the tree in figure \ref{fig:dist_tree}. Then the
sequences of parents associated with that tree are
$$s_1 : \rho, n_1, t_1$$
$$s_2 : \rho, n_1, t_2$$
$$s_3 : \rho, n_2, t_3$$
$$s_4 : \rho, n_2, t_4.$$
Notice that for taxa $t_1$ and $t_2$, the sequences differ only at $t_1$ and
$t_2$. So backing up a step, we have the path $t_1, n_1, t_2$. For the taxa
$t_1$ and $t_4$, the node before the difference is $\rho$. So, the path between
$t_1$ and $t_4$ is $t_1, n_1, \rho, n_2, t_4$.

![An example tree \label{fig:dist_tree}](./figs/distances.png){width=30%}

###Complexity

The complexity of finding the distance between two nodes for this is $\OO(n)$,
where $n$ is the number of taxa on the tree. First, consider the steps involved
with this algorithm. They are

-   Make a list from taxa to root,
-   Walk along that list,
-   Create new path,
-   Sum edge weights over path.

Note that the number of nodes on a rooted tree with $n$ taxa is $2n-3$. Since
a tree is acyclic, each node in the path is visited at most once. Therefore the
size of the lists from taxa to root is at most $\OO(n)$. Since the lists are
bounded by $n$, walking along them is also bounded by $n$. Creating a new
path can only be as bad as concatenating two lists together. Since the two
lists are of $\OO(n)$, concatenating them is only as bad as $2\OO(n) = \OO(n)$.
Finally, summing over the final list of nodes requires iterating through the
final path. Since the final path is $\OO(n)$, this step is also, and therefore
this method is as well.

STAR
-------------------------------------------------------------------------------

###Specification

The algorithm requires a set of rooted trees $|T| = m$. STAR involves three
major steps:

-   Calculate the distance table for each tree in the set
-   Average the distance tables
-   Run Neighbor Joining on the new table.

Calculating the distance table requires finding the pairwise distance between
all the taxa. To do this, we create a matrix with the taxa labeling the rows
and columns. For each entry in the matrix, we find the distance between two
taxa using the algorithm specified above. We need to do this for each tree in
$T$.

We now have a distance matrix that records the pairwise distance for each taxa
in the tree, for each tree in $T$. To average, we simply add each matrix, and
then divide by $m$. Once we have this new distance matrix, we pass it to
Neighbor Joining and get a new tree.

###Complexity

The complexity of STAR is $\OO(mn^3)$. Calculating the distance matrix requires
a call to a $\OO(n)$ algorithm, and there are $n^2$ elements in the matrix.
Therefore $n^2\OO(n) = \OO(n^3)$. We need compute this distance matrix for each
of the $m$ trees in $T$, so $m\OO(n^3) = \OO(mn^3)$.

Computing the average of the distance matrices involves a series of matrix
additions ($mn^2$) and a scalar matrix operation ($n^2$). These are less than
the time required to simply compute a single distance matrix, and so they don't
change the overall $\OO(mn^3)$. Similarly with Neighbor Joining, which has
a complexity of $\OO(n^3)$, but only needs to be performed once.

Generalized STAR
-------------------------------------------------------------------------------

As discussed above Generalized STAR works by varying the schedule of weights on
a set of trees. To implement this algorithm, we simply produce different
schedules and then run STAR on the tree sets with the different weights. There
are two strategies by which weight schedules are generated.

###Default Schedule

This schedule starts by assigning 1 to every edge, and adjusts leaf edges to be
ultrametric. Then, all possible combinations of weights with at least one
weight being 1 are computed[^binary]. By the results from generalized STAR, if
there is no error then the tree should come out the same, despite the 0 weights
on the edges. But, We believe that it is likely that introducing zeros into the
system will maximize the exposer of any errors, and therefore will cause
instabilities.

[^binary]: This should be thought of as counting in binary from 1 to the
height of the tree.

###Randomized Schedule

The other strategy to generate weight schedules is to generate them randomly.
This strategy works by generating a list of numbers from a uniform distribution
from 0 to 1. These numbers are then used as the schedule for a run of STAR.

Implementation
===============================================================================

\texttt{tree.h}
-------------------------------------------------------------------------------

\texttt{newick.h}{#newick-imp}
-------------------------------------------------------------------------------

\texttt{nj.h}
-------------------------------------------------------------------------------

\texttt{star.h}
-------------------------------------------------------------------------------

\texttt{gstar.h}
-------------------------------------------------------------------------------

Conclusion
===============================================================================

Future Work
===============================================================================

Optimizations
-------------------------------------------------------------------------------

Investigations
-------------------------------------------------------------------------------


Code
===============================================================================

The entire code base, including this write up, is located on
[GitHub](https://github.com/computations/SunStar)

References
===============================================================================
