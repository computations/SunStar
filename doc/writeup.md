---
author: Ben Bettisworth
title: "SunStar: an Implementation of the Generalized STAR Method"
date: 2017-03-13
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

[^SunStar_name]: The name is a slight pun on GSTAR. The sun that earth orbits
is a g-class star, therefore SunStar.

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

A tree is said to be in _Newick Notation_ when it is represented as a string
according to the standard found
[here](http://evolution.genetics.washington.edu/phylip/newicktree.html)
[@newick_format]

We will use the notation $F(n) = \OO(f(n))$ to denote the order of a function
$F$. Normally, $F$ is not reported, but instead referred to by this notation.
By convention, this will roughly the smallest $f$ that satisfies the
relationship.

Background
-------------------------------------------------------------------------------

Phylogenetics is primarily concerned with the inference of species trees, but
in general we can only observe gene trees. Since organisms are made of genes,
it follows that the gene trees should be related to the species tree. This
naturally leads to the problem of combining gene trees together to make a
species tree. This is the goal of STAR .

Newick
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
-   `t` `--trials`: The number of trials using random weights. If this flag
is not present, then a default schedule will be used.
-   `-l` `--logfile`: Filename to log the sequences that are generated by the
`-t` option
-   `-o` `--outgroup`: The name of the outgroup taxa. This taxa must be present
on all of the trees.

###Default Schedule

There is a default schedule for the weights. This schedule starts by assigning
1 to every edge, and adjusts leaf edges to be ultrametric. Then, all possible
combinations of weights with at least one weight being 1 are computed[^binary].
By the results from generalized STAR, if there is no error then the tree should
come out the same, despite the 0 weights on the edges. But, We believe that it
is likely that introducing zeros into the system will maximize the exposer of
any errors, and therefore will cause instabilities.

[^binary]: This should be thought of as counting in binary from 1 to the
height of the tree.

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

![An example tree \label{fig:dist_tree}](./figs/distances.png){width=50%}

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

Neighbor Joining
-------------------------------------------------------------------------------

###Specification

###Complexity

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
change the overall $\OO(mn^3)$.


Generalized STAR
-------------------------------------------------------------------------------

Implementation
===============================================================================

\texttt{tree.h}
-------------------------------------------------------------------------------

\texttt{newick.h}
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
