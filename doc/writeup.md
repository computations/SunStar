---
author: Ben Bettisworth, University of Alaska Fairbanks
title: "SunStar: An Implementation of the Generalized STAR Method \\ DRAFT 2"
subtitle: DRAFT 1
date: 2017-03-13
documentclass: article
bibliography: bib.yaml
numbersections: true
geometry: margin=1in
toc: true
abstract: STAR[@liu09] is a method of computing species trees from gene trees.
    Later, STAR was generalized and proven to be statistically consistent given
    a few conditions [@rhodes_star]. Using these conditions, it is possible to
    investigate robustness in the species tree inference process, which will
    produce instabilities in the tree resulting from STAR. We have developed
    a piece of software that does this called \texttt{SunStar}.
header-includes:
    -   \usepackage[linesnumbered,lined,ruled,vlined]{algorithm2e}
    -   \usepackage{nicefrac}
...

\newcommand{\OO}{\mathcal{O}}
\newcommand{\SunStar}{{\tt SunStar}}
\newcommand{\margalg}[1]{\IncMargin{0em}\begin{algorithm}[H] #1
\end{algorithm}\DecMargin{0em}}

Introduction
===============================================================================

\texttt{SunStar}[^SunStar_name] is a program designed to investigate the
support for a species tree that has been inferred from gene trees via the STAR
method.

[^SunStar_name]: The name is a slight pun on G-STAR, the original name for this
project. The sun that earth orbits is a g-class star, therefore \SunStar.

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
as a proxy for species can be misleading. More advanced methods of inferring
species trees get around the ILS difficulties by using information from multiple
genes.

An example of this is the multispecies coalescent model, which is a probabilistic
model of the ILS process. Using this model we can describe the way gene trees
form from species trees. The methods that are discussed in this paper are shown
to be statistically consistent under this model. In this context, statistically
constistent means that we can make the probability of inferring the correct
species tree equal to 1 if we make the sample large enough. Informally, given
enough perfect data, we always get the right tree.

###Prior Work

Previous software that infer species trees from from gene tree summaries
includes ASTRAL [@astral], NJst [@njst], STAR[@liu09]  and ASTRID[@astrid]. All
of these programs will compute species trees from existing gene trees. ASTRAL
is the exception in this group, in that it does not compute a distance table
from the set of gene trees passed to it. For the rest of the software, they all
compute a version of a distance table.

In particular, STAR and ASTRID are a very similar method to the method used in
\SunStar, with one important difference, they do not attempt to infer the
support for the tree reported.  \texttt{SunStar} does, using the results from
generalized STAR [@rhodes_star].

Notation, Conventions and Definitions
-------------------------------------------------------------------------------

A _tree_ is a set of vertices (or nodes) and a set of edges that connect them
with no cycles. A vertex with only one edge incident is called a _leaf_, plural
_leaves_. A tree can have a special node that is designated the _root_. The
root of a tree is often labeled as $\rho$. A node is called an _interior_ node
if it is not a leaf.

Edges can have have tags which are numbers. These are typically called
_weights_ or _lengths_, and are often a measure of a notion of distance.

A tree is _trivalent_  when all vertices of the tree have degree 3 or 1.
A rooted trivalent tree is a tree that is allowed to have a single vertex, the
root, with degree 2. All trees discussed here will be trivalent unless otherwise
noted.

A _gene tree_ is a rooted or unrooted tree that has been inferred from gene
sequences, and relates the divergence of those genes. A _species tree_ is a tree
that relates the divergence of species. An unrooted tree is often made rooted
by using an _outgroup_. This is a taxa that is intensionally distant from the
other taxa, ensuring that the root is connected directly to the outgroup.

An _ultrametric_ tree is a tree where all the distance to the root vertex to
the leaves are the same.

We will use the notation $F(n) = \OO(f(n))$ to denote the order of a function
$F$. Normally, $F$ is not reported, but instead referred to by this notation.
By convention, the $f$ reported will be close to as small as possible.


Newick Notation
-------------------------------------------------------------------------------

Newick Notation is a format for specifying trees. An example of a Newick tree
is

```
((a:.1,b:.3):1.0, (c:.43,d:.5):1.0);
```

The '`:`' signifies a weight on an edge leading up from the node towards the
root, and is optional.

One advantage of Newick notation is that the grammar[^newick_ref] for this
language is quite simple, and requires only a few productions. This makes it
easy to write a parser to recognize a Newick string. Unfortunately there are
quite a few disadvantages, one of which is the non-uniqueness of the string for
a tree. Specifically, there are many strings which represent the same tree. The
problem is due, largely to the ordering of subtrees. Note that these two trees
are the same

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
is as follows: Given a set of gene rooted trees, do the following for each tree.

-   Set all the edge weights
    -   If the edge is not connected to a leaf, set the weight to 1.
    -   If the edge is connected to a leaf, set it so that the tree will be
    ultrametric[^ultrametric_note].
-   Calculate pairwise distance, put them into a table.

Once all the tables have been calculated, do an entry wise sum, and divide each
entry by the number of gene trees. Informally, think of this as taking the
average of the distance table. Once we have the "average" distance table, run
a tree building algorithm on it, such as Neighbor Joining.

[^1]: Here the term rank is used to mean edge length or weight, and will be
  referred to as edge length in the future.

[^ultrametric_note]: Practically this means that we can take the number of taxa
    and subtract the depth, then use that for the edge weight.

Generalized STAR
-------------------------------------------------------------------------------

Generalized star [@rhodes_star] is version of STAR where the edge weight
schedule is allowed to vary. In this version, edge weights of an equal distance
from the root (using graph-theoretic distance of counting edges in a path) are
given the same weight. Weights need only be non-negative except for one.
Additionally, the edges connected to leaves are adjusted to make the tree
ultrametric, as in STAR. Other than weight assignment, algorithm proceeds the
same.

Additionally, we can attempt to infer the stability of the species tree
produced by STAR by varying the weights and reporting the different trees that
come out of the tree creation step. This is the primary goal of \SunStar.

\texttt{SunStar}
-------------------------------------------------------------------------------

\SunStar currently has only a few options. Here they are:

-   `-f` `--filename`: Filename for a set of Newick trees. These trees must all
together be rooted or unrooted. If they are all unrooted trees, then the `-o`
is required.
-   `-t` `--trials`: The number of trials using randomly genereated schedule of
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
ancestors for each node. Once we have this list, we walk the list, starting from
the rootside, until we find a difference. Then, the remaining nodes in the
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

![An example tree \label{fig:dist_tree}](./figs/distances.pdf)

###Complexity

The complexity of finding the distance between two nodes for this is $\OO(n)$,
where $n$ is the number of taxa on the tree. First, consider the steps involved
with this algorithm. They are

-   Make a list from taxa to root,
-   Walk along that list,
-   Create new path,
-   Sum edge weights over path.

Note that the number of nodes on a rooted tree with $n$ taxa is $2n-1$. Since
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

The algorithm requires a set of $m$ rooted trees. STAR involves three
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

The complexity of STAR is $\OO(mn^3)$. Calculating a single entry in the
distance matrix requires a call to a $\OO(n)$ algorithm, and there are $n^2$
elements in the matrix. Therefore the complexity of computing each table is
$n^2\OO(n) = \OO(n^3)$. We need compute this distance matrix for each of the
$m$ trees in $T$, so $m\OO(n^3) = \OO(mn^3)$.

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
there is no a very large sample of gene trees with error then the tree should
come out the same, despite the 0 weights on the edges. But, we believe that it
is likely that introducing zeros into the system will maximize the exposer of
any errors, and therefore will cause instabilities.

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

At a high level, phylogenetic trees are implemented as doubly linked binary
trees. This is to say, each node on the tree has a parent, and two children.
This also means that there is a sense of directionality towards the root. This
makes finding paths between nodes easier, as we don't have to search in all
directions to find the other node.

There is also an odd concept of an `unroot`. The `unroot` is the special node
that is the root, for the sake of starting to find things on the tree, of an
unrooted tree. It is implemented as a vector of pointers to the node class.

The `node` class contains the information relating to each individual
node in the tree. This includes data including the label of the node, the
weight of the edge towards the parent, and pointers to the children and parent.

The `tree` class is the main workhorse of the tree implementation. It contains
the tree and the `unroot`. Since these trees don't have the normal CRUD
operations, we can predetermine the size of the tree. This allows us to
preallocate a contiguous section of memory to store the `node`s that make up
the tree, as opposed to each node being potentially distant from other nodes.
By doing this, we can take advantage of cache locality, which can be
a potentially large speed up on some systems. This technique I call tree
packing, which is discussed in section \ref{tree-packing}.

The `tree` class is also responsible for setting the weights of the tree. This
can be done few ways. The first is to pass the `set_weight` function a function
that takes a `size_t` or equivalent, and returns a `double`. Another is to
pass it a vector of weights which is indexed by depth. Yet another way is to
pass a constant to the function, which will set the weights to be the constant,
with the exception of the leaves. In all these cases, the ultrametricity of the
tree is preserved, i.e. the leaf edges are adjusted to ensure it.

###Tree Packing

\margalg{
\caption{Tree Packing}
\KwData{List of pointers to tree roots $r$}
\KwResult{A packed representation of a tree in memory}
stack $s$\;
queue $q$\;

\For{$n$ in $r$}{
    push $n$ onto $s$\;
    push $n$ onto $q$\;
}

\While{$s$ is not empty}{
    $n \leftarrow$ top of $s$\;
    pop top of $s$\;
    \If{n has children}{
        push children of $n$ onto $s$\;
        push children of $n$ onto $q$\;
    }
}

$tree \leftarrow$ array the size of $q$\;
\For{n in q}{
    push $n$ onto $tree$\;
    update pointers of $n$\;
}
}


###Set Root

The final major role of the tree is to set and reset the root. Sometimes,
unrooted trees might be passed to this function, but the algorithms only work
on rooted trees. So we must be able to root a tree. To do this, the following
algorithm is used.

\begin{function}
\caption{SwapParent(node $n$, node $p$)}
    \If{p is $n$'s left child}{
        swap $n$'s parent and left child\;
        \SwapParent (left child, n)\;
    }
    \ElseIf{p is $n$'s left child}{
        swap $n$'s parent and right child\;
        \SwapParent (right child, n)\;
    }
\end{function}

\margalg{
\caption{Set Root}
\KwData{List of pointers to roots of subtrees $r$ representing the unroot of an
unrooted tree, outgroup pointer $o$}
\KwResult{New root of a rooted tree, with $o$ as an outgroup}
\If(\tcp*[h]{Unroot is where it should be, just need to root the tree}){$o$ in $r$}{
    make new node $nr$\;
    set $nr$'s children to be the two other nodes in $r$\;
    clear $r$\;
    add $o$ to $r$\;
    add $nr$ to $r$\;
}
\Else{
    make new node $n$\;
    assign $n$'s parents and children to the three nodes of $r$\;
    clear $r$\;
    add $o$ to $r$\;

    \For{child, parent $c$ of $n$}{
        set $c$'s parent to $n$
    }

    $p \leftarrow$ $o$'s parent\;
    \If{$o$ is not the left child of $p$}{
        swap $p$'s left child and right child\;
    }

    set $p$'s left child to \texttt{null}\;
    \SwapParent(p, \texttt{null})
}
}

\texttt{newick.h}{#newick-imp}
-------------------------------------------------------------------------------

In order to work with gene trees, I need to be able to parse Newick Notation.
So, the specification of the parser is as follows.

###Grammar

The grammar used has two 3 productions and 2 terminal lexemes. The terminal
lexemes are labels and weights. The regular expression for labels and weight are

```
    weight = [0-9]+('.'[0-9]*)?('e'[0-9]+)
    label  = [A-z][A-z0-9]
```

In other words, floating point numbers with possible exponential notation, and
c identifiers. Every other literal is punctuation, and is there for structure
only. Below is the grammar that is in the parser

```
    tree    = '(' subtree ',' subtree [',' subtree] ')' [';']
    subtree = label [':' weight]
            | '(' subtree ',' subtree ')' [':' weight]
```

###Parser/Lexer

Since the language is so small, We have made the parser and lexer tightly
integrated. This is to say, there really is no lexer. Instead, characters are
read and dealt with immediately. I read floats with the `stod` function in the
standard `c++` library.

The parser itself is simple top down recursive decent parser. The parser
returns a list of nodes pointers to subtrees that make up the tree. This list
is intended to be the unroot.

\texttt{nj.h}
-------------------------------------------------------------------------------

This header contains my implementation of the Neighbor Joining[@saito87]
algorithm. I need to implement the algorithm, because it ends up being tightly
bound to whatever tree class is used. This makes finding fast libraries hard.
There will not be much discussion of the theory of Neighbor Joining here, just
the details of this implementation.

This implementation starts by putting all the nodes into a list, which should
be thought of as the unroot. It then join pairs until there are only 3 elements
in the unroot left. To join pairs, it calls the find pair routine. The find
pair routine finds[^q_store] the lowest entry of a matrix $Q$. $Q$ is
calculated based of the distance matrix of all the nodes in the unroot. The
lowest element of $Q$ is identified by its row and column (We don't really
care about the value of that element, just that its the lowest or tied for it).
The row and columns of $Q$ are indexes into the unroot, and sow the row column
pair from this step identify a specific two nodes in the unroot, which are to
be joined. So the find pair routine returns this pair.

With this found pair, the implementation joins the pair by removing these
elements from the unroot, making them the children of a new node, and then
putting that new node back into the unroot. We calculate the distance from the
new node to the pair by the three point formula[^three_point], and do the same
for the distance between the new node and all the rest of the nodes in the
unroot. The net effect of this is that the size of unroot shrinks by one every
time the implementation joins a pair.

When the unroot size shrinks to 3, we then then take the remaining three and
use the three point formula to calculate distances between these final three
nodes. With that, complete, we have a complete unrooted tree.

[^q_store]: I don't actually store the matrix Q. Its only ever used to find the
smallest element, so I just calculate it and store the smallest so far. This
is a pretty minor optimization, but saves a $n^2$ memory. Furthermore, $Q$ is
symmetric, so we only compute the top half of the matrix, and skip the diagonal
for another minor optimization.

[^three_point]: In a 4 node tree (trivalent), with leaves $x,y,z$ and central
node $r$, we can calculate and edge length given the distances $d(x,y), d(x,z),
d(y,z)$ by using the three point formula: $d(r,x) = \frac{1}{2}[d(y,x) + d(x,z)
- d(y,z)]$

\texttt{star.h}
-------------------------------------------------------------------------------

The STAR algorithm is straight forward. Most of the complexity of the code
comes from the tree itself, but once that is implemented, the algorithm is
simple. Please refer to algorithm \ref{alg:STAR} for the details.

\margalg{
\label{alg:STAR}
\caption{STAR}
\KwData{A set of trees $T$, with the same set of taxa $X$ where $|X| = n$, with
    weights set.}
\KwResult{A single tree $s$}
$D\leftarrow 0^{n \times n}$\;
\For{$t$ in $T$}{
    $d \leftarrow$ distance table of $t$\;
    $D \leftarrow D + d$\;
}
$D \leftarrow \nicefrac{1}{n} \cdot D$\;
$s \leftarrow$ tree from neighbor joining with $D$\;
\Return{$s$}
}

\texttt{gstar.h}
-------------------------------------------------------------------------------

As stated before, there are two schedules for \SunStar. The first is the
default schedule, which is run when the `--trials` option is not specified. The
other is the randomized schedule, which is run when `--trials` is specified.

The default schedule is an attempt to find an error by brute force. We treat
each level of edges as either having a weight of 0 or 1. We then try every
possible assignment of 0s and 1s.

\margalg{
\caption{GSTAR-default}
\KwData{A set of $n$ trees, $T$, which all have the same set of taxa.}
\KwResult{A set of trees and associated frequencies.}
$R \leftarrow$ map of trees to integers\;
$h \leftarrow \max_{t\in T}(\text{height of }t)$\;
$S \leftarrow$ list of 1's of length $h$\;
\For{$i$ in $1:2^n$}{
    \For{$t$ in $T$}{
        set weights of $t$ by schedule $S$\;
    }
    $t_s \leftarrow$ STAR on $T$\;
    \If{$t_s \not\in R$}{
        $R(t_s) \leftarrow 0$\;
    }
    $R(t_s) \leftarrow R(t_s) + 1$\;
    \tcc{Here we treat S as a binary number, and decrement it}
    decrement $S$\; 
}
\tcc{Normalize the counts into proportions}
\For{$r$ in $R$}{
    $R(r) \leftarrow \nicefrac{1}{2^n} \cdot R(r)$\;
}
\Return{$R$}
}

The random schedule proceeds by randomly generating weights, and assigning them
to the trees in the set. 

\margalg{
\caption{GSTAR-random}
\KwData{A set of $n$ trees, $T$, which all have the same set of taxa. A
number of trials $t$.}
\KwResult{A set of trees and associated frequencies.}
$R \leftarrow$ map of trees to integers\;
$h \leftarrow \max_{t\in T}(\text{height of }t)$\;
$S \leftarrow$ a list of $h$ random numbers, uniformly distributed on $[0,1]$\;
\For{$i$ in $t$}{
    \For{$t$ in $T$}{
        set weights of $t$ by schedule $S$\;
    }
    $t_s \leftarrow$ STAR on $T$\;
    \If{$t_s \not\in R$}{
        $R(t_s) \leftarrow 0$\;
    }
    $R(t_s) \leftarrow R(t_s) + 1$\;
    $S \leftarrow$ a list of $h$ random numbers, uniformly distributed on $[0,1]$\;
}
\tcc{Normalize the counts into proportions}
\For{$r$ in $R$}{
    $R(r) \leftarrow \nicefrac{1}{2^n} \cdot R(r)$\;
}
\Return{$R$}
}

Early tests have shown that the default and random schedules are about
equivalent in detecting error, though the proportions differ.

Conclusion
===============================================================================

We have created a software package, \SunStar, that will use Generalized STAR to
infer the stability of a species tree inferred by STAR. Furthermore, we have
found early results that suggest that this technique might be useful. 

Future Work
===============================================================================

There are many optimizations, enhancements, and investigations we would like to
do with \SunStar. Like many software projects \SunStar could be developed
forever barring resources. What follows is a short list of priority items that
should be done before other options.

Optimizations
-------------------------------------------------------------------------------

There are many optimizations that we would like to implement. They are in no
particular order:

-   Refactor distance table data structure, so that we can take advantage of
    symmetry.
-   Refactor the node class, investigate making the class smaller in memory
-   Implement OpenMP routines on the GSTAR level of computation (Give each STAR
    inference a thread).
-   Optimize the default schedule, find unnecessary steps

Enhancements
-------------------------------------------------------------------------------

There are a few enhancements that are possible for \SunStar. Some of these are
only necessary if the technique turns out to be very successful, but are
included as a roadmap nonetheless.

-   More distributions/schedules for weights in GSTAR.
-   Add MPI support
-   Add "Full Pipeline" inference. Specifically, take in gene sequences and
    return a set of trees with support.

Investigations
-------------------------------------------------------------------------------

To really investigate this method of detecting errors, we would like a study
that

1.   Models gene trees from a species tree using the multispecies coalescent
    model
2.   Generate gene sequences from those gene trees
3.   Infer new gene trees with error from the generated sequences
4.   Run \SunStar and examine the error rate.

Most of these computations are relatively fast, except for step 3. Inferring
gene trees from sequences can take quite a while. Therefore, this investigation
will likely take several months of work.


Code
===============================================================================

The entire code base, including this write up, is located on
[GitHub](https://github.com/computations/SunStar)

References
===============================================================================
