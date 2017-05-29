//nj.cpp
//Ben Bettisworth
//neighbor joiner for gstar
//given a distance table, will create a pretty good guess at a tree
//wikipedia has a pretty good working man's explaination of the alg.
//  https://en.wikipedia.org/wiki/Neighbor_joining

#include "nj.h"
#include "tree.h"
#include "debug.h"

#include <vector>
using std::vector;
#include <string>
using std::string;
#include <queue>
using std::queue;
#include <stack>
using std::stack;
#include <unordered_map>
using std::unordered_map;
#include <limits>
//for std::limits
#include <utility>
//for std::swap
#include <cmath>
#include <iostream>

typedef std::vector<std::vector<double>> d2vector_t;

void delete_rowcol(d2vector_t& m, size_t r){
    d2vector_t tmp_m(m.size()-1, vector<double>(m.size()-1,0.0));
    for(size_t i=0;i<m.size();++i){
        if(i==r) continue;
        size_t _i = i;
        if(i>r) _i--;
        for(size_t j=0;j<m.size();++j){
            if(j==r) continue;
            size_t _j = j;
            if(j>r) _j--;
            tmp_m[_i][_j] = m[i][j];
        }
    }
    m.swap(tmp_m);
}

void add_rowcol(d2vector_t& m){
    for(size_t i = 0;i<m.size();++i){
        m[i].push_back(0.0);
    }
    m.push_back(vector<double>(m.size()+1, 0.0));
}

/*
 * To find the pair to join, we need to calculate a matrix with the values
 *      M[i][j] = (SIZE-2)*dists[i][j] - R[i] - R{j]
 *  Where dists is the distance table. R is an array with the values
 *      R[i] = Sum of dists over row i
 *  and then pick the smallest value of M. Since M will simply be destroyed
 *  after this function, we 
 */
std::pair<size_t, size_t> find_pair(const d2vector_t& dists){
    //calculate R
    size_t row_size = dists.size();
    vector<double> R(row_size, 0);
    for(size_t i = 0;i<row_size;++i){
        for(size_t j=0;j<row_size;++j){
            R[i] += dists[i][j];
        }
    }
    size_t _i=0, _j=0;
    //double lowest=std::numeric_limits<double>::max();
    double lowest = 0.0;
    for(size_t i = 0;i<row_size;++i){
        for(size_t j=i+1;j<row_size;++j){
            double tmp = (row_size-2)*dists[i][j] -R[i] - R[j];
            debug_print("current M[%lu][%lu] = %f, lowest = %f", i, j, tmp, 
                    lowest);
            if(tmp <= lowest){
                _i = i; _j = j;
                lowest = tmp;
                debug_print("new lowest: %f, _i:%lu, _j:%lu", lowest, _i, _j);
            }
        }
    }
    return std::make_pair(_i, _j);
}

/*
 * A slave function of nj. Role is to calculate which pair to join, and to join
 * them. The role of calculating which pair to join is left to find_pair, and
 * then this one will will actually join the pair. Finally, the dists will need
 * to be updated.
 */
void join_pair(d2vector_t& dists, vector<node_t*>& unroot, 
        vector<node_t*>& tree){

    auto p = find_pair(dists);

    debug_print("pair found: (%lu, %lu)", p.first, p.second);
    /*
     * join the last 3
     *  to do that we need to use the three distance formulas
     *    for a graph like
     *        x
     *        |
     *        r
     *       / \
     *      y   z
     *  We can calculate the x-r (d_xr) distance by calculating the following
     *     d_xr = (d_yx + d_xz - d_yz)/2
     *  and we can calculate the other d_ir for i in {x,y,z} the same way
     *  In this case, the three taxa are lchild = x, rchild = y, v = r.
     */

    node_t* lchild = unroot[p.first];
    node_t* rchild = unroot[p.second];
    node_t* v = node_factory(lchild, rchild);
    tree.push_back(v);

    double die = 0.0, dje = 0.0;
    for(size_t k = 0;k<dists.size();++k){
        if(k==p.first || k==p.second) continue;
        die += dists[p.first][k];
    }
    for(size_t k = 0;k<dists.size();++k){
        if(k==p.first || k==p.second) continue;
        dje += dists[p.second][k];
    }
    die/=dists.size();
    dje/=dists.size();
    lchild->_weight = (dists[p.first][p.second] + die - dje)/2.0;
    rchild->_weight = (dists[p.second][p.first] + dje - die)/2.0;

    //remove the two nodes from the unroot, add the new v to the root.
    
    vector<node_t*> tmp_unroot;
    tmp_unroot.reserve(unroot.size()-2);
    for(size_t i=0;i<unroot.size();++i){
        if(i==p.first ||i == p.second) continue;
        tmp_unroot.push_back(unroot[i]);
    }
    tmp_unroot.push_back(v);
    unroot.swap(tmp_unroot);
    //done with the tree update, need to make a new distance table.
    //Doing this with the three point distance table.
    auto tmp_dists = dists;
    delete_rowcol(tmp_dists,p.first);
    auto second_del = p.second;
    if(p.first< p.second) second_del--;
    delete_rowcol(tmp_dists,second_del);
    add_rowcol(tmp_dists);
    /*
     * Here we are goin to use the three point formulas again but a bit
     * differently. Its still the case that for the tree.
     *
     *        x
     *        |
     *        .
     *        .
     *        .
     *        |
     *        r
     *       / \
     *      y   z
     *
     * The following distance relation holds
     *
     *      d_xr = (d_yx + d_xz - d_yz)/2
     *
     * But in this case its dispite the possibly longer path from r to x. Since
     * the only nodes left in the table are "working" nodes, ie they are yet to
     * be joined, we can think of them all begina directly adjacent to r.
     */

    for(size_t i = 0;i<tmp_dists.size()-1;++i){
        size_t cur_i = i;
        if(p.first <= cur_i){cur_i++;}
        if(p.second <= cur_i){cur_i++;}
        tmp_dists[i].back() = (dists[p.first][cur_i] + dists[p.second][cur_i] - dists[p.first][p.second])/2.0;
    }

    for(size_t i = 0;i<tmp_dists.size()-1;++i){
        size_t cur_i = i;
        if(p.first <= cur_i){cur_i++;}
        if(p.second <= cur_i){cur_i++;}
        tmp_dists.back()[i] = (dists[p.first][cur_i] + dists[p.second][cur_i] - dists[p.first][p.second])/2.0;
    }

    dists.swap(tmp_dists);
    debug_d2vector_t("dists, after swap", dists);
}

/*
 * Neighbor Joining aka nj
 * params:
 *  dists:  "Square" vector of distances. actually flat
 *  labels: A list of labels. The idea is that each label's index corrisponds
 *          with the index in the dists. Should be made by inverting the label
 *          map. It is assumed that every taxa (i.e. leaf) has a label.
 * For an introduction to the NJ algorithm, please refer to
 *      https://en.wikipedia.org/wiki/Neighbor_joining
 * and Saito, 1987.
 * Plan:
 *  -   join pairs until 3
 *      -   calculate q
 *      -   pick the smallest value of q (there will be two, pick the top one)
 *      -   make a new node with the chosen pair as children
 *      -   update dists, using the three distsance formula
 *      -   repeat
 *  -   at 3 nodes left, do a final join that just sets the last weights.
 *  -   return created tree
 * Varaibles:
 *  unroot: The root of the unrooted tree. Because NJ makes and deals with
 *          unrooted trees, we can't actually designate a root. But we can make
 *          a special interior node that is the "start". I call that node the
 *          unroot, because its cute.
 */
tree_t nj(const vector<double>& d, const vector<string>& labels){
    size_t row_size = labels.size();
    d2vector_t dists(row_size, vector<double>(row_size));
    //convert the dists to a 2d vector
    for(size_t i = 0; i<row_size;++i){
        for(size_t j = 0; j<row_size;++j){
            dists[i][j] = d[i*row_size+j];
        }
    }
    debug_d2vector_t("dists, after conversion", dists);

    vector<node_t*> unroot;
    vector<node_t*> tree;
    for(auto && l: labels){
        debug_print("making a new node with label: %s", l.c_str());
        auto tn = new node_t(l);
        unroot.push_back(tn);
        tree.push_back(tn);
    }
    while(unroot.size() > 3){
        join_pair(dists, unroot, tree);
    }
    /*
     * Time for the FINAL JOIN (DU DU DU DUUU)!
     * to do this, we have to take the final three taxa in the unroot and set
     * their distances based on the three point formula. The tree looks like
     *              0
     *              |
     *              u
     *             / \
     *            1   2
     */ 
    if(unroot.size() == 3){
        unroot[0]->_weight = (dists[0][1] + dists[0][2] - dists[1][2])/2.0;
        unroot[1]->_weight = (dists[1][0] + dists[1][2] - dists[0][2])/2.0;
        unroot[2]->_weight = (dists[2][0] + dists[2][1] - dists[0][1])/2.0;
    }
    if(unroot.size() == 2){
        unroot[0]->_weight = dists[0][1]/2.0;
        unroot[1]->_weight = dists[0][1]/2.0;
    }
    tree_t ret(unroot);
    for(auto&& n: tree){
        delete n;
    }
    return ret;
}
