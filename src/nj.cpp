#include "nj.h"
#include "tree.h"
#include "debug.h"
#include <vector>
using std::vector;
#include <string>
using std::string;
#include <cmath>

nj_t::nj_t(const vector<float>& dists, const vector<string>& labels){
    _dists = dists;
    //because dists is a square matrix, we need to compute the row size.
    //furthemore, since it is a square matrix, the sqrt of the size should an integer
    //there might be some round off, but I tested this on a skylake up to 2^30ish
    //so any number that I would need to calculate this for should work
    _row_size = sqrt(dists.size());
    _tree.resize(_row_size);
    for(size_t i=0;i<_row_size;++i){
        _tree[i]->_label = labels[i];
    }

    while(_row_size>=3){
        join_pair();
    }
}

void nj_t::compute_r(){
    //if the row size is the same as the vec size, then we have already calculated this vector
    //so skip it
    if(_row_size == _r_vec.size()) return;
    _r_vec.resize(_row_size, 0.0);
    for(size_t i=0;i<_row_size;++i){
        for(size_t j=0;j<_row_size;++j){
            _r_vec[i] += _dists[i*_row_size+j];
        }
    }
}

void nj_t::compute_q(){
    //need the r_vec to compute q_vec, so make sure its updated
    compute_r();

    //now to compute the new matrix of values to determine cherry picking
    _q_vec.resize(_row_size*_row_size, 0.0);
    for(size_t i=0;i<_row_size;++i){
        for(size_t j=0;j<_row_size;++j){
            _q_vec[i*_row_size+j] = (_row_size-2)*_dists[i*_row_size+j] - _r_vec[i] - _r_vec[j];
        }
    }
}

void nj_t::find_pair(){
    //compute the matrix Q, which is put into a private data member
    compute_q();

    //find the smallest entry in Q
    //that i,j is the pair we join
    //the way this loop is structured, i>j
    size_t low_i = 0;
    size_t low_j = 0;
    for(size_t i=0;i<_row_size;++i){
        for(size_t j=0;j<i;++j){
            if(_q_vec[i*_row_size+j] <= _q_vec[low_i*_row_size+low_j]){
                low_i = i;
                low_j = j;
            }
        }
    }
    _i = low_i;
    _j = low_j;
}

//in this funciton, I modify the size of _tree and _dists.
//furthermore, _row_size's value changes to keep up.
void nj_t::join_pair(){
    find_pair();
    //make a temp vector
    std::vector<node_t*> tmp_tree;
    for(size_t i=0;i<_row_size;++i){
        if(i==_i || i==_j) continue;
        tmp_tree.push_back(_tree[i]);
    }
    //join nodes and push onto the vector
    node_t* tmp = new node_t;
    tmp->_lchild = _tree[_i];
    tmp->_rchild = _tree[_j];
    tmp_tree.push_back(tmp);
    _tree.swap(tmp_tree);

    //need to calculate new distances for the new node
    //equation boosted from https://en.wikipedia.org/wiki/Neighbor_joining
    //actually, could be simplier
    tmp->_lchild->_weight = .5*_dists[_i*_row_size+_j] + 1/(2*(2*_row_size-2)) * 
        (_r_vec[_i] - _r_vec[_j]);

    tmp->_rchild->_weight = _dists[_i*_row_size+_j] - tmp->_lchild->_weight;

    //integrate the new node into the distance table
    std::vector<float> tmp_dists((_row_size-1)*(_row_size-1),0.0);

    for(size_t i=0;i<_row_size;++i){
        if(i==_i || i==_j) continue;
        size_t cur_i = i;
        if(i>_i) cur_i--;
        if(i>_j) cur_i--;
        for(size_t j=0;j<_row_size;++j){
            if(j==_j || j==_i) continue;
            size_t cur_j = j;
            if(j>_j) cur_j--;
            if(j>_i) cur_j--;
            tmp_dists[cur_i*(_row_size-1)+cur_j] = _dists[i*_row_size-j];
        }
    }
    
    for(size_t i=0;i<_row_size-1;++i){
        tmp_dists[i*(_row_size-1) + (_row_size-2)] = .5 * 
            (_dists[i*_row_size+_i] + _dists[i*_row_size+_j] - _dists[_i*_row_size + _j]);
    }

    _dists.swap(tmp_dists);
}
