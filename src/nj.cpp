#include "nj.h"
#include "tree.h"
#include "debug.h"
#include <vector>
using std::vector;
#include <cmath>

nj_t::nj_t(std::vector<float> dists){
    _dists = dists;
    //because dists is a square matrix, we need to compute the row size.
    //furthemore, since it is a square matrix, the sqrt of the size should an integer
    //there might be some round off, but I tested this on a skylake up to 2^30ish
    //so any number that I would need to calculate this for should work
    _row_size = sqrt(dists.size());
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

void find_pair(){
}
