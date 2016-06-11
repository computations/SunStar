#include "tree.h"

tree_t::tree_t(const tree_t& t){
    _size = t._size;
    _tree = new node_t[_size];
}

tree_t::~tree_t(){
    if(_tree)
        delete[] _tree;
}

tree_t& tree_t::operator=(const tree_t& t){
    if(_tree)
        delete[] _tree;
    
    _size = t._size;
    _tree = new node_t[_size];
    for(size_t i=0;i<_size;++i){
        _tree[i] = t._tree[i];
    }
    return *this;
}
