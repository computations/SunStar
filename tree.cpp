#include "tree.h"
#include "newick.h"
#include <string>
using std::string;

tree_t::tree_t(const tree_t& t){
    _size = t._size;
    _tree = new node_t[_size];
}

//takes a tree specified by newick notation, and makes a tree_t
tree_t::tree_t(const string& newick){
    _tree = make_tree_from_newick(newick, _size);
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
