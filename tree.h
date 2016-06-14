//tree.h
//Ben Bettisworth
//2016-06-10
//Intedned to be a packed tree for the purposes of phylogenetics
//specificially, the tree has these atributes:
//  saturated binary tree
//  fixed number of nodes, usually given not by the nodes, but the number of leaves
//  branch lengths,
//  maybe no root?

#pragma once

#include <cstdlib>
#include <string>

class node_t{
    public:
        node_t(): _parent(0), _weight(0.0) {};
        node_t(size_t s, float w): _parent(s), _weight(w){};
        std::string _label;
        size_t _parent;
        float _weight;
        size_t _lchild;
        size_t _rchild;
        bool _children;
};

class tree_t{
    public:
        tree_t(): _tree(0), _size(0){};
        tree_t(size_t s): _size(s){
            _tree = new node_t[_size];
        };
        tree_t(const tree_t&);
        tree_t(const std::string&);
        ~tree_t();
        tree_t& operator=(const tree_t&);
    private:
        node_t* _tree;
        size_t _size;
};
