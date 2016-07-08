//tree.h
//Ben Bettisworth
//2016-06-10
//Intedned to be a packed tree for the purposes of phylogenetics
//specificially, the tree has these atributes:
//  saturated binary tree
//  fixed number of nodes, usually given not by the nodes, but the number of leaves
//  branch lengths,
//  maybe no root?
//      this changed to using an ''unroot''
//      the unroot is just a collection of three rooted trees
//      which get joined at the unroot, to make the three rooted trees an unrooted tree

#pragma once

#include <cstdlib>
#include <string>
#include <unordered_map>
#include <vector>
#include <functional>

class node_t{
    public:
        node_t(): _parent(0), _weight(0.0), _lchild(0), _rchild(0), _children(false) {};
        node_t(node_t* s, float w): _parent(s), _weight(w){};


        std::string to_string();

        void set_weights(std::function<float(size_t)>, size_t);
        void update_children(std::unordered_map<node_t*, node_t*>);
        std::string sort();

        size_t count_nodes();

        std::string _label;
        node_t* _parent;
        float _weight;
        node_t* _lchild;
        node_t* _rchild;
        bool _children;
};

class tree_t{
    public:
        tree_t(): _tree(0), _size(0){};
        tree_t(size_t s): _size(s){
            _tree = new node_t[_size];
        };
        tree_t(const std::vector<node_t*>&);
        //tree_t(node_t*,size_t, const std::vector<node_t*>&);
        tree_t(const tree_t&);
        tree_t(const std::string&);
        ~tree_t();
        tree_t& operator=(tree_t);

        std::vector<node_t*> get_parents_of(node_t*);

        std::string to_string() const;
        std::string print_labels() const;

        std::unordered_map<std::string, size_t> make_label_map();

        float* calc_distance_matrix(const std::unordered_map<std::string, size_t>&);
        void calc_distance_matrix(const std::unordered_map<std::string, size_t>&, float*);
        float calc_distance(node_t*, node_t*);
        float parent_distance(node_t* child, node_t* parent);
        void set_weights(const std::vector<float>&);
        void set_weights(std::function<float(size_t)>);
        void sort();

    private:
        
        void make_flat_tree(const std::vector<node_t*>&);

        //for unrooted trees, we can thing of them as being 3 unrooted trees
        //we join them together in a vector
        //this intersection is the new root of an unrooted tree
        //hence, the unroot
        std::vector<node_t*> _unroot;
        node_t* _tree;
        size_t _size;
};
