/*
 * newick.h
 * Written by Ben Bettisworth
 */
#pragma once

#include "tree.h"
#include <string>
#include <vector>

/*
 * Function to parse a newick string and create a tree. There are two
 * parameters
 *
 *  newick_string    The newick string which needs to be parsed. 
 *
 *  tree_size        Return parameter. The size of the tree.
 *
 * This function is really only intended to be used by the tree class. Since
 * this is the case, it returns a vector of nodes with the internal topology of
 * the newick string.
 */
std::vector<node_t*> make_tree_from_newick(const std::string& newick_string, 
        size_t& tree_size);
