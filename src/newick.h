/*
 * newick.h
 * Written by Ben Bettisworth
 */
#pragma once

#include "tree.h"
#include <string>
#include <vector>

//tree_size is a return param
std::vector<node_t*> make_tree_from_newick(const std::string& s, size_t& tree_size);
