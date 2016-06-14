#pragma once

#include "tree.h"
#include <string>

node_t* make_tree_from_newick(const std::string& s, size_t& tree_size);
