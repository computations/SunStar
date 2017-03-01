#pragma once

#include "tree.h"
#include <string>
#include <vector>

std::vector<node_t*> make_tree_from_newick(const std::string& s, size_t& tree_size);
