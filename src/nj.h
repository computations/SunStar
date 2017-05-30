//nj.h
//Ben Bettisworth
#pragma once

#include "tree.h"
#include <vector>
#include <string>

/*
 * Makes a new tree given a distance table. The vector of strings are labels
 * for the rows of the distance table. As such, they need to be matched up. The
 * first parameter is a distance table or matrix. As such, while it is only one
 * dimensional, it will be treated as a 2 dimensional vector. Therefore, the
 * size of the vector needs to be a perfect square. If it really is a distance
 * table, then there should be no problem.
 */
tree_t nj(const std::vector<double>&, const std::vector<std::string>&);
