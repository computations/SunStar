#define CATCH_CONFIG_MAIN
#include "catch.hpp"

#include "../src/debug.h"
bool __PROGRESS_BAR_FLAG = false;
#include "../src/tree.h"
#include "../src/newick.h"
#include "../src/star.h"
#include "../src/nj.h"

#include <iostream>
using std::cout;
using std::endl;
#include <vector>
using std::vector;
#include <string>
using std::string;
