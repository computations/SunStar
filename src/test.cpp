#define CATCH_CONFIG_MAIN
#include "catch.hpp"

#include "tree.h"
#include "star.h"
#include "debug.h"
#include "nj.h"

#include <iostream>
using std::cout;
using std::endl;
#include <vector>
using std::vector;
#include <string>
using std::string;


TEST_CASE("making trees from newick notation", "[tree]"){
    vector<string> tree_strings = 
        {"((a:1.0,b:1.0):1.0,(c:1.0,d:1.0):1.0);",
         "((b:1.0,c:1.0):1.0,(a:1.0,d:1.0):1.0);"};

    tree_t t1(tree_strings[0]);
    REQUIRE(t1.to_string() == tree_strings[0]);
}
