#include "catch.hpp"
#include "../src/tree.cpp"

TEST_CASE("tree, default construct", "[tree]"){
    tree_t t;
    REQUIRE(t.to_string() == "");
}
