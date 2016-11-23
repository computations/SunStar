#include "catch.hpp"
#include "../src/newick.cpp"

TEST_CASE("newick label parser", "[newick]"){
    string l = "abcd";
    size_t idx = 0;
    auto ret = parse_label(l, idx);
    REQUIRE(ret == l);
}
