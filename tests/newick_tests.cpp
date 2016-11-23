#include "catch.hpp"
#include "../src/newick.cpp"

TEST_CASE("newick label parser", "[newick]"){
    string l = "abcd";
    size_t idx = 0;
    auto ret = parse_label(l, idx);
    REQUIRE(ret == l);
    REQUIRE(idx == l.size());
    idx = 0;
    l = "abdc:";
    ret = parse_label(l, idx);
    REQUIRE(ret == l.substr(0,4));
    REQUIRE(idx == l.size()-1);

    idx = 1;
    l = "(a,b);";
    ret = parse_label(l, idx);
    REQUIRE(ret == "a");
    REQUIRE(idx == 2);
}
