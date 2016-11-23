#include "catch.hpp"
#include "../src/newick.cpp"

TEST_CASE("newick label parser, simple", "[newick]"){
    string l = "abcd";
    size_t idx = 0;
    auto ret = parse_label(l, idx);
    REQUIRE(ret == l);
    REQUIRE(idx == l.size());
}

TEST_CASE("newick label parser, with colon at end", "[newick]"){
    size_t idx = 0;
    string l = "abdc:";
    auto ret = parse_label(l, idx);
    REQUIRE(ret == l.substr(0,4));
    REQUIRE(idx == l.size()-1);
}

TEST_CASE("newick label parser, label in newick", "[newick]"){
    size_t idx = 1;
    string l = "(acc,b);";
    auto ret = parse_label(l, idx);
    REQUIRE(ret == "acc");
    REQUIRE(idx == 4);
}

TEST_CASE("newick weight parser, simple int", "[newick]"){
    string l = "1";
    size_t idx = 0;
    auto ret = parse_weight(l, idx);
    REQUIRE(ret == 1);
    REQUIRE(idx == 1);
}

TEST_CASE("newick weight parser, simple float", "[newick]"){
    string l = "1.0";
    size_t idx = 0;
    auto ret = parse_weight(l, idx);
    REQUIRE(ret == 1.0);
    REQUIRE(idx == 3);
}

TEST_CASE("newick weight parser, float in newick", "[newick]"){
    string l = "(a:1.0,b);";
    size_t idx = 3;
    auto ret = parse_weight(l, idx);
    REQUIRE(ret == 1.0);
    REQUIRE(idx == 6);
}

TEST_CASE("newick parser, simple tree", "[newick][tree]"){
    string t_string = "(a,b);";
    string expected = "(a:1.0,b:1.0);";
    tree_t t(t_string);
    t.sort();
    REQUIRE(t.to_string() == expected);
}
