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
    label_set.clear();
    auto ret = parse_label(l, idx);
    REQUIRE(ret == "acc");
    REQUIRE(idx == 4);
}

TEST_CASE("newick label parser, label in newick 2", "[newick]"){
    size_t idx = 5;
    string l = "(acc,b);";
    label_set.clear();
    auto ret = parse_label(l, idx);
    REQUIRE(ret == "b");
    REQUIRE(idx == 6);
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

TEST_CASE("newick parser, simple tree with weights", "[newick][tree]"){
    string t_string = "(a:1.0,b:1.0);";
    string expected = "(a:1.0,b:1.0);";
    tree_t t(t_string);
    t.sort();
    REQUIRE(t.to_string() == expected);
}

TEST_CASE("newick parser, tree with spaces", "[newick][tree]"){
    string t_string = "(a:1.0, b:1.0);";
    string expected = "(a:1.0,b:1.0);";
    tree_t t(t_string);
    t.sort();
    REQUIRE(t.to_string() == expected);
}

TEST_CASE("newick parser, complicated tree", "[newick][tree]"){
    string t_string = "((a:1.0,b:1.0):1.0, (d:1.0,e:1.0):1.0);";
    string expected = "((a:1.0,b:1.0):1.0,(d:1.0,e:1.0):1.0);";
    tree_t t(t_string);
    t.sort();
    REQUIRE(t.to_string() == expected);
}

TEST_CASE("newick parser, complicated tree without weights","[newick][tree]"){
    string t_string = "((a,b),(c,(d,e)));";
    string expected = "((a:1.0,b:1.0):1.0,(c:1.0,(d:1.0,e:1.0):1.0):1.0);";
    tree_t t(t_string);
    t.sort();
    REQUIRE(t.to_string() == expected);
}

TEST_CASE("newick parser, complicated tree with weights", "[newick][tree]"){
    string t_string = "((a:2.3,b:1.4):1.3,(c:3.2,(d:1.2,e:3.1):4.2):9.1);";
    string expected = "((a:2.3,b:1.4):1.3,(c:3.2,(d:1.2,e:3.1):4.2):9.1);";
    tree_t t(t_string);
    t.sort();
    REQUIRE(t.to_string() == expected);
}

TEST_CASE("newick parser edge case 1", "[newick][tree]"){
    string t_string = "((a,((b,c),k)),e);";
    string expected = "((a,((b,c),k)),e);";
    tree_t t(t_string);
    t.clear_weights().sort();
    REQUIRE(expected == t.to_string());
}

TEST_CASE("newick parser edge case 2", "[newick][tree]"){
    string t_string = "((a,( ( b,c),k ) ), e);";
    string expected = "((a,((b,c),k)),e);";
    tree_t t(t_string);
    t.clear_weights().sort();
    REQUIRE(expected == t.to_string());
}
