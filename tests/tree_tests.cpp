#include "catch.hpp"
#include "../src/tree.cpp"

std::vector<std::string> tree_strings {
    "(a,b);",
    "((a,b),(c,d));",
    "((a,b),(c,(d,e)));"
};

std::vector<std::string> result_strings {
    "(a:1.0,b:1.0);",
    "((a:1.0,b:1.0):1.0,(c:1.0,d:1.0):1.0);",
    "((a:1.0,b:1.0):1.0,(c:1.0,(d:1.0,e:1.0):1.0):1.0);"
};

TEST_CASE("tree, default construct", "[tree]"){
    tree_t t;
    REQUIRE(t.to_string() == "");
}

TEST_CASE("tree, string constructor, complicated tree", "[tree]"){
    tree_t t(tree_strings[2]);
    t.sort();
    REQUIRE(t.to_string() == result_strings[2]);
}


TEST_CASE("tree, copy constructor", "[tree]"){
    tree_t t1(tree_strings[0]);
    tree_t t2(t1);
    t1.sort();
    t2.sort();
    REQUIRE(t1.to_string() == result_strings[0]);
    REQUIRE(t1.to_string() == t2.to_string());
}

TEST_CASE("tree, assignment operator", "[tree]"){
    tree_t t1(tree_strings[0]);
    tree_t t2;
    t2=t1;
    t1.sort();
    t2.sort();
    REQUIRE(t1.to_string() == result_strings[0]);
    REQUIRE(t1.to_string() == t2.to_string());
}


TEST_CASE("tree, make label map", "[tree]"){
    tree_t t1(tree_strings[0]);
    auto m = t1.make_label_map();
    REQUIRE(m.size()==2);
    REQUIRE(m.count("a")==1);
    REQUIRE(m.count("b")==1);
}

TEST_CASE("tree, calculate simple distance matrix", "[tree]"){
    tree_t t1(tree_strings[0]);
    auto f = t1.calc_distance_matrix();
    REQUIRE(f[0] == 0.0);
    REQUIRE(f[1] == 2.0);
    REQUIRE(f[2] == 2.0);
    REQUIRE(f[3] == 0.0);
}

TEST_CASE("tree, calculate larger distance matrix", "[tree]"){
    tree_t t1(tree_strings[1]);
    auto f = t1.calc_distance_matrix();
    std::vector<float> r = {
        0,2,4,4,
        2,0,4,4,
        4,4,0,2,
        4,4,2,0
    };
    REQUIRE(r.size() == f.size());
    for(size_t i=0;i<r.size();++i){
        REQUIRE(r[i] == f[i]);
    }
}

TEST_CASE("tree, calculate complicated, not ultrametric distance matrix", "[tree]"){
    tree_t t1(tree_strings[2]);
    auto f = t1.calc_distance_matrix();

    std::vector<float> r = {
        0,2,4,5,5,
        2,0,4,5,5,
        4,4,0,3,3,
        5,5,3,0,2,
        5,5,3,2,0
    };
    REQUIRE(r.size() == f.size());
    for(size_t i=0;i<r.size();++i){
        REQUIRE(r[i] == f[i]);
    }
}

TEST_CASE("tree, setting weights","[tree]"){
    tree_t t1(tree_strings[0]);
    t1.set_weights([](size_t) -> float {return 1.0;});
    t1.sort();
    auto s = t1.to_string();
    REQUIRE(s == "(a:1.0,b:1.0);");
}
