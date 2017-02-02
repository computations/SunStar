#include "catch.hpp"
#include "../src/tree.cpp"

std::vector<std::string> tree_strings {
    "(a,b);",
    "((a,b),(c,d));",
    "((a,b),(c,(d,e)));",
    "(((a,b),(c,(d,e))),f);",
    "(((a,b),(c,d)),(((e,f),g),h));"
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


TEST_CASE("tree, make label map, small", "[tree]"){
    tree_t t1(tree_strings[0]);
    auto m = t1.make_label_map();
    REQUIRE(m.size()==2);
    REQUIRE(m.count("a")==1);
    REQUIRE(m.count("b")==1);
}

TEST_CASE("tree, make label map, medium", "[tree]"){
    tree_t t1(tree_strings[1]);
    auto m = t1.make_label_map();
    REQUIRE(m.size()==4);
    REQUIRE(m.count("a")==1);
    REQUIRE(m.count("b")==1);
    REQUIRE(m.count("c")==1);
    REQUIRE(m.count("d")==1);
}

TEST_CASE("tree, make label map, large", "[tree]"){
    tree_t t1(tree_strings[2]);
    auto m = t1.make_label_map();
    REQUIRE(m.size()==5);
    REQUIRE(m.count("a")==1);
    REQUIRE(m.count("b")==1);
    REQUIRE(m.count("c")==1);
    REQUIRE(m.count("d")==1);
    REQUIRE(m.count("e")==1);
}

TEST_CASE("tree, make label map, larger", "[tree]"){
    tree_t t1(tree_strings[3]);
    auto m = t1.make_label_map();
    REQUIRE(m.size()==6);
    REQUIRE(m.count("a")==1);
    REQUIRE(m.count("b")==1);
    REQUIRE(m.count("c")==1);
    REQUIRE(m.count("d")==1);
    REQUIRE(m.count("e")==1);
    REQUIRE(m.count("f")==1);
}

TEST_CASE("tree, make label map, largest", "[tree]"){
    tree_t t1(tree_strings[4]);
    auto m = t1.make_label_map();
    REQUIRE(m.size()==8);
    REQUIRE(m.count("a")==1);
    REQUIRE(m.count("b")==1);
    REQUIRE(m.count("c")==1);
    REQUIRE(m.count("d")==1);
    REQUIRE(m.count("e")==1);
    REQUIRE(m.count("f")==1);
    REQUIRE(m.count("g")==1);
    REQUIRE(m.count("h")==1);
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

TEST_CASE("tree, calculate large tree distance matrix", "[tree]"){
    tree_t t1(tree_strings[3]);
    t1.set_weights(1.0);
    t1.sort();
    auto f = t1.calc_distance_matrix();

    std::vector<double> r = {
        0, 4, 6, 6, 6, 8,
        4, 0, 6, 6, 6, 8,
        6, 6, 0, 4, 4, 8,
        6, 6, 4, 0, 2, 8,
        6, 6, 4, 2, 0, 8,
        8, 8, 8, 8, 8, 0
    };
    REQUIRE(r.size() == f.size());
    for(size_t i=0;i<r.size();++i){
        REQUIRE(r[i] == f[i]);
    }
}

TEST_CASE("tree, calculate larger tree distance matrix", "[tree]"){
    tree_t t1(tree_strings[4]);
    t1.set_weights(1.0);
    t1.sort();
    auto f = t1.calc_distance_matrix();

    std::vector<double> r = {
        0.0, 2.0, 4.0, 4.0, 7.0, 7.0, 6.0, 5.0,
        2.0, 0.0, 4.0, 4.0, 7.0, 7.0, 6.0, 5.0,
        4.0, 4.0, 0.0, 2.0, 7.0, 7.0, 6.0, 5.0,
        4.0, 4.0, 2.0, 0.0, 7.0, 7.0, 6.0, 5.0,
        7.0, 7.0, 7.0, 7.0, 0.0, 2.0, 3.0, 4.0,
        7.0, 7.0, 7.0, 7.0, 2.0, 0.0, 3.0, 4.0,
        6.0, 6.0, 6.0, 6.0, 3.0, 3.0, 0.0, 3.0,
        5.0, 5.0, 5.0, 5.0, 4.0, 4.0, 3.0, 0.0
    };
    REQUIRE(r.size() == f.size());
    for(size_t i=0;i<r.size();++i){
        REQUIRE(r[i] == f[i]);
    }
}

TEST_CASE("tree, setting weights with double","[tree]"){
    tree_t t1(tree_strings[0]);
    t1.set_weights(3.2, 3.2);
    t1.sort();
    auto s = t1.to_string();
    REQUIRE(s == "(a:3.2,b:3.2);");
}

TEST_CASE("tree, setting weights with a vector","[tree]"){
    tree_t t1(tree_strings[0]);
    std::vector<double> ws{4.4, 2.0, 3.0};
    t1.set_weights(ws, 4.4);
    t1.sort();
    auto s = t1.to_string();
    REQUIRE(s == "(a:4.4,b:4.4);");
}

TEST_CASE("tree, setting weights with a function","[tree]"){
    tree_t t1(tree_strings[0]);
    t1.set_weights([](size_t d) -> double {return d*2.2;}, 2.2);
    t1.sort();
    auto s = t1.to_string();
    REQUIRE(s == "(a:2.2,b:2.2);");
}
