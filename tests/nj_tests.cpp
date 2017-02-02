#include "catch.hpp"
#include "../src/nj.cpp"

TEST_CASE("nj with simple distance table", "[nj]"){
    std::vector<double> d = {0.0,1.0,
                            1.0,0.0};
    std::vector<std::string> l {"a", "b"};
    nj_t n(d,l);
    auto nj_tree = n.get_tree();
    nj_tree.sort();
    REQUIRE(nj_tree.to_string() == "(a:0.5,b:0.5);");
}

TEST_CASE("nj with larger distance table", "[nj]"){
    std::vector<double> d = {0.0, 1.0, 1.0,
                            1.0, 0.0, 1.0,
                            1.0, 1.0, 0.0};
    std::vector<std::string> l {"a", "b", "c"};
    nj_t n(d,l);
    REQUIRE(n.get_tree().to_string() == "(a:0.5,b:0.5,c:0.5);");
}

TEST_CASE("nj with largest distance table", "[nj][regression]"){
    std::vector<double> d = {0.0, 1.0, 2.5, 2.5,
                            1.0, 0.0, 2.5, 2.5,
                            2.5, 2.5, 0.0, 1.0,
                            2.5, 2.5, 1.0, 0.0};
    std::vector<std::string> l {"a", "b", "c", "d"};
    nj_t n(d,l);
    REQUIRE(n.get_tree().to_string() == "(a:0.5,b:0.5,(d:0.5,c:0.5):1.5);");
}

TEST_CASE("nj and then setting weights", "[nj]"){
    std::vector<double> d = {0.0,1.0,
                            1.0,0.0};
    std::vector<std::string> l {"a", "b"};
    nj_t n(d,l);
    auto nj_tree = n.get_tree();
    nj_tree.set_weights([](size_t) -> double {return 1.0;});
    nj_tree.sort();
    REQUIRE(nj_tree.to_string() == "(a:1.0,b:1.0);");
}

TEST_CASE("nj, larger distance table, complicated tree", "[nj][regression]"){
    std::vector<double> d = {0.0, 0.5, 1.0, 2.75, 2.5, 2.75, 
                             0.5, 0.0, 1.0, 2.75, 2.5, 2.75,
                             1.0, 1.0, 0.0, 2.75, 2.5, 2.75,
                             2.75, 2.75, 2.75, 0.0, 0.75, 2.0,
                             2.25, 2.25, 2.5, 0.75, 0.0, 1.75,
                             2.75, 2.75, 2.75, 2.0, 1.75, 0.0};
    
    std::vector<std::string> l {"a", "b", "c", "d", "e", "f"};
    nj_t n(d,l);
    REQUIRE(n.get_tree().to_string() == "(a:0.2,b:0.2,(((e:0.4,d:0.4),f):1.2,c:1.2):0.2);");
}

TEST_CASE("nj, from tree distance matrix", "[nj]"){
    std::string newick_string = "(((a,b),(c,d)),(((e,f),g),h));";
    tree_t t(newick_string);
    auto dists = t.calc_distance_matrix();
    auto lm = t.make_label_map();
    nj_t n(dists, lm);
    debug_print(n.get_tree().to_string().c_str());

}
