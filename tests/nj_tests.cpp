#include "catch.hpp"
#include "../src/nj.cpp"

TEST_CASE("nj with simple distance table", "[nj]"){
    std::vector<float> d = {0.0,1.0,
                            1.0,0.0};
    std::vector<std::string> l {"a", "b"};
    nj_t n(d,l);
    auto nj_tree = n.get_tree();
    nj_tree.sort();
    REQUIRE(nj_tree.to_string() == "(a:0.5,b:0.5);");
}

TEST_CASE("nj with larger distance table", "[nj]"){
    std::vector<float> d = {0.0, 1.0, 2.0,
                            1.0, 0.0, 1.0,
                            2.0, 1.0, 0.0};
    std::vector<std::string> l {"a", "b", "c"};
    nj_t n(d,l);
    REQUIRE(n.get_tree().to_string() == "((a,b),c);");
}

TEST_CASE("nj with largest distance table", "[nj]"){
    std::vector<float> d = {0.0, 1.0, 2.0, 3.0,
                            1.0, 0.0, 1.0, 2.0,
                            2.0, 1.0, 0.0, 1.0, 
                            3.0, 2.0, 1.0, 0.0};
    std::vector<std::string> l {"a", "b", "c", "d"};
    nj_t n(d,l);
    REQUIRE(n.get_tree().to_string() == "((a,b),c);");
}
