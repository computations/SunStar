#include "catch.hpp"
#include "../src/star.cpp"

TEST_CASE("star, one tree","[star]"){
    std::string t1="(a:1.0,b:1.0);";
    std::vector<std::string> vt;
    vt.push_back(t1);
    star_t s(vt);
    auto star_tree = s.get_tree();
    star_tree.sort();
    REQUIRE(star_tree.to_string()== t1);
}
