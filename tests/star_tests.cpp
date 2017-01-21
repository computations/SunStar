#include "catch.hpp"
#include "../src/star.cpp"

TEST_CASE("star, one tree","[star]"){
    std::string t1="(a,b);";
    std::vector<std::string> vt;
    vt.push_back(t1);
    star_t s(vt);
    REQUIRE(s.get_tree().to_string() == t1);
}
