#include "catch.hpp"
#include "../src/gstar.cpp"

TEST_CASE("gstar, functionality test","[gstar]"){
    std::string ns = "(a,b);";
    auto c = gstar({ns});
}

TEST_CASE("gstar, functionality test, larger class", "[gstar]"){
    std::string ns1 = "(a,(b,(c,d)));";
    std::string ns2 = "((a,b),(c,d));";
    auto c = gstar({ns1,ns2});
}

TEST_CASE("gstar, functional test, larger class", "[gstar]"){
    std::string ns = "(((a,b),(c,d)),(((e,f),g),h));";
    auto c = gstar({ns});
}
