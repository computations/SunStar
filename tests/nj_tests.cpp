#include "catch.hpp"
#include "../src/nj.cpp"

TEST_CASE("nj, column and row delete", "[nj]"){
    const d2vector_t d {{4.4, 3.3, 32.0}, 
                        {2.1, 45.0, 22.0}, 
                        {55, 49, 31}};
    auto d1 = d;
    auto d2 = d;
    auto d3 = d;
    delete_rowcol(d1, 0);
    REQUIRE(d1[0][0]==45);
    REQUIRE(d1[0][1]==22.0);
    REQUIRE(d1[1][0]==49);
    REQUIRE(d1[1][1]==31);

    delete_rowcol(d2, 1);
    REQUIRE(d2[0][0]==4.4);
    REQUIRE(d2[0][1]==32.0);
    REQUIRE(d2[1][0]==55);
    REQUIRE(d2[1][1]==31);

    delete_rowcol(d3, 2);
    REQUIRE(d3[0][0]==4.4);
    REQUIRE(d3[0][1]==3.3);
    REQUIRE(d3[1][0]==2.1);
    REQUIRE(d3[1][1]==45.0);
}

TEST_CASE("nj, add row and col", "[nj]"){
    d2vector_t d {{4.4}};
    add_rowcol(d);
    REQUIRE(d[1][0]==0.0);
    REQUIRE(d[1][1]==0.0);
    REQUIRE(d[0][1]==0.0);
}

TEST_CASE("new nj with simple distance table", "[nj]"){
    std::vector<double> d = {0.0,1.0,
                            1.0,0.0};
    std::vector<std::string> l {"a", "b"};
    auto t = nj(d, l);
    t.sort();
    REQUIRE(t.to_string() == "(a:0.5,b:0.5);");
}

TEST_CASE("new nj with large enough distance table", "[nj][regression]"){
    std::vector<double> d = {0.0, 1.0, 2.5, 2.5,
                            1.0, 0.0, 2.5, 2.5,
                            2.5, 2.5, 0.0, 1.0,
                            2.5, 2.5, 1.0, 0.0};
    std::vector<std::string> l {"a", "b", "c", "d"};
    auto n = nj(d,l);
    n.sort();
    REQUIRE(n.to_string() =="((a:0.5,b:0.5):3.5,c:0.5,d:0.5);");
}

TEST_CASE("new nj, tree from wikipedia", "[nj][wiki]"){
    std::vector<double> d = { 0.0, 5.0, 9.0, 9.0, 8.0,
                              5.0, 0.0, 10.0, 10.0, 9.0,
                              9.0, 10.0, 0.0, 8.0, 7.0,
                              9.0, 10.0, 8.0, 0.0, 3.0,
                              8.0, 9.0, 7.0, 3.0, 0.0 };
    std::vector<std::string> l {"a", "b", "c", "d", "e"};
    auto tree = nj(d,l);
    tree.sort();
    tree.clear_weights();
    REQUIRE(tree.to_string() == "(((a,b),c),d,e);");
}
