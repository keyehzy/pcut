#include <pcut/coefficients.hpp>
#include <catch2/catch_test_macros.hpp>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <string>

TEST_CASE("Universal coefficients match every supplied Appendix C entry", "[coefficients][literature]") {
    const pcut::Coefficients c({-2,-1,0,1,2},6);
    std::ifstream file(PCUT_TEST_DATA "/coefficients.tsv");
    REQUIRE(file.good());
    std::string line;
    std::size_t checked = 0;
    while (std::getline(file,line)) {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream input(line);
        std::string word_text, value;
        input >> word_text >> value;
        std::replace(word_text.begin(),word_text.end(),',',' ');
        std::istringstream words(word_text);
        pcut::Word word;
        int m;
        while (words >> m) word.push_back(m);
        INFO(line);
        REQUIRE(c.at(word) == pcut::Rational(value));
        ++checked;
    }
    REQUIRE(checked > 400);
    for (const auto& [word,value] : c.terms()) {
        auto adjoint = word;
        std::reverse(adjoint.begin(),adjoint.end());
        for (auto& m : adjoint) m = -m;
        REQUIRE(c.at(adjoint) == value);
    }
}
TEST_CASE("General charge alphabet and coefficient input validation", "[coefficients]") {
    const pcut::Coefficients c({-3,0,3},4);
    REQUIRE(c.at({-3,3}) == pcut::Rational(-1)/3);
    REQUIRE(c.at({3,-3}) == pcut::Rational(1)/3);
    REQUIRE(c.at({0}) == 1);
    REQUIRE(c.at({0,0}) == 0);
    REQUIRE_THROWS_AS(c.at({1}),std::invalid_argument);
    REQUIRE_THROWS_AS(c.at({}),std::out_of_range);
    REQUIRE_THROWS_AS(pcut::Coefficients({0},65),std::invalid_argument);
    REQUIRE(pcut::Coefficients({},0).terms().empty());
}
