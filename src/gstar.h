#include "star.h"
#include <string>
#include <vector>
#include <utility>

std::vector<std::pair<std::string, double>> gstar(const std::vector<std::string>&,
        std::string = "");
std::vector<std::pair<std::string, double>> gstar(const std::vector<std::string>&,
        const std::string&, size_t, std::string="");
