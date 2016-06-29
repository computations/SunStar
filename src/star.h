#pragma once
#include "tree.h"
#include <vector>
#include <string>


class star_t{
    public:
        star_t(const std::vector<std::string>&);
    private:
        void  calc_average_distances(std::vector<tree_t>&);

        std::vector<float> _avg_dists;
};
