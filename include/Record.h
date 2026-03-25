#pragma once

#include <cstddef>
#include <string>
#include <vector>

struct Record {
    std::size_t rowIndex = 0;
    std::string junctionAA;
    std::string vGene;
    std::string jGene;
    std::string epitope;
    std::string species;
    std::string chain;
    std::vector<std::string> rawFields;
};
