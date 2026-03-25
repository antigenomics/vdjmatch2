#pragma once

#include <optional>
#include <string>

struct CliConfig {
    std::string queryPath;
    std::string targetPath;
    std::string outPath = "match_result.tsv";

    int maxSub = 1;
    int maxIns = 0;
    int maxDel = 0;
    int maxEdits = 1;

    std::optional<std::string> matrixPath = std::nullopt;
    float maxCost = 6.0f;

    bool matchV = false;
    bool matchJ = false;
    bool writeAlignment = false;

    std::string gene = "";
    std::string species = "HomoSapiens";
    std::optional<std::string> epitope = std::nullopt;

    int threads = 4;

    std::optional<std::string> junctionCol = std::nullopt;
    std::optional<std::string> vCol = std::nullopt;
    std::optional<std::string> jCol = std::nullopt;
    std::string epitopeCol = "antigen.epitope";
    std::string speciesCol = "species";
    std::string chainCol = "gene";
};

CliConfig ParseCli(int argc, char** argv);
void PrintUsage();
