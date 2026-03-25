#include "CliConfig.h"

#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>

namespace {
int ParseInt(const std::string& value, const std::string& flag) {
    try {
        return std::stoi(value);
    } catch (...) {
        throw std::runtime_error("Invalid integer for " + flag + ": " + value);
    }
}

float ParseFloat(const std::string& value, const std::string& flag) {
    try {
        return std::stof(value);
    } catch (...) {
        throw std::runtime_error("Invalid float for " + flag + ": " + value);
    }
}
}

void PrintUsage() {
    std::cerr
        << "Usage: vdjmatch2 <query.tsv> <target.tsv> [options]\n"
        << "Options:\n"
        << "  --out <path>\n"
        << "  --max-sub <int>\n"
        << "  --max-ins <int>\n"
        << "  --max-del <int>\n"
        << "  --max-edits <int>\n"
        << "  --matrix-path <path>\n"
        << "  --max-cost <float>\n"
        << "  --match-v\n"
        << "  --match-j\n"
        << "  --align\n"
        << "  --gene <value>\n"
        << "  --species <value>\n"
        << "  --epitope <value>\n"
        << "  --threads <int>\n"
        << "  --junction-col <name>\n"
        << "  --v-col <name>\n"
        << "  --j-col <name>\n"
        << "  --epitope-col <name>\n"
        << "  --species-col <name>\n"
        << "  --chain-col <name>\n";
}

CliConfig ParseCli(int argc, char** argv) {
    if (argc < 3) {
        PrintUsage();
        throw std::runtime_error("Two positional TSV inputs are required");
    }

    CliConfig config;
    config.queryPath = argv[1];
    config.targetPath = argv[2];

    for (int i = 3; i < argc; ++i) {
        const std::string arg = argv[i];

        auto requireValue = [&](const std::string& flag) -> std::string {
            if (i + 1 >= argc) {
                throw std::runtime_error("Missing value for " + flag);
            }
            return argv[++i];
        };

        if (arg == "--out") {
            config.outPath = requireValue(arg);
        } else if (arg == "--max-sub") {
            config.maxSub = ParseInt(requireValue(arg), arg);
        } else if (arg == "--max-ins") {
            config.maxIns = ParseInt(requireValue(arg), arg);
        } else if (arg == "--max-del") {
            config.maxDel = ParseInt(requireValue(arg), arg);
        } else if (arg == "--max-edits") {
            config.maxEdits = ParseInt(requireValue(arg), arg);
        } else if (arg == "--matrix-path") {
            config.matrixPath = requireValue(arg);
        } else if (arg == "--max-cost") {
            config.maxCost = ParseFloat(requireValue(arg), arg);
        } else if (arg == "--match-v") {
            config.matchV = true;
        } else if (arg == "--match-j") {
            config.matchJ = true;
        } else if (arg == "--align") {
            config.writeAlignment = true;
        } else if (arg == "--gene") {
            config.gene = requireValue(arg);
        } else if (arg == "--species") {
            config.species = requireValue(arg);
        } else if (arg == "--epitope") {
            config.epitope = requireValue(arg);
        } else if (arg == "--threads") {
            config.threads = ParseInt(requireValue(arg), arg);
        } else if (arg == "--junction-col") {
            config.junctionCol = requireValue(arg);
        } else if (arg == "--v-col") {
            config.vCol = requireValue(arg);
        } else if (arg == "--j-col") {
            config.jCol = requireValue(arg);
        } else if (arg == "--epitope-col") {
            config.epitopeCol = requireValue(arg);
        } else if (arg == "--species-col") {
            config.speciesCol = requireValue(arg);
        } else if (arg == "--chain-col") {
            config.chainCol = requireValue(arg);
        } else if (arg == "--help" || arg == "-h") {
            PrintUsage();
            std::exit(0);
        } else {
            throw std::runtime_error("Unknown argument: " + arg);
        }
    }

    if (config.threads <= 0) {
        throw std::runtime_error("--threads must be positive");
    }

    if (config.maxSub < 0 || config.maxIns < 0 || config.maxDel < 0) {
        throw std::runtime_error("Edit limits must be non-negative");
    }

    return config;
}
