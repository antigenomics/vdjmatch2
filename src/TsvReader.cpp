#include "TsvReader.h"

#include <algorithm>
#include <cctype>
#include <fstream>
#include <stdexcept>
#include <string_view>
#include <unordered_map>

namespace {
std::string Trim(std::string s) {
    while (!s.empty() && std::isspace(static_cast<unsigned char>(s.back()))) s.pop_back();
    std::size_t start = 0;
    while (start < s.size() && std::isspace(static_cast<unsigned char>(s[start]))) ++start;
    return s.substr(start);
}

void TrimCR(std::string& s) {
    if (!s.empty() && s.back() == '\r') s.pop_back();
}

std::vector<std::string> SplitTSV(const std::string& line) {
    std::vector<std::string> fields;
    std::size_t start = 0;
    while (start <= line.size()) {
        std::size_t end = line.find('\t', start);
        if (end == std::string::npos) end = line.size();
        fields.emplace_back(line.substr(start, end - start));
        start = end + 1;
        if (end == line.size()) break;
    }
    return fields;
}

std::string NormalizeSequence(std::string sequence) {
    for (char& c : sequence) {
        c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
    }
    return sequence;
}

bool IsValidSequence(const std::string& sequence) {
    if (sequence.empty()) return false;
    for (char c : sequence) {
        if (c < 'A' || c > 'Z') return false;
    }
    return true;
}

int ResolveColumn(const std::unordered_map<std::string, int>& columns,
                  const std::optional<std::string>& explicitName,
                  const std::vector<std::string>& fallbackNames,
                  bool required) {
    if (explicitName) {
        auto it = columns.find(*explicitName);
        if (it == columns.end()) {
            if (required) {
                throw std::runtime_error("Missing required column: " + *explicitName);
            }
            return -1;
        }
        return it->second;
    }

    for (const auto& name : fallbackNames) {
        auto it = columns.find(name);
        if (it != columns.end()) return it->second;
    }

    if (required) {
        throw std::runtime_error("Missing required junction column");
    }
    return -1;
}

std::string FieldAt(const std::vector<std::string>& fields, int idx) {
    if (idx < 0 || idx >= static_cast<int>(fields.size())) return {};
    return fields[idx];
}
}

TsvReadResult ReadTsv(const std::string& path, const CliConfig& config) {
    std::ifstream input(path);
    if (!input) {
        throw std::runtime_error("Failed to open TSV: " + path);
    }

    TsvReadResult result;
    std::string headerLine;
    if (!std::getline(input, headerLine)) {
        throw std::runtime_error("Empty TSV: " + path);
    }
    TrimCR(headerLine);
    result.header = SplitTSV(headerLine);

    std::unordered_map<std::string, int> columns;
    for (int i = 0; i < static_cast<int>(result.header.size()); ++i) {
        columns.emplace(result.header[i], i);
    }

    const int junctionCol = ResolveColumn(columns, config.junctionCol, {"junction_aa", "cdr3"}, true);
    const int vCol = ResolveColumn(columns, config.vCol, {"v_call", "v.segm"}, false);
    const int jCol = ResolveColumn(columns, config.jCol, {"j_call", "j.segm"}, false);
    const int epitopeCol = ResolveColumn(columns, std::optional<std::string>(config.epitopeCol), {}, false);
    const int speciesCol = ResolveColumn(columns, std::optional<std::string>(config.speciesCol), {}, false);
    const int chainCol = ResolveColumn(columns, std::optional<std::string>(config.chainCol), {}, false);

    std::string line;
    std::size_t rowIndex = 0;
    while (std::getline(input, line)) {
        TrimCR(line);
        if (line.empty()) continue;

        auto fields = SplitTSV(line);
        auto junction = NormalizeSequence(Trim(FieldAt(fields, junctionCol)));
        if (junction.empty()) {
            ++result.skippedMissingJunction;
            continue;
        }
        if (!IsValidSequence(junction)) {
            ++result.skippedInvalidSequence;
            continue;
        }

        Record record;
        record.rowIndex = rowIndex++;
        record.junctionAA = std::move(junction);
        record.vGene = Trim(FieldAt(fields, vCol));
        record.jGene = Trim(FieldAt(fields, jCol));
        record.epitope = Trim(FieldAt(fields, epitopeCol));
        record.species = Trim(FieldAt(fields, speciesCol));
        record.chain = Trim(FieldAt(fields, chainCol));
        record.rawFields = std::move(fields);

        if (!config.gene.empty() && record.chain != config.gene) {
            ++result.skippedByFilter;
            continue;
        }
        if (!config.species.empty() && record.species != config.species) {
            ++result.skippedByFilter;
            continue;
        }
        if (config.epitope && epitopeCol >= 0 && record.epitope != *config.epitope) {
            ++result.skippedByFilter;
            continue;
        }

        result.records.push_back(std::move(record));
    }

    return result;
}
