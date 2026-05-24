#include "Trie.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <tuple>

namespace {
    constexpr float kInf = std::numeric_limits<float>::infinity();
    constexpr float kEps = 1e-6f;

    int LetterIndex(char c) {
        if (c < 'A' || c > 'Z') throw std::runtime_error("Invalid amino acid letter");
        return c - 'A';
    }

    int OpPriority(AlignmentOpType op) {
        if (op == AlignmentOpType::Match || op == AlignmentOpType::Substitution) return 0;
        if (op == AlignmentOpType::Deletion) return 1;
        return 2;
    }

    struct MatrixCell {
        float cost = kInf;
        int prevI = -1;
        int prevJ = -1;
        AlignmentOpType op = AlignmentOpType::Match;
        bool reachable = false;
    };

    bool Better(float newCost, AlignmentOpType newOp, const MatrixCell& current) {
        if (!current.reachable) return true;
        if (newCost + kEps < current.cost) return true;
        return std::fabs(newCost - current.cost) <= kEps && OpPriority(newOp) < OpPriority(current.op);
    }

    struct EditCounts {
        int sub = 0;
        int ins = 0;
        int del = 0;
        int total() const { return sub + ins + del; }
    };

    bool Dominates(const EditCounts& a, const EditCounts& b) {
        return a.sub <= b.sub && a.ins <= b.ins && a.del <= b.del &&
               (a.sub < b.sub || a.ins < b.ins || a.del < b.del);
    }

    void AddState(std::vector<EditCounts>& states, EditCounts candidate) {
        for (const auto& state : states) {
            if ((state.sub == candidate.sub && state.ins == candidate.ins && state.del == candidate.del) || Dominates(state, candidate)) {
                return;
            }
        }
        states.erase(std::remove_if(states.begin(), states.end(), [&](const EditCounts& state) {
            return Dominates(candidate, state);
        }), states.end());
        states.push_back(candidate);
    }
}

Trie::Trie(const std::vector<std::string>& sequences,
           const std::vector<std::string>& vGenes,
           const std::vector<std::string>& jGenes)
    : root_(std::make_unique<Node>()), sequences_(sequences), vGenes_(vGenes), jGenes_(jGenes) {
    BuildTrie();
}

void Trie::LoadSubstitutionMatrix(const std::string& matrixPath) {
    matrix_.Load(matrixPath);
}

void Trie::BuildTrie() {
    for (int index = 0; index < static_cast<int>(sequences_.size()); ++index) {
        Node* node = root_.get();
        for (char c : sequences_[index]) {
            int childIndex = LetterIndex(c);
            if (!node->children[childIndex]) {
                node->children[childIndex] = std::make_unique<Node>();
                node->children[childIndex]->letter = c;
            }
            node = node->children[childIndex].get();
        }
        node->indices.push_back(index);
    }
}

std::vector<std::pair<std::size_t, int>> Trie::SearchIndices(
    const std::string& query,
    int maxSubstitution,
    int maxInsertion,
    int maxDeletion,
    int maxEdits,
    const std::optional<std::string>& vGeneFilter,
    const std::optional<std::string>& jGeneFilter) const {

    std::vector<std::pair<std::size_t, int>> results;
    for (std::size_t index = 0; index < sequences_.size(); ++index) {
        if (vGeneFilter && vGenes_[index] != *vGeneFilter) continue;
        if (jGeneFilter && jGenes_[index] != *jGeneFilter) continue;
        auto distance = BoundedEditDistance(query, sequences_[index], maxSubstitution, maxInsertion, maxDeletion, maxEdits);
        if (distance) results.emplace_back(index, *distance);
    }
    return results;
}

std::optional<int> Trie::BoundedEditDistance(const std::string& query,
                                             const std::string& target,
                                             int maxSubstitution,
                                             int maxInsertion,
                                             int maxDeletion,
                                             int maxEdits) const {
    int qLen = static_cast<int>(query.size());
    int tLen = static_cast<int>(target.size());
    std::vector<std::vector<std::vector<EditCounts>>> dp(qLen + 1, std::vector<std::vector<EditCounts>>(tLen + 1));
    dp[0][0].push_back({0, 0, 0});

    for (int i = 0; i <= qLen; ++i) {
        for (int j = 0; j <= tLen; ++j) {
            auto states = dp[i][j];
            for (const auto& state : states) {
                if (i < qLen) {
                    EditCounts next = state;
                    ++next.del;
                    if (next.del <= maxDeletion && next.total() <= maxEdits) AddState(dp[i + 1][j], next);
                }
                if (j < tLen) {
                    EditCounts next = state;
                    ++next.ins;
                    if (next.ins <= maxInsertion && next.total() <= maxEdits) AddState(dp[i][j + 1], next);
                }
                if (i < qLen && j < tLen) {
                    EditCounts next = state;
                    if (query[i] != target[j]) ++next.sub;
                    if (next.sub <= maxSubstitution && next.total() <= maxEdits) AddState(dp[i + 1][j + 1], next);
                }
            }
        }
    }

    int best = std::numeric_limits<int>::max();
    for (const auto& state : dp[qLen][tLen]) {
        if (state.sub <= maxSubstitution && state.ins <= maxInsertion && state.del <= maxDeletion && state.total() <= maxEdits) {
            best = std::min(best, state.total());
        }
    }
    if (best == std::numeric_limits<int>::max()) return std::nullopt;
    return best;
}

std::vector<std::pair<std::size_t, float>> Trie::SearchIndicesWithMatrix(
    const std::string& query,
    float maxCost,
    const std::optional<std::string>& vGeneFilter,
    const std::optional<std::string>& jGeneFilter) const {

    if (!matrix_.Loaded()) throw std::runtime_error("No substitution matrix is loaded");

    std::vector<float> initialRow(query.size() + 1, 0.0f);
    for (std::size_t i = 1; i <= query.size(); ++i) {
        initialRow[i] = initialRow[i - 1] + matrix_.Cost(query[i - 1], '-');
    }

    std::vector<std::pair<std::size_t, float>> results;
    SearchMatrixRecursive(*root_, query, maxCost, initialRow, results, vGeneFilter, jGeneFilter);
    return results;
}

void Trie::SearchMatrixRecursive(const Node& node,
                                 const std::string& query,
                                 float maxCost,
                                 const std::vector<float>& prevRow,
                                 std::vector<std::pair<std::size_t, float>>& results,
                                 const std::optional<std::string>& vGeneFilter,
                                 const std::optional<std::string>& jGeneFilter) const {
    if (!node.indices.empty() && prevRow.back() <= maxCost) {
        for (int index : node.indices) {
            if (vGeneFilter && vGenes_[index] != *vGeneFilter) continue;
            if (jGeneFilter && jGenes_[index] != *jGeneFilter) continue;
            results.emplace_back(static_cast<std::size_t>(index), prevRow.back());
        }
    }

    for (const auto& childPtr : node.children) {
        if (!childPtr) continue;
        const Node& child = *childPtr;
        std::vector<float> nextRow(query.size() + 1, 0.0f);
        nextRow[0] = prevRow[0] + matrix_.Cost('-', child.letter);
        float minValue = nextRow[0];

        for (std::size_t j = 1; j <= query.size(); ++j) {
            char queryChar = query[j - 1];
            float consumeTarget = prevRow[j] + matrix_.Cost('-', child.letter);
            float consumeQuery = nextRow[j - 1] + matrix_.Cost(queryChar, '-');
            float consumeBoth = prevRow[j - 1] + matrix_.Cost(queryChar, child.letter);
            nextRow[j] = std::min({consumeTarget, consumeQuery, consumeBoth});
            minValue = std::min(minValue, nextRow[j]);
        }

        if (minValue <= maxCost) {
            SearchMatrixRecursive(child, query, maxCost, nextRow, results, vGeneFilter, jGeneFilter);
        }
    }
}

std::optional<AlignmentResult> Trie::AlignIndexHit(const std::string& query,
                                                   std::size_t targetIndex,
                                                   int maxSubstitution,
                                                   int maxInsertion,
                                                   int maxDeletion,
                                                   int maxEdits) const {
    if (targetIndex >= sequences_.size()) return std::nullopt;
    auto distance = BoundedEditDistance(query, sequences_[targetIndex], maxSubstitution, maxInsertion, maxDeletion, maxEdits);
    if (!distance) return std::nullopt;
    AlignmentResult result;
    result.distance = static_cast<float>(*distance);
    return result;
}

std::optional<AlignmentResult> Trie::AlignIndexHitWithMatrix(const std::string& query,
                                                             std::size_t targetIndex,
                                                             std::optional<float> maxCost) const {
    if (!matrix_.Loaded()) throw std::runtime_error("No substitution matrix is loaded");
    if (targetIndex >= sequences_.size()) return std::nullopt;

    const std::string& target = sequences_[targetIndex];
    int qLen = static_cast<int>(query.size());
    int tLen = static_cast<int>(target.size());
    std::vector<std::vector<MatrixCell>> dp(qLen + 1, std::vector<MatrixCell>(tLen + 1));
    dp[0][0].reachable = true;
    dp[0][0].cost = 0.0f;

    for (int i = 0; i <= qLen; ++i) {
        for (int j = 0; j <= tLen; ++j) {
            if (!dp[i][j].reachable) continue;
            if (i < qLen) {
                float nextCost = dp[i][j].cost + matrix_.Cost(query[i], '-');
                if (Better(nextCost, AlignmentOpType::Deletion, dp[i + 1][j])) {
                    dp[i + 1][j] = {nextCost, i, j, AlignmentOpType::Deletion, true};
                }
            }
            if (j < tLen) {
                float nextCost = dp[i][j].cost + matrix_.Cost('-', target[j]);
                if (Better(nextCost, AlignmentOpType::Insertion, dp[i][j + 1])) {
                    dp[i][j + 1] = {nextCost, i, j, AlignmentOpType::Insertion, true};
                }
            }
            if (i < qLen && j < tLen) {
                auto op = query[i] == target[j] ? AlignmentOpType::Match : AlignmentOpType::Substitution;
                float nextCost = dp[i][j].cost + matrix_.Cost(query[i], target[j]);
                if (Better(nextCost, op, dp[i + 1][j + 1])) {
                    dp[i + 1][j + 1] = {nextCost, i, j, op, true};
                }
            }
        }
    }

    if (!dp[qLen][tLen].reachable) return std::nullopt;
    if (maxCost && dp[qLen][tLen].cost > *maxCost + kEps) return std::nullopt;

    AlignmentResult result;
    result.distance = dp[qLen][tLen].cost;
    std::string queryRev;
    std::string targetRev;
    int i = qLen;
    int j = tLen;

    while (i > 0 || j > 0) {
        MatrixCell cell = dp[i][j];
        if (cell.op == AlignmentOpType::Deletion) {
            ++result.deletions;
            queryRev.push_back(query[i - 1]);
            targetRev.push_back('-');
        } else if (cell.op == AlignmentOpType::Insertion) {
            ++result.insertions;
            queryRev.push_back('-');
            targetRev.push_back(target[j - 1]);
        } else {
            if (cell.op == AlignmentOpType::Substitution) ++result.substitutions;
            queryRev.push_back(query[i - 1]);
            targetRev.push_back(target[j - 1]);
        }
        i = cell.prevI;
        j = cell.prevJ;
    }

    std::reverse(queryRev.begin(), queryRev.end());
    std::reverse(targetRev.begin(), targetRev.end());
    result.queryAligned = std::move(queryRev);
    result.targetAligned = std::move(targetRev);
    return result;
}
