#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

class Trie {
public:
    static constexpr int kAlphabetSize = 26;
    static constexpr int kMaxPareto = 64;

    struct TrieNode {
        std::array<TrieNode*, kAlphabetSize> children{};
        std::vector<int> indices;
    };

    struct EditState {
        int16_t sub;
        int16_t ins;
        int16_t del;
        int total() const { return sub + ins + del; }
    };

    struct StatCell {
        EditState data[kMaxPareto];
        int size = 0;

        void clear() { size = 0; }
        bool empty() const { return size == 0; }
        void push_back(EditState s) {
            if (size < kMaxPareto) data[size++] = s;
        }
    };

    enum class AlignmentOpType {
        Match,
        Substitution,
        Insertion,
        Deletion
    };

    struct AlignmentOp {
        AlignmentOpType type;
        int queryPos;
        char queryChar;
        char targetChar;
    };

    struct AlignmentResult {
        std::string queryAligned;
        std::string targetAligned;
        int substitutions = 0;
        int insertions = 0;
        int deletions = 0;
        float distance = 0.0f;
        std::vector<AlignmentOp> ops;
    };

    explicit Trie();
    explicit Trie(const std::vector<std::string>& sequences,
                  const std::vector<std::string>& vGenes,
                  const std::vector<std::string>& jGenes);

    Trie(const Trie& other);
    Trie& operator=(const Trie& other);
    Trie(Trie&& other) noexcept;
    Trie& operator=(Trie&& other) noexcept;
    ~Trie();

    void LoadSubstitutionMatrix(const std::string& matrixPath,
                                const std::string& delimiter = "",
                                float gapFactor = 1.5f);

    std::vector<std::pair<size_t, int>> SearchIndices(
        const std::string& query,
        int maxSubstitution = 0,
        int maxInsertion = 0,
        int maxDeletion = 0,
        std::optional<int> maxEdits = std::nullopt,
        const std::optional<std::string>& vGeneFilter = std::nullopt,
        const std::optional<std::string>& jGeneFilter = std::nullopt);

    std::vector<std::pair<size_t, float>> SearchIndicesWithMatrix(
        const std::string& query,
        float maxCost,
        const std::optional<std::string>& vGeneFilter = std::nullopt,
        const std::optional<std::string>& jGeneFilter = std::nullopt);

    std::optional<AlignmentResult> AlignQueryToTarget(
        const std::string& query,
        const std::string& target,
        std::optional<int> maxSubstitution = std::nullopt,
        std::optional<int> maxInsertion = std::nullopt,
        std::optional<int> maxDeletion = std::nullopt,
        std::optional<int> maxEdits = std::nullopt);

    std::optional<AlignmentResult> AlignIndexHit(
        const std::string& query,
        size_t targetIndex,
        std::optional<int> maxSubstitution = std::nullopt,
        std::optional<int> maxInsertion = std::nullopt,
        std::optional<int> maxDeletion = std::nullopt,
        std::optional<int> maxEdits = std::nullopt);

    std::optional<AlignmentResult> AlignQueryToTargetWithMatrix(
        const std::string& query,
        const std::string& target,
        std::optional<float> maxCost = std::nullopt);

    std::optional<AlignmentResult> AlignIndexHitWithMatrix(
        const std::string& query,
        size_t targetIndex,
        std::optional<float> maxCost = std::nullopt);

private:
    bool useSubstitutionMatrix_ = false;
    std::unordered_map<char, std::unordered_map<char, float>> substitutionMatrix_;
    TrieNode* root_ = nullptr;
    std::vector<std::string> sequences_;
    std::vector<std::string> vGenes_;
    std::vector<std::string> jGenes_;

    void BuildTrie();
    void DeleteTrie(TrieNode* node);
    TrieNode* CopyTrie(const TrieNode* node);

    void SearchSubstitutionOnlyIDs(
        const std::string& query,
        int maxSub,
        TrieNode* node,
        int depth,
        int mismatches,
        int queryLength,
        std::vector<std::pair<size_t, int>>& results,
        const std::optional<std::string>& vGeneFilter,
        const std::optional<std::string>& jGeneFilter);

    void SearchRecursiveIDs(
        const std::string& query,
        int maxEdits,
        TrieNode* node,
        const int* prevRow,
        int queryLength,
        std::vector<std::pair<size_t, int>>& results,
        const std::optional<std::string>& vGeneFilter,
        const std::optional<std::string>& jGeneFilter);

    template <typename ResultType, typename EmitFunc>
    void SearchRecursiveDetailed(
        const std::string& query,
        int maxEdits,
        int maxSub,
        int maxIns,
        int maxDel,
        TrieNode* node,
        const int* prevRowSimple,
        const StatCell* prevRowDetailed,
        int queryLength,
        std::vector<ResultType>& results,
        const std::optional<std::string>& vGeneFilter,
        const std::optional<std::string>& jGeneFilter,
        EmitFunc emitFunc);

    static void PrunePareto(StatCell& cell);

    void SearchRecursiveCostIDs(
        const std::string& query,
        float maxCost,
        TrieNode* node,
        const float* prevRow,
        int queryLength,
        std::vector<std::pair<size_t, float>>& results,
        const std::optional<std::string>& vGeneFilter,
        const std::optional<std::string>& jGeneFilter);
};
