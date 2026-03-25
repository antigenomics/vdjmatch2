#include "RepertoireMatcher.h"

#include "MatchWriter.h"
#include "Trie.h"
#include "TsvReader.h"

#include <atomic>
#include <chrono>
#include <iomanip>
#include <iostream>
#include <optional>
#include <sstream>
#include <thread>
#include <vector>

namespace {
std::string EscapeField(const std::string& value) {
    return value;
}

std::string ToStringFloat(float value) {
    std::ostringstream out;
    out << std::fixed << std::setprecision(6) << value;
    return out.str();
}

std::string BuildLine(const Record& query,
                      const Record& target,
                      const std::optional<Trie::AlignmentResult>& alignment,
                      bool withAlignment,
                      bool matrixMode,
                      float scoreOrDistance) {
    std::ostringstream out;
    out << query.rowIndex << '\t'
        << EscapeField(query.junctionAA) << '\t'
        << EscapeField(query.vGene) << '\t'
        << EscapeField(query.jGene) << '\t'
        << EscapeField(query.epitope) << '\t'
        << EscapeField(query.species) << '\t'
        << EscapeField(query.chain) << '\t'
        << target.rowIndex << '\t'
        << EscapeField(target.junctionAA) << '\t'
        << EscapeField(target.vGene) << '\t'
        << EscapeField(target.jGene) << '\t'
        << EscapeField(target.epitope) << '\t'
        << EscapeField(target.species) << '\t'
        << EscapeField(target.chain) << '\t';

    if (alignment) {
        out << ToStringFloat(alignment->distance) << '\t'
            << alignment->substitutions << '\t'
            << alignment->insertions << '\t'
            << alignment->deletions << '\t';
    } else {
        out << "\t\t\t\t";
    }

    if (matrixMode) {
        out << ToStringFloat(scoreOrDistance);
    } else {
        out << ToStringFloat(scoreOrDistance);
    }

    if (withAlignment) {
        out << '\t';
        if (alignment) out << alignment->queryAligned;
        out << '\t';
        if (alignment) out << alignment->targetAligned;
    }

    return out.str();
}
}

RepertoireMatcher::RepertoireMatcher(CliConfig config) : config_(std::move(config)) {}

int RepertoireMatcher::Run() {
    const auto queryResult = ReadTsv(config_.queryPath, config_);
    const auto targetResult = ReadTsv(config_.targetPath, config_);

    if (queryResult.records.empty()) {
        throw std::runtime_error("No query records remained after filtering");
    }
    if (targetResult.records.empty()) {
        throw std::runtime_error("No target records remained after filtering");
    }

    std::vector<std::string> sequences;
    std::vector<std::string> vGenes;
    std::vector<std::string> jGenes;
    sequences.reserve(targetResult.records.size());
    vGenes.reserve(targetResult.records.size());
    jGenes.reserve(targetResult.records.size());
    for (const auto& record : targetResult.records) {
        sequences.push_back(record.junctionAA);
        vGenes.push_back(record.vGene);
        jGenes.push_back(record.jGene);
    }

    auto start = std::chrono::high_resolution_clock::now();
    Trie trie(sequences, vGenes, jGenes);
    const bool matrixMode = config_.matrixPath.has_value();
    if (matrixMode) {
        trie.LoadSubstitutionMatrix(*config_.matrixPath);
    }

    MatchWriterQueue queue;
    MatchWriter writer(config_.outPath, config_.writeAlignment);
    writer.WriteHeader();

    std::thread writerThread([&]() {
        std::vector<std::string> batch;
        while (queue.Pop(batch)) {
            writer.WriteBatch(batch);
            batch.clear();
        }
    });

    constexpr std::size_t kQueryBatchSize = 1000;

    std::atomic<std::size_t> next{0};
    std::atomic<std::size_t> written{0};
    std::vector<std::thread> workers;
    workers.reserve(static_cast<std::size_t>(config_.threads));

    for (int t = 0; t < config_.threads; ++t) {
        workers.emplace_back([&]() {
            while (true) {
                const std::size_t batchBegin = next.fetch_add(kQueryBatchSize);
                if (batchBegin >= queryResult.records.size()) {
                    break;
                }

                const std::size_t batchEnd = std::min(batchBegin + kQueryBatchSize, queryResult.records.size());
                std::vector<std::string> outputBatch;

                for (std::size_t idx = batchBegin; idx < batchEnd; ++idx) {
                    const auto& query = queryResult.records[idx];
                    const auto vFilter = config_.matchV ? std::optional<std::string>(query.vGene) : std::nullopt;
                    const auto jFilter = config_.matchJ ? std::optional<std::string>(query.jGene) : std::nullopt;

                    if (matrixMode) {
                        auto hits = trie.SearchIndicesWithMatrix(query.junctionAA, config_.maxCost, vFilter, jFilter);
                        outputBatch.reserve(outputBatch.size() + hits.size());
                        for (const auto& [targetIndex, cost] : hits) {
                            std::optional<Trie::AlignmentResult> alignment = std::nullopt;
                            if (config_.writeAlignment) {
                                alignment = trie.AlignIndexHitWithMatrix(query.junctionAA, targetIndex, config_.maxCost);
                            }
                            outputBatch.push_back(BuildLine(query,
                                                            targetResult.records[targetIndex],
                                                            alignment,
                                                            config_.writeAlignment,
                                                            true,
                                                            cost));
                        }
                    } else {
                        auto hits = trie.SearchIndices(query.junctionAA,
                                                       config_.maxSub,
                                                       config_.maxIns,
                                                       config_.maxDel,
                                                       config_.maxEdits,
                                                       vFilter,
                                                       jFilter);
                        outputBatch.reserve(outputBatch.size() + hits.size());
                        for (const auto& [targetIndex, distance] : hits) {
                            std::optional<Trie::AlignmentResult> alignment = std::nullopt;
                            if (config_.writeAlignment) {
                                alignment = trie.AlignIndexHit(query.junctionAA,
                                                               targetIndex,
                                                               config_.maxSub,
                                                               config_.maxIns,
                                                               config_.maxDel,
                                                               config_.maxEdits);
                            }
                            outputBatch.push_back(BuildLine(query,
                                                            targetResult.records[targetIndex],
                                                            alignment,
                                                            config_.writeAlignment,
                                                            false,
                                                            static_cast<float>(distance)));
                        }
                    }
                }

                written.fetch_add(outputBatch.size());
                if (!outputBatch.empty()) {
                    queue.Push(std::move(outputBatch));
                }
            }
        });
    }

    for (auto& worker : workers) {
        worker.join();
    }
    queue.Close();
    writerThread.join();

    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::cout << "query_records=" << queryResult.records.size()
              << " target_records=" << targetResult.records.size()
              << " matches_written=" << written.load() << '\n';

    std::cout << "Time: " << duration.count() << " ms" << std::endl;
    return 0;
}
