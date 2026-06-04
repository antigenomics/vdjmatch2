#include "RepertoireMatcher.h"

#include "Alignment.h"
#include "MatchWriter.h"
#include "RegionTrie.h"
#include "Trie.h"
#include "TsvReader.h"

#include <atomic>
#include <mutex>
#include <thread>
#include <exception>
#include <algorithm>
#include <chrono>
#include <iomanip>
#include <iostream>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
    std::string ToStringFloat(float value) {
        std::ostringstream out;
        out << std::fixed << std::setprecision(6) << value;
        return out.str();
    }

    std::string MaybeInt(const std::optional<AlignmentResult>& alignment, int value) {
        if (!alignment) return "";
        return std::to_string(value);
    }

    std::string BuildLine(const Record& query,
                          const Record& target,
                          const std::optional<AlignmentResult>& alignment,
                          bool withAlignment,
                          OutputMode mode,
                          float score) {
        std::ostringstream out;
        out << query.rowIndex << '\t'
            << query.junctionAA << '\t'
            << query.vGene << '\t'
            << query.jGene << '\t'
            << query.epitope << '\t'
            << query.species << '\t'
            << query.chain << '\t'
            << target.rowIndex << '\t'
            << target.junctionAA << '\t'
            << target.vGene << '\t'
            << target.jGene << '\t'
            << target.epitope << '\t'
            << target.species << '\t'
            << target.chain << '\t';

        if (alignment) {
            out << ToStringFloat(alignment->distance) << '\t'
                << alignment->substitutions << '\t'
                << alignment->insertions << '\t'
                << alignment->deletions << '\t'
                << ToStringFloat(score);
        } else {
            out << ToStringFloat(score) << "\t\t\t\t" << ToStringFloat(score);
        }

        if (mode == OutputMode::RegionalMatrix) {
            if (alignment) {
                out << '\t' << ToStringFloat(alignment->regionCosts[RegionIndex(Region::V)])
                    << '\t' << ToStringFloat(alignment->regionCosts[RegionIndex(Region::NDN)])
                    << '\t' << ToStringFloat(alignment->regionCosts[RegionIndex(Region::J)]);
            } else {
                out << "\t\t\t";
            }
        }

        if (withAlignment) {
            out << '\t';
            if (alignment) out << alignment->queryAligned;
            out << '\t';
            if (alignment) out << alignment->targetAligned;
            if (mode == OutputMode::RegionalMatrix) {
                out << '\t';
                if (alignment) out << alignment->queryRegionsAligned;
                out << '\t';
                if (alignment) out << alignment->targetRegionsAligned;
            }
        }

        return out.str();
    }

    void PrintReadStats(const std::string& title, const TsvReadResult& result) {
        std::cerr << title
                  << " records=" << result.records.size()
                  << " skipped_missing_cdr3=" << result.skippedMissingJunction
                  << " skipped_invalid_cdr3=" << result.skippedInvalidSequence
                  << " skipped_missing_region=" << result.skippedMissingRegion
                  << " skipped_invalid_region=" << result.skippedInvalidRegion
                  << " skipped_by_filter=" << result.skippedByFilter
                  << '\n';
    }

    std::optional<std::string> BuildVFilter(const CliConfig& config, const Record& query) {
        if (config.matchV && !query.vGene.empty()) return query.vGene;
        return std::nullopt;
    }

    std::optional<std::string> BuildJFilter(const CliConfig& config, const Record& query) {
        if (config.matchJ && !query.jGene.empty()) return query.jGene;
        return std::nullopt;
    }
}

RepertoireMatcher::RepertoireMatcher(CliConfig config) : config_(std::move(config)) {}

int RepertoireMatcher::Run() {
    auto start = std::chrono::high_resolution_clock::now();

    const auto queryResult = ReadTsv(config_.queryPath, config_);
    const auto targetResult = ReadTsv(config_.targetPath, config_);

    PrintReadStats("query", queryResult);
    PrintReadStats("target", targetResult);

    if (queryResult.records.empty()) throw std::runtime_error("No query records remained after filtering");
    if (targetResult.records.empty()) throw std::runtime_error("No target records remained after filtering");

    std::vector<std::string> targetSequences;
    std::vector<std::vector<Region>> targetRegions;
    std::vector<std::string> targetVGenes;
    std::vector<std::string> targetJGenes;

    targetSequences.reserve(targetResult.records.size());
    targetRegions.reserve(targetResult.records.size());
    targetVGenes.reserve(targetResult.records.size());
    targetJGenes.reserve(targetResult.records.size());

    for (const auto& record : targetResult.records) {
        targetSequences.push_back(record.junctionAA);
        targetRegions.push_back(record.regions);
        targetVGenes.push_back(record.vGene);
        targetJGenes.push_back(record.jGene);
    }

    OutputMode outputMode = OutputMode::Edit;
    if (config_.IsRegionalMatrixMode()) outputMode = OutputMode::RegionalMatrix;
    else if (config_.IsMatrixMode()) outputMode = OutputMode::Matrix;

    MatchWriterQueue queue;
    MatchWriter writer(config_.outPath, config_.writeAlignment, outputMode);
    writer.WriteHeader();

    std::exception_ptr firstException;
    std::mutex exceptionMutex;
    std::atomic<bool> stopRequested{false};

    auto storeException = [&](std::exception_ptr ex) {
        {
            std::lock_guard<std::mutex> lock(exceptionMutex);
            if (!firstException) {
                firstException = ex;
            }
        }
        stopRequested.store(true);
        queue.Close();
    };

    std::thread writerThread([&]() {
        try {
            std::vector<std::string> batch;
            while (!stopRequested.load()) {
                if (!queue.Pop(batch)) {
                    break;
                }
                writer.WriteBatch(batch);
                batch.clear();
            }
        } catch (...) {
            storeException(std::current_exception());
        }
    });

    constexpr std::size_t kQueryBatchSize = 1000;

    std::atomic<std::size_t> next{0};
    std::atomic<std::size_t> written{0};
    std::vector<std::thread> workers;
    workers.reserve(static_cast<std::size_t>(config_.threads));

    auto workerLoop = [&](auto& trieInstance) {
        return [&]() {
            try {
                while (!stopRequested.load()) {
                    const std::size_t batchBegin = next.fetch_add(kQueryBatchSize);
                    if (batchBegin >= queryResult.records.size()) {
                        break;
                    }

                    const std::size_t batchEnd =
                            std::min(batchBegin + kQueryBatchSize, queryResult.records.size());

                    std::vector<std::string> outputBatch;

                    for (std::size_t idx = batchBegin; idx < batchEnd; ++idx) {
                        if (stopRequested.load()) break;

                        const auto& query = queryResult.records[idx];

                        if constexpr (std::is_same_v<std::decay_t<decltype(trieInstance)>, RegionTrie>) {
                            auto hits = trieInstance.SearchIndices(query.junctionAA, query.regions, config_.maxCost, BuildVFilter(config_, query), BuildJFilter(config_, query));
                            outputBatch.reserve(outputBatch.size() + hits.size());
                            for (const auto& [targetIndex, cost] : hits) {
                                std::optional<AlignmentResult> alignment = std::nullopt;
                                if (config_.writeAlignment) {
                                    alignment = trieInstance.AlignIndexHit(query.junctionAA, query.regions, targetIndex, config_.maxCost);
                                }
                                outputBatch.push_back(BuildLine(query, targetResult.records[targetIndex], alignment, config_.writeAlignment, outputMode, cost));
                            }
                        } else {
                            if (config_.IsMatrixMode()) {
                                auto hits = trieInstance.SearchIndicesWithMatrix(query.junctionAA, config_.maxCost, BuildVFilter(config_, query), BuildJFilter(config_, query));
                                outputBatch.reserve(outputBatch.size() + hits.size());
                                for (const auto& [targetIndex, cost] : hits) {
                                    std::optional<AlignmentResult> alignment = std::nullopt;
                                    if (config_.writeAlignment) alignment = trieInstance.AlignIndexHitWithMatrix(query.junctionAA, targetIndex, config_.maxCost);
                                    outputBatch.push_back(BuildLine(query, targetResult.records[targetIndex], alignment, config_.writeAlignment, outputMode, cost));
                                }
                            } else {
                                auto hits = trieInstance.SearchIndices(query.junctionAA, config_.maxSub, config_.maxIns, config_.maxDel, config_.maxEdits, BuildVFilter(config_, query), BuildJFilter(config_, query));
                                outputBatch.reserve(outputBatch.size() + hits.size());
                                for (const auto& [targetIndex, distance] : hits) {
                                    std::optional<AlignmentResult> alignment = std::nullopt;
                                    if (config_.writeAlignment) alignment = trieInstance.AlignIndexHit(query.junctionAA, targetIndex, config_.maxSub, config_.maxIns, config_.maxDel, config_.maxEdits);
                                    outputBatch.push_back(BuildLine(query, targetResult.records[targetIndex], alignment, config_.writeAlignment, outputMode, static_cast<float>(distance)));
                                }
                            }
                        }
                    }

                    written.fetch_add(outputBatch.size());
                    if (!outputBatch.empty() && !stopRequested.load()) {
                        queue.Push(std::move(outputBatch));
                    }
                }
            } catch (...) {
                storeException(std::current_exception());
            }
        };
    };

    if (config_.IsRegionalMatrixMode()) {
        RegionTrie trie(targetSequences, targetRegions, targetVGenes, targetJGenes);
        trie.LoadRegionMatrices(*config_.matrixVPath, *config_.matrixNDNPath, *config_.matrixJPath);
        for (int t = 0; t < config_.threads; ++t) {
            workers.emplace_back(workerLoop(trie));
        }
        for (auto& worker : workers) {
            if (worker.joinable()) worker.join();
        }
    } else {
        Trie trie(targetSequences, targetVGenes, targetJGenes);
        if (config_.IsMatrixMode()) trie.LoadSubstitutionMatrix(*config_.matrixPath);
        for (int t = 0; t < config_.threads; ++t) {
            workers.emplace_back(workerLoop(trie));
        }
        for (auto& worker : workers) {
            if (worker.joinable()) worker.join();
        }
    }

    queue.Close();

    if (writerThread.joinable()) {
        writerThread.join();
    }

    if (firstException) {
        std::rethrow_exception(firstException);
    }

    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

    std::cout << "query_records=" << queryResult.records.size()
              << " target_records=" << targetResult.records.size()
              << " matches_written=" << written.load() << '\n';
    std::cout << "Time: " << duration.count() << " ms" << std::endl;
    return 0;
}
