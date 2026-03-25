#include "MatchWriter.h"

#include <stdexcept>

void MatchWriterQueue::Push(std::vector<std::string> batch) {
    {
        std::lock_guard<std::mutex> lock(mutex_);
        if (closed_) {
            throw std::runtime_error("Attempt to push to closed MatchWriterQueue");
        }
        queue_.push_back(std::move(batch));
    }
    cv_.notify_one();
}

bool MatchWriterQueue::Pop(std::vector<std::string>& batch) {
    std::unique_lock<std::mutex> lock(mutex_);
    cv_.wait(lock, [&]() { return closed_ || !queue_.empty(); });
    if (queue_.empty()) {
        return false;
    }
    batch = std::move(queue_.front());
    queue_.pop_front();
    return true;
}

void MatchWriterQueue::Close() {
    {
        std::lock_guard<std::mutex> lock(mutex_);
        closed_ = true;
    }
    cv_.notify_all();
}

MatchWriter::MatchWriter(const std::string& path, bool withAlignment)
    : out_(path), withAlignment_(withAlignment) {
    if (!out_) {
        throw std::runtime_error("Failed to open output file: " + path);
    }
}

void MatchWriter::WriteHeader() {
    out_ << "query_row\tquery_cdr3\tquery_v\tquery_j\tquery_epitope\tquery_species\tquery_chain"
         << "\ttarget_row\ttarget_cdr3\ttarget_v\ttarget_j\ttarget_epitope\ttarget_species\ttarget_chain"
         << "\tdistance\tsubstitutions\tinsertions\tdeletions\tscore";
    if (withAlignment_) {
        out_ << "\tquery_aligned\ttarget_aligned";
    }
    out_ << '\n';
}

void MatchWriter::WriteBatch(const std::vector<std::string>& batch) {
    for (const auto& line : batch) {
        out_ << line << '\n';
    }
}
