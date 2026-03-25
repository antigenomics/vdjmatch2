#include "CliConfig.h"
#include "RepertoireMatcher.h"
#include "TsvReader.h"

#include <filesystem>
#include <fstream>
#include <functional>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace fs = std::filesystem;

struct TestFailure : std::runtime_error {
    using std::runtime_error::runtime_error;
};

#define EXPECT_TRUE(cond) do { if (!(cond)) throw TestFailure(std::string("EXPECT_TRUE failed: ") + #cond + " at " + __FILE__ + ":" + std::to_string(__LINE__)); } while (0)
#define EXPECT_EQ(a, b) do { auto _a = (a); auto _b = (b); if (!((_a) == (_b))) { std::ostringstream _oss; _oss << "EXPECT_EQ failed: " << #a << " != " << #b << " at " << __FILE__ << ":" << __LINE__ << " (" << _a << " vs " << _b << ")"; throw TestFailure(_oss.str()); } } while (0)
#define EXPECT_NE(a, b) do { auto _a = (a); auto _b = (b); if (((_a) == (_b))) { std::ostringstream _oss; _oss << "EXPECT_NE failed: " << #a << " == " << #b << " at " << __FILE__ << ":" << __LINE__; throw TestFailure(_oss.str()); } } while (0)

static fs::path MakeTempDir(const std::string& name) {
    auto dir = fs::temp_directory_path() / ("vdjmatch2_tests_" + name);
    fs::remove_all(dir);
    fs::create_directories(dir);
    return dir;
}

static void WriteFile(const fs::path& path, const std::string& content) {
    std::ofstream out(path);
    if (!out) throw std::runtime_error("Failed to open file for writing: " + path.string());
    out << content;
}

static std::vector<std::string> ReadLines(const fs::path& path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("Failed to open file for reading: " + path.string());
    std::vector<std::string> lines;
    std::string line;
    while (std::getline(in, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        lines.push_back(line);
    }
    return lines;
}

static void TestParseCliDefaults() {
    const char* argv[] = {
        "vdjmatch2", "query.tsv", "target.tsv"
    };
    CliConfig cfg = ParseCli(3, const_cast<char**>(argv));
    EXPECT_EQ(cfg.queryPath, std::string("query.tsv"));
    EXPECT_EQ(cfg.targetPath, std::string("target.tsv"));
    EXPECT_EQ(cfg.outPath, std::string("match_result.tsv"));
    EXPECT_TRUE(cfg.gene.empty());
    EXPECT_EQ(cfg.species, std::string("HomoSapiens"));
    EXPECT_TRUE(!cfg.matrixPath.has_value());
    EXPECT_EQ(cfg.maxEdits, -1);
}

static void TestReadTsvOptionalColumns() {
    auto dir = MakeTempDir("optional_columns");
    auto path = dir / "input.tsv";
    WriteFile(path,
              "junction_aa\tspecies\n"
              "cassl\tHomoSapiens\n");

    CliConfig cfg;
    cfg.species = "HomoSapiens";
    cfg.gene = "";
    cfg.epitope = std::nullopt;

    auto result = ReadTsv(path.string(), cfg);
    EXPECT_EQ(result.records.size(), static_cast<std::size_t>(1));
    EXPECT_EQ(result.records[0].junctionAA, std::string("CASSL"));
    EXPECT_TRUE(result.records[0].epitope.empty());
    EXPECT_TRUE(result.records[0].chain.empty());
}

static void TestReadTsvPreFiltersApplied() {
    auto dir = MakeTempDir("prefilters");
    auto path = dir / "input.tsv";
    WriteFile(path,
              "junction_aa\tv_call\tj_call\tantigen.epitope\tspecies\tgene\n"
              "CASSL\tTRBV1\tTRBJ1\tEPI\tHomoSapiens\tTRB\n"
              "CASST\tTRBV2\tTRBJ2\tEPI\tMusMusculus\tTRB\n"
              "CASSQ\tTRAV1\tTRAJ1\tEPI\tHomoSapiens\tTRA\n"
              "CASSP\tTRBV3\tTRBJ3\tOTHER\tHomoSapiens\tTRB\n");

    CliConfig cfg;
    cfg.species = "HomoSapiens";
    cfg.gene = "TRB";
    cfg.epitope = std::string("EPI");

    auto result = ReadTsv(path.string(), cfg);
    EXPECT_EQ(result.records.size(), static_cast<std::size_t>(1));
    EXPECT_EQ(result.records[0].junctionAA, std::string("CASSL"));
    EXPECT_EQ(result.skippedByFilter, static_cast<std::size_t>(3));
}

static void TestRepertoireMatcherEditModeWithAlignment() {
    auto dir = MakeTempDir("edit_mode");
    auto queryPath = dir / "query.tsv";
    auto targetPath = dir / "target.tsv";
    auto outPath = dir / "out.tsv";

    WriteFile(queryPath,
              "junction_aa\tv_call\tj_call\tspecies\tgene\n"
              "CASSL\tTRBV1\tTRBJ1\tHomoSapiens\tTRB\n"
              "AAAAA\tTRBV9\tTRBJ9\tHomoSapiens\tTRB\n");
    WriteFile(targetPath,
              "junction_aa\tv_call\tj_call\tspecies\tgene\n"
              "CASSQ\tTRBV1\tTRBJ1\tHomoSapiens\tTRB\n"
              "AAA\tTRBV8\tTRBJ8\tHomoSapiens\tTRB\n");

    CliConfig cfg;
    cfg.queryPath = queryPath.string();
    cfg.targetPath = targetPath.string();
    cfg.outPath = outPath.string();
    cfg.species = "HomoSapiens";
    cfg.gene = "TRB";
    cfg.writeAlignment = true;
    cfg.threads = 2;
    cfg.maxSub = 0;
    cfg.maxIns = 0;
    cfg.maxDel = 2;
    cfg.maxEdits = 2;

    RepertoireMatcher matcher(cfg);
    EXPECT_EQ(matcher.Run(), 0);

    auto lines = ReadLines(outPath);
    EXPECT_EQ(lines.size(), static_cast<std::size_t>(2));
    EXPECT_NE(lines[0].find("query_aligned"), std::string::npos);
    EXPECT_NE(lines[1].find("CASSL"), std::string::npos);
    EXPECT_NE(lines[1].find("CASSQ"), std::string::npos);
}

static void TestRepertoireMatcherMatchVFiltersPairing() {
    auto dir = MakeTempDir("match_v");
    auto queryPath = dir / "query.tsv";
    auto targetPath = dir / "target.tsv";
    auto outPath = dir / "out.tsv";

    WriteFile(queryPath,
              "junction_aa\tv_call\tj_call\tspecies\tgene\n"
              "CASSL\tTRBV1\tTRBJ1\tHomoSapiens\tTRB\n");
    WriteFile(targetPath,
              "junction_aa\tv_call\tj_call\tspecies\tgene\n"
              "CASSL\tTRBV1\tTRBJ9\tHomoSapiens\tTRB\n"
              "CASSL\tTRBV2\tTRBJ1\tHomoSapiens\tTRB\n");

    CliConfig cfg;
    cfg.queryPath = queryPath.string();
    cfg.targetPath = targetPath.string();
    cfg.outPath = outPath.string();
    cfg.species = "HomoSapiens";
    cfg.gene = "TRB";
    cfg.matchV = true;
    cfg.threads = 1;
    cfg.maxSub = 0;
    cfg.maxIns = 0;
    cfg.maxDel = 0;
    cfg.maxEdits = 0;

    RepertoireMatcher matcher(cfg);
    EXPECT_EQ(matcher.Run(), 0);

    auto lines = ReadLines(outPath);
    EXPECT_EQ(lines.size(), static_cast<std::size_t>(2));
    EXPECT_NE(lines[1].find("TRBV1"), std::string::npos);
    EXPECT_TRUE(lines[1].find("TRBV2") == std::string::npos);
}

int main() {
    const std::vector<std::pair<std::string, std::function<void()>>> tests = {
        {"ParseCliDefaults", TestParseCliDefaults},
        {"ReadTsvOptionalColumns", TestReadTsvOptionalColumns},
        {"ReadTsvPreFiltersApplied", TestReadTsvPreFiltersApplied},
        {"RepertoireMatcherEditModeWithAlignment", TestRepertoireMatcherEditModeWithAlignment},
        {"RepertoireMatcherMatchVFiltersPairing", TestRepertoireMatcherMatchVFiltersPairing},
    };

    int failed = 0;
    for (const auto& [name, test] : tests) {
        try {
            test();
            std::cout << "[PASS] " << name << '\n';
        } catch (const std::exception& e) {
            ++failed;
            std::cerr << "[FAIL] " << name << ": " << e.what() << '\n';
        }
    }

    if (failed != 0) {
        std::cerr << failed << " test(s) failed\n";
        return 1;
    }

    std::cout << "All tests passed\n";
    return 0;
}
