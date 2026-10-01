// Data-only fragment-body overlap counting. Python guards source identities.
#include <algorithm>
#include <array>
#include <charconv>
#include <chrono>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <fcntl.h>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <unistd.h>
#include <zlib.h>

using U64 = std::uint64_t;
constexpr U64 RECORD_BOUND = 1091589686;
constexpr std::size_t PEAKS = 306706;
constexpr const char* DUPLICATE_FAILURE = "duplicate fragment key within cell/chromosome/start/end (completion refused)";
void require(bool ok, const std::string& text) {
    if (!ok) throw std::runtime_error(text);
}
std::vector<std::string> fields(const std::string& line) {
    std::vector<std::string> out;
    std::size_t first = 0;
    while (true) {
        auto last = line.find('\t', first);
        out.push_back(line.substr(first, last == std::string::npos ? last : last - first));
        if (last == std::string::npos) return out;
        first = last + 1;
    }
}
U64 number(const std::string& value) {
    U64 out = 0;
    auto parsed = std::from_chars(value.data(), value.data() + value.size(), out);
    require(!value.empty() && parsed.ec == std::errc() &&
            parsed.ptr == value.data() + value.size(), "invalid unsigned integer");
    return out;
}
std::vector<std::string> contigs() {
    std::vector<std::string> out;
    for (int i = 1; i <= 22; ++i) out.push_back("chr" + std::to_string(i));
    out.insert(out.end(), {"chrX", "chrY"});
    std::sort(out.begin(), out.end());
    return out;
}
struct Peak { U64 start, end; std::size_t index; };
struct Index {
    std::vector<Peak> peaks;
    std::vector<U64> max_end;
    U64 build(std::size_t node, std::size_t lo, std::size_t hi) {
        if (hi - lo == 1) return max_end[node] = peaks[lo].end;
        auto mid = lo + (hi - lo) / 2;
        return max_end[node] = std::max(build(node * 2, lo, mid), build(node * 2 + 1, mid, hi));
    }
    void prepare() {
        std::sort(peaks.begin(), peaks.end(), [](const Peak& a, const Peak& b) {
            return a.start < b.start || (a.start == b.start && a.index < b.index);
        });
        max_end.resize(4 * peaks.size() + 4);
        if (!peaks.empty()) build(1, 0, peaks.size());
    }
    template<class F> void visit(std::size_t node, std::size_t lo, std::size_t hi,
                                std::size_t limit, U64 start, F& hit) const {
        if (lo >= limit || max_end[node] <= start) return;
        if (hi - lo == 1) { hit(peaks[lo].index); return; }
        auto mid = lo + (hi - lo) / 2;
        visit(node * 2, lo, mid, limit, start, hit);
        visit(node * 2 + 1, mid, hi, limit, start, hit);
    }
    template<class F> void overlaps(U64 start, U64 end, F& hit) const {
        auto stop = std::lower_bound(peaks.begin(), peaks.end(), end,
                                    [](const Peak& a, U64 b) { return a.start < b; });
        if (!peaks.empty()) visit(1, 0, peaks.size(), stop - peaks.begin(), start, hit);
    }
};
using Geometry = std::map<std::string, Index>;
Geometry geometry(const std::filesystem::path& path) {
    std::ifstream input(path);
    require(bool(input), "geometry input unavailable");
    std::string line;
    std::getline(input, line);
    if (!line.empty() && line.back() == '\r') line.pop_back();
    require(line == "peak_id\tchromosome\tsource_start_1based_closed\tsource_end_1based_closed\tbed_start_0based\tbed_end_half_open",
            "geometry header differs");
    Geometry out;
    std::unordered_set<std::string> ids;
    std::size_t count = 0;
    while (std::getline(input, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        auto row = fields(line);
        require(row.size() == 6, "geometry width differs");
        auto start = number(row[2]), end = number(row[3]);
        auto bed_start = number(row[4]), bed_end = number(row[5]);
        require(start > 0 && end >= start && bed_start == start - 1 && bed_end == end,
                "native geometry BED conversion differs");
        require((row[0] == row[1] + ":" + row[2] + "-" + row[3] ||
                 row[0] == row[1] + "-" + row[2] + "-" + row[3]) && ids.insert(row[0]).second,
                "native peak identifier differs or is duplicated");
        require(row[1] != "chrX" && row[1] != "chrY", "nonautosomal native target");
        auto chromosome = row[1].substr(0, 3) == "chr" ? number(row[1].substr(3)) : 0;
        require(chromosome >= 1 && chromosome <= 22, "native chromosome differs");
        out[row[1]].peaks.push_back({bed_start, bed_end, count++});
    }
    require(input.eof() && count == PEAKS && out.size() == 22, "native axis dimensions differ");
    for (auto& entry : out) entry.second.prepare();
    return out;
}
struct Cell { std::size_t index, half; };
using Cells = std::unordered_map<std::string, Cell>;
Cells cells(const std::filesystem::path& path) {
    std::ifstream input(path);
    require(bool(input), "cell identity input unavailable");
    std::string line;
    std::getline(input, line);
    require(line == "cell_id\tnucleus_half", "cell identity header differs");
    Cells out;
    while (std::getline(input, line)) {
        auto row = fields(line);
        require(row.size() == 2 && (row[1] == "A" || row[1] == "B"), "half identity differs");
        require(out.emplace(row[0], Cell{out.size(), row[1] == "A" ? 1UL : 2UL}).second,
                "duplicate cell identity");
    }
    require(input.eof() && !out.empty(), "empty or unreadable cell identity");
    return out;
}
struct DuplicateKey {
    U64 end; std::size_t cell;
    bool operator==(const DuplicateKey& other) const { return end == other.end && cell == other.cell; }
};
struct DuplicateHash {
    std::size_t operator()(const DuplicateKey& key) const {
        return std::hash<U64>{}(key.end) ^ (std::hash<std::size_t>{}(key.cell) + 0x9e3779b9);
    }
};
struct Result {
    std::array<std::vector<U64>, 3> counts;
    std::array<U64, 3> primary{}, autosomal{}, union_overlap{}, support{}, summed{};
    U64 records = 0;
    std::size_t seen_cells = 0, maximum_start_bucket_keys = 0;
    double elapsed_seconds = 0;
};
Result count(const std::filesystem::path& source, const Geometry& axis, const Cells& membership,
             std::size_t n_peaks, const std::string& tag) {
    const auto began = std::chrono::steady_clock::now();
    Result out;
    for (auto& values : out.counts) values.resize(n_peaks);
    auto chromosomes = contigs();
    std::unordered_map<std::string, int> rank;
    for (std::size_t i = 0; i < chromosomes.size(); ++i) rank[chromosomes[i]] = int(i);
    std::vector<bool> seen(membership.size(), false);
    std::unordered_set<DuplicateKey, DuplicateHash> bucket;
    int prior_chromosome = -1;
    U64 prior_start = 0;
    gzFile input = gzopen(source.c_str(), "rb");
    require(input != nullptr, "fragment input unavailable");
    gzbuffer(input, 1 << 20);
    try {
        char buffer[4096];
        while (gzgets(input, buffer, sizeof(buffer)) != nullptr) {
            std::size_t length = std::strlen(buffer);
            require(length > 0 && buffer[length - 1] == '\n', "fragment line truncated or oversized");
            std::string line(buffer, length - 1);
            if (!line.empty() && line.back() == '\r') line.pop_back();
            auto row = fields(line);
            require(row.size() == 5, "selected fragment schema differs (expected no comment header)");
            auto found_chromosome = rank.find(row[0]);
            require(found_chromosome != rank.end(), "nonprimary selected fragment chromosome");
            int chromosome = found_chromosome->second;
            U64 start = number(row[1]), end = number(row[2]), support = number(row[4]);
            require(end > start && support > 0, "invalid fragment coordinates or readSupport");
            require(chromosome > prior_chromosome ||
                    (chromosome == prior_chromosome && start >= prior_start), "fragment start order differs");
            // End order is deliberately unrestricted within each start bucket.
            if (chromosome != prior_chromosome || start != prior_start) bucket.clear();
            prior_chromosome = chromosome;
            prior_start = start;
            auto found_cell = membership.find(row[3]);
            require(found_cell != membership.end(), "fragment cell is outside native unit membership");
            const auto cell = found_cell->second;
            require(bucket.insert({end, cell.index}).second,
                    DUPLICATE_FAILURE);
            out.maximum_start_bucket_keys = std::max(out.maximum_start_bucket_keys, bucket.size());
            require(bucket.size() <= 10000000, "start bucket exceeds bounded duplicate-check memory");
            if (!seen[cell.index]) { seen[cell.index] = true; ++out.seen_cells; }
            require(++out.records <= RECORD_BOUND, "fragment record bound exceeded");
            ++out.primary[0]; ++out.primary[cell.half];
            require(support <= std::numeric_limits<U64>::max() - out.support[0], "readSupport total overflow");
            out.support[0] += support; out.support[cell.half] += support;
            bool overlapped = false;
            const bool autosomal = row[0] != "chrX" && row[0] != "chrY";
            if (autosomal) {
                ++out.autosomal[0]; ++out.autosomal[cell.half];
                auto index = axis.find(row[0]);
                if (index != axis.end()) {
                    auto hit = [&](std::size_t peak) {
                        ++out.counts[0][peak]; ++out.counts[cell.half][peak];
                        ++out.summed[0]; ++out.summed[cell.half];
                        overlapped = true;
                    };
                    index->second.overlaps(start, end, hit);
                }
            }
            if (overlapped) { ++out.union_overlap[0]; ++out.union_overlap[cell.half]; }
            if (out.records % 10000000 == 0)
                std::cerr << tag << " records=" << out.records << '\n';
        }
        int status = Z_OK;
        const char* text = gzerror(input, &status);
        require(status == Z_OK || status == Z_STREAM_END, "gzip CRC/read failure: " + std::string(text));
        require(gzeof(input) != 0, "gzip stream did not reach EOF");
        int close_status = gzclose(input);
        input = nullptr;
        require(close_status == Z_OK, "gzip close/CRC failure");
    } catch (...) { if (input != nullptr) gzclose(input); throw; }
    require(out.seen_cells == membership.size(), "selected fragment source misses native nuclei");
    for (std::size_t i = 0; i < n_peaks; ++i) {
        require(out.counts[0][i] == out.counts[1][i] + out.counts[2][i], "per-peak full=A+B differs");
        require(out.counts[0][i] <= out.records && out.counts[0][i] < (U64(1) << 31), "int32 cast unsafe");
    }
    out.elapsed_seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - began).count();
    return out;
}
void write_row(const std::filesystem::path& path, U64 offset, const std::vector<U64>& values) {
    require(std::filesystem::is_regular_file(path) && std::filesystem::file_size(path) >= offset + values.size() * 4,
            "preallocated NPY row exceeds file size");
    std::vector<unsigned char> bytes(values.size() * 4);
    for (std::size_t i = 0; i < values.size(); ++i) {
        require(values[i] < (U64(1) << 31), "int32 cast unsafe");
        auto value = std::uint32_t(values[i]);
        for (int j = 0; j < 4; ++j) bytes[i * 4 + j] = (value >> (j * 8)) & 255;
    }
    std::fstream output(path, std::ios::in | std::ios::out | std::ios::binary);
    require(bool(output), "NPY output unavailable");
    output.seekp(offset);
    output.write(reinterpret_cast<const char*>(bytes.data()), bytes.size());
    output.flush();
    require(bool(output), "NPY row write failed");
}
void receipt(const std::filesystem::path& path, const Result& out) {
    require(!std::filesystem::exists(path), "refusing receipt overwrite");
    std::ofstream file(path);
    file << "{\n\"records\":" << out.records << ",\n\"seen_cells\":" << out.seen_cells
         << ",\n\"duplicate_keys\":0,\n\"maximum_start_bucket_keys\":" << out.maximum_start_bucket_keys
         << ",\n\"gzip_crc_and_eof_verified\":true,\n\"full_equals_halves_per_peak\":true,\n\"zlib_runtime_version\":\""
         << zlibVersion() << "\",\n\"elapsed_seconds\":" << out.elapsed_seconds;
    const std::array<std::string, 5> names{"primary_fragment_records", "autosomal_fragment_records",
                                        "union_overlapped_fragments", "readSupport_diagnostic", "summed_peak_counts"};
    const std::array<std::array<U64, 3>, 5> values{out.primary, out.autosomal, out.union_overlap, out.support, out.summed};
    for (std::size_t i = 0; i < names.size(); ++i)
        file << ",\n\"" << names[i] << "\":[" << values[i][0] << ',' << values[i][1] << ',' << values[i][2] << ']';
    file << "\n}\n";
    file.flush();
    require(bool(file), "receipt write failed");
}
void gzip_fixture(const std::filesystem::path& path, const std::string& content) {
    require(!std::filesystem::exists(path), "refusing synthetic fixture overwrite");
    gzFile output = gzopen(path.c_str(), "wb");
    require(output != nullptr, "synthetic gzip output unavailable");
    require(gzwrite(output, content.data(), content.size()) == int(content.size()), "synthetic gzip write failed");
    require(gzclose(output) == Z_OK, "synthetic gzip close failed");
}
void self_checks(const std::filesystem::path& output) {
    Geometry axis;
    axis["chr1"].peaks = {{10,20,0}, {12,14,1}, {0,100,2}, {10,12,3}};
    axis["chr1"].prepare();
    Cells membership{{"well1_a", {0,1}}, {"well1_b", {1,2}}};
    std::string rows = "chr1\t0\t10\twell1_a\t1\nchr1\t10\t15\twell1_a\t99\n"
                       "chr1\t10\t15\twell1_b\t1\nchr1\t10\t12\twell1_b\t7\n"
                       "chr1\t20\t25\twell1_b\t1\n";
    auto good = output / "boundary_nested_distinct_cell.tsv.gz";
    gzip_fixture(good, rows);
    auto result = count(good, axis, membership, 4, "synthetic");
    require(result.counts[0] == std::vector<U64>({3,2,5,3}) &&
            result.counts[1] == std::vector<U64>({1,1,2,1}) &&
            result.counts[2] == std::vector<U64>({2,1,3,2}), "hand-derived overlap expectations differ");
    require(result.primary == std::array<U64,3>{5,2,3} && result.union_overlap == result.primary &&
            result.support[0] == 109 && result.summed[0] == 13, "hand-derived denominator expectations differ");
    auto bad = output / "duplicate_nonadjacent_end_order.tsv.gz";
    gzip_fixture(bad, "chr1\t10\t15\twell1_a\t1\nchr1\t10\t12\twell1_b\t1\nchr1\t10\t15\twell1_a\t2\n");
    bool refused = false;
    try { count(bad, axis, membership, 4, "synthetic-duplicate"); }
    catch (const std::runtime_error& error) { refused = std::string(error.what()).find("duplicate fragment key") != std::string::npos; }
    require(refused, "duplicate-key refusal failed");
    auto truncated = output / "truncated_crc.tsv.gz";
    std::ifstream source(good, std::ios::binary);
    std::string compressed((std::istreambuf_iterator<char>(source)), std::istreambuf_iterator<char>());
    require(compressed.size() > 8 && !std::filesystem::exists(truncated), "synthetic CRC fixture invalid");
    std::ofstream corrupt(truncated, std::ios::binary);
    corrupt.write(compressed.data(), compressed.size() - 8);
    corrupt.close();
    refused = false;
    try { count(truncated, axis, membership, 4, "synthetic-CRC"); }
    catch (const std::runtime_error& error) { refused = std::string(error.what()).find("gzip") != std::string::npos; }
    require(refused, "truncated-gzip refusal failed");
    std::ofstream evidence(output / "self_checks.json");
    evidence << "{\"boundary_touch_excluded\":true,\"nested_peaks_all_counted\":true,"
                "\"distinct_cells_identical_coordinates_counted_separately\":true,"
                "\"end_order_not_required\":true,\"readSupport_not_weighted\":true,"
                "\"duplicate_key_refused\":true,\"truncated_gzip_refused\":true}\n";
    evidence.flush();
    require(bool(evidence), "self-check receipt write failed");
}
int main(int argc, char** argv) {
    try {
        require(std::string(zlibVersion()) == ZLIB_VERSION, "zlib compile/runtime versions differ");
        if (argc == 3 && std::string(argv[1]) == "--self-check") {
            self_checks(argv[2]); return 0;
        }
        require(argc == 10, "usage: counter source geometry cells full.npy.partial A.npy.partial B.npy.partial offset expected_records receipt.json");
        auto axis = geometry(argv[2]);
        auto membership = cells(argv[3]);
        auto result = count(argv[1], axis, membership, PEAKS, argv[3]);
        require(result.records == number(argv[8]), "streamed fragment records differ from frozen manifest");
        U64 offset = number(argv[7]);
        for (std::size_t i = 0; i < 3; ++i) write_row(argv[4+i], offset, result.counts[i]);
        receipt(argv[9], result);
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "counter failure: " << error.what() << '\n';
        if (argc == 10) {
            const bool duplicate = std::string(error.what()) == DUPLICATE_FAILURE;
            const std::string path = std::string(argv[9]) + ".failure.json";
            const std::string evidence = std::string("{\"completed\":false,\"categorical_reason\":\"") +
                (duplicate ? "duplicate_fragment_key" : "counter_error") +
                "\",\"duplicate_key_refused\":" + (duplicate ? "true" : "false") + "}\n";
            // No raw fields or exception text enter this bounded receipt. Never overwrite.
            const int descriptor = ::open(path.c_str(), O_WRONLY | O_CREAT | O_EXCL, 0640);
            if (descriptor < 0) {
                std::cerr << "exclusive counter failure receipt could not be created\n";
            } else {
                const auto written = ::write(descriptor, evidence.data(), evidence.size());
                if (written != static_cast<ssize_t>(evidence.size()) || ::fsync(descriptor) != 0)
                    std::cerr << "counter failure receipt write failed\n";
                if (::close(descriptor) != 0) std::cerr << "counter failure receipt close failed\n";
            }
        }
        return 1;
    }
}
