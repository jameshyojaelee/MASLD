#include <zlib.h>

#include <algorithm>
#include <array>
#include <cctype>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <string_view>
#include <unordered_map>
#include <utility>
#include <vector>

namespace {

constexpr std::size_t kOligoLength = 107;
constexpr std::size_t kBarcodeLength = 20;
constexpr std::string_view kRead1Adapter =
    "GGCCTAACTGGCCGGTACCCTGAGTACTGTATGGGCGGGTACC";
constexpr std::string_view kRead2Adapter = "TCACCATGGTGGCTTTACCAACAG";

class ReconstructionError : public std::runtime_error {
 public:
  using std::runtime_error::runtime_error;
};

class GzipLineReader {
 public:
  explicit GzipLineReader(const std::string& path)
      : path_(path), file_(gzopen(path.c_str(), "rb")), buffer_(1U << 20) {
    if (file_ == nullptr) {
      throw ReconstructionError("cannot open gzip input: " + path);
    }
  }

  ~GzipLineReader() {
    if (file_ != nullptr) {
      gzclose(file_);
    }
  }

  GzipLineReader(const GzipLineReader&) = delete;
  GzipLineReader& operator=(const GzipLineReader&) = delete;

  bool GetLine(std::string& output) {
    output.clear();
    while (true) {
      if (position_ == available_) {
        const int observed =
            gzread(file_, buffer_.data(), static_cast<unsigned int>(buffer_.size()));
        if (observed < 0) {
          int error_number = 0;
          const char* message = gzerror(file_, &error_number);
          throw ReconstructionError("gzip read failed for " + path_ + ": " +
                                    (message == nullptr ? "unknown" : message));
        }
        position_ = 0;
        available_ = static_cast<std::size_t>(observed);
        if (available_ == 0) {
          return !output.empty();
        }
      }
      const auto begin = buffer_.begin() + static_cast<std::ptrdiff_t>(position_);
      const auto end = buffer_.begin() + static_cast<std::ptrdiff_t>(available_);
      const auto newline = std::find(begin, end, '\n');
      output.append(begin, newline);
      position_ += static_cast<std::size_t>(newline - begin);
      if (newline != end) {
        ++position_;
        if (!output.empty() && output.back() == '\r') {
          output.pop_back();
        }
        return true;
      }
    }
  }

 private:
  std::string path_;
  gzFile file_ = nullptr;
  std::vector<char> buffer_;
  std::size_t position_ = 0;
  std::size_t available_ = 0;
};

struct FastqRecord {
  std::string name;
  std::string sequence;
  std::string quality;
};

bool ReadFastq(GzipLineReader& reader, FastqRecord& record) {
  std::string plus;
  if (!reader.GetLine(record.name)) {
    return false;
  }
  if (!reader.GetLine(record.sequence) || !reader.GetLine(plus) ||
      !reader.GetLine(record.quality)) {
    throw ReconstructionError("truncated FASTQ record");
  }
  if (record.name.empty() || record.name.front() != '@' || plus.empty() ||
      plus.front() != '+' || record.sequence.size() != record.quality.size()) {
    throw ReconstructionError("FASTQ structure differs");
  }
  return true;
}

std::string_view ReadNameKey(const std::string& header) {
  std::size_t end = header.find_first_of(" \t");
  if (end == std::string::npos) {
    end = header.size();
  }
  std::string_view key(header.data() + 1, end - 1);
  if (key.size() > 2 && key[key.size() - 2] == '/' &&
      (key.back() == '1' || key.back() == '2')) {
    key.remove_suffix(2);
  }
  return key;
}

int BaseIndex(char base) {
  switch (std::toupper(static_cast<unsigned char>(base))) {
    case 'A':
      return 0;
    case 'C':
      return 1;
    case 'G':
      return 2;
    case 'T':
      return 3;
    default:
      return -1;
  }
}

bool EncodeBarcode(std::string_view barcode, std::uint64_t& encoded) {
  if (barcode.size() != kBarcodeLength) {
    return false;
  }
  encoded = 0;
  for (char base : barcode) {
    const int value = BaseIndex(base);
    if (value < 0) {
      return false;
    }
    encoded = (encoded << 2U) | static_cast<std::uint64_t>(value);
  }
  return true;
}

bool StartsWith(const std::string& sequence, std::string_view prefix) {
  return sequence.size() >= prefix.size() &&
         std::equal(prefix.begin(), prefix.end(), sequence.begin());
}

std::vector<std::string_view> SplitFirstFields(const std::string& line,
                                               std::size_t count) {
  std::vector<std::string_view> output;
  output.reserve(count);
  std::size_t start = 0;
  while (output.size() < count) {
    const std::size_t end = line.find('\t', start);
    if (end == std::string::npos) {
      if (output.size() + 1 != count) {
        throw ReconstructionError("count row has too few columns");
      }
      output.emplace_back(line.data() + start, line.size() - start);
      break;
    }
    output.emplace_back(line.data() + start, end - start);
    start = end + 1;
  }
  return output;
}

struct Construct {
  std::string name;
  std::uint64_t barcode_count = 0;
  std::uint64_t mapped_reads = 0;
  std::array<std::array<std::uint32_t, 4>, kOligoLength> base_counts{};
  std::string consensus;
  std::uint64_t exact_consensus_reads = 0;
  std::uint64_t minimum_position_coverage = 0;
  double minimum_consensus_fraction = 0.0;
  double exact_consensus_fraction = 0.0;
  std::string status;
};

struct PassStats {
  std::uint64_t paired_reads = 0;
  std::uint64_t exact_adapters = 0;
  std::uint64_t invalid_barcode = 0;
  std::uint64_t unmapped_barcode = 0;
  std::uint64_t mapped_barcode = 0;
  std::uint64_t oligos_with_non_acgt = 0;
};

struct Arguments {
  std::string mapping_counts;
  std::string fastq_r1;
  std::string fastq_r2;
  std::string output_tsv;
  std::string output_stats;
  double minimum_consensus_fraction = 0.99;
  double minimum_exact_fraction = 0.50;
  std::uint64_t minimum_position_coverage = 20;
};

Arguments ParseArguments(int argc, char** argv) {
  Arguments arguments;
  for (int index = 1; index < argc; index += 2) {
    if (index + 1 >= argc) {
      throw ReconstructionError("command-line option lacks a value");
    }
    const std::string option = argv[index];
    const std::string value = argv[index + 1];
    if (option == "--mapping-counts") {
      arguments.mapping_counts = value;
    } else if (option == "--fastq-r1") {
      arguments.fastq_r1 = value;
    } else if (option == "--fastq-r2") {
      arguments.fastq_r2 = value;
    } else if (option == "--output-tsv") {
      arguments.output_tsv = value;
    } else if (option == "--output-stats") {
      arguments.output_stats = value;
    } else if (option == "--minimum-consensus-fraction") {
      arguments.minimum_consensus_fraction = std::stod(value);
    } else if (option == "--minimum-exact-fraction") {
      arguments.minimum_exact_fraction = std::stod(value);
    } else if (option == "--minimum-position-coverage") {
      arguments.minimum_position_coverage = std::stoull(value);
    } else {
      throw ReconstructionError("unknown command-line option: " + option);
    }
  }
  if (arguments.mapping_counts.empty() || arguments.fastq_r1.empty() ||
      arguments.fastq_r2.empty() || arguments.output_tsv.empty() ||
      arguments.output_stats.empty() ||
      arguments.minimum_consensus_fraction <= 0.0 ||
      arguments.minimum_consensus_fraction > 1.0 ||
      arguments.minimum_exact_fraction <= 0.0 ||
      arguments.minimum_exact_fraction > 1.0 ||
      arguments.minimum_position_coverage == 0) {
    throw ReconstructionError("command-line contract differs");
  }
  return arguments;
}

void LoadBarcodeMap(const std::string& path,
                    std::vector<Construct>& constructs,
                    std::unordered_map<std::uint64_t, std::uint32_t>& barcode_map,
                    std::uint64_t& mapping_rows,
                    std::uint64_t& invalid_mapping_barcodes) {
  GzipLineReader reader(path);
  std::unordered_map<std::string, std::uint32_t> construct_ids;
  construct_ids.reserve(12000);
  barcode_map.reserve(25000000);
  barcode_map.max_load_factor(0.80F);
  std::string line;
  while (reader.GetLine(line)) {
    const auto fields = SplitFirstFields(line, 3);
    const std::string construct_name(fields[0]);
    const std::string_view barcode = fields[1];
    if (construct_name.empty()) {
      throw ReconstructionError("empty construct identifier");
    }
    auto construct_match = construct_ids.find(construct_name);
    std::uint32_t construct_id = 0;
    if (construct_match == construct_ids.end()) {
      if (constructs.size() >= std::numeric_limits<std::uint32_t>::max()) {
        throw ReconstructionError("construct count overflows identifier type");
      }
      construct_id = static_cast<std::uint32_t>(constructs.size());
      construct_ids.emplace(construct_name, construct_id);
      Construct construct;
      construct.name = construct_name;
      constructs.push_back(std::move(construct));
    } else {
      construct_id = construct_match->second;
    }
    ++mapping_rows;
    std::uint64_t encoded = 0;
    if (!EncodeBarcode(barcode, encoded)) {
      ++invalid_mapping_barcodes;
      continue;
    }
    const auto insertion = barcode_map.emplace(encoded, construct_id);
    if (!insertion.second) {
      throw ReconstructionError("barcode is duplicated in the mapping count file");
    }
    ++constructs[construct_id].barcode_count;
  }
  std::size_t alternative_constructs = 0;
  for (const Construct& construct : constructs) {
    if (construct.name.size() >= 4 &&
        construct.name.compare(construct.name.size() - 4, 4, "_Mut") == 0) {
      ++alternative_constructs;
    }
  }
  if (constructs.size() != 10797 || alternative_constructs != 5361 ||
      constructs.size() - alternative_constructs != 5436 ||
      mapping_rows != barcode_map.size() + invalid_mapping_barcodes) {
    throw ReconstructionError("current mapping count construct census differs");
  }
}

template <typename Callback>
PassStats ScanPairedFastq(
    const std::string& fastq_r1, const std::string& fastq_r2,
    const std::unordered_map<std::uint64_t, std::uint32_t>& barcode_map,
    Callback callback) {
  GzipLineReader r1_reader(fastq_r1);
  GzipLineReader r2_reader(fastq_r2);
  FastqRecord r1;
  FastqRecord r2;
  PassStats stats;
  while (true) {
    const bool has_r1 = ReadFastq(r1_reader, r1);
    const bool has_r2 = ReadFastq(r2_reader, r2);
    if (has_r1 != has_r2) {
      throw ReconstructionError("paired FASTQ record census differs");
    }
    if (!has_r1) {
      break;
    }
    ++stats.paired_reads;
    if (ReadNameKey(r1.name) != ReadNameKey(r2.name)) {
      throw ReconstructionError("paired FASTQ identifiers differ");
    }
    if (!StartsWith(r1.sequence, kRead1Adapter) ||
        !StartsWith(r2.sequence, kRead2Adapter) ||
        r1.sequence.size() < kRead1Adapter.size() + kOligoLength ||
        r2.sequence.size() < kRead2Adapter.size() + kBarcodeLength) {
      continue;
    }
    ++stats.exact_adapters;
    const std::string_view barcode(
        r2.sequence.data() + static_cast<std::ptrdiff_t>(kRead2Adapter.size()),
        kBarcodeLength);
    std::uint64_t encoded = 0;
    if (!EncodeBarcode(barcode, encoded)) {
      ++stats.invalid_barcode;
      continue;
    }
    const auto match = barcode_map.find(encoded);
    if (match == barcode_map.end()) {
      ++stats.unmapped_barcode;
      continue;
    }
    ++stats.mapped_barcode;
    std::string oligo = r1.sequence.substr(kRead1Adapter.size(), kOligoLength);
    bool all_acgt = true;
    for (char& base : oligo) {
      base = static_cast<char>(std::toupper(static_cast<unsigned char>(base)));
      if (BaseIndex(base) < 0) {
        all_acgt = false;
      }
    }
    if (!all_acgt) {
      ++stats.oligos_with_non_acgt;
    }
    callback(match->second, oligo);
  }
  return stats;
}

void BuildConsensus(std::vector<Construct>& constructs) {
  for (Construct& construct : constructs) {
    construct.consensus.reserve(kOligoLength);
    construct.minimum_position_coverage =
        std::numeric_limits<std::uint64_t>::max();
    construct.minimum_consensus_fraction = 1.0;
    for (const auto& counts : construct.base_counts) {
      const std::uint64_t coverage = static_cast<std::uint64_t>(counts[0]) +
                                     counts[1] + counts[2] + counts[3];
      const auto maximum = std::max_element(counts.begin(), counts.end());
      const bool tied =
          std::count(counts.begin(), counts.end(), *maximum) != 1;
      if (coverage == 0 || tied) {
        construct.consensus.push_back('N');
        construct.minimum_position_coverage = 0;
        construct.minimum_consensus_fraction = 0.0;
        continue;
      }
      const std::size_t index =
          static_cast<std::size_t>(maximum - counts.begin());
      construct.consensus.push_back(std::string_view("ACGT")[index]);
      construct.minimum_position_coverage =
          std::min(construct.minimum_position_coverage, coverage);
      construct.minimum_consensus_fraction = std::min(
          construct.minimum_consensus_fraction,
          static_cast<double>(*maximum) / static_cast<double>(coverage));
    }
    if (construct.minimum_position_coverage ==
        std::numeric_limits<std::uint64_t>::max()) {
      construct.minimum_position_coverage = 0;
    }
  }
}

void FinalizeStatuses(std::vector<Construct>& constructs,
                      const Arguments& arguments) {
  for (Construct& construct : constructs) {
    construct.exact_consensus_fraction =
        construct.mapped_reads == 0
            ? 0.0
            : static_cast<double>(construct.exact_consensus_reads) /
                  static_cast<double>(construct.mapped_reads);
    if (construct.minimum_position_coverage <
        arguments.minimum_position_coverage) {
      construct.status = "below_minimum_position_coverage";
    } else if (construct.minimum_consensus_fraction <
               arguments.minimum_consensus_fraction) {
      construct.status = "below_minimum_consensus_fraction";
    } else if (construct.exact_consensus_fraction <
               arguments.minimum_exact_fraction) {
      construct.status = "below_minimum_exact_consensus_fraction";
    } else if (construct.consensus.size() != kOligoLength ||
               construct.consensus.find('N') != std::string::npos) {
      construct.status = "consensus_unresolved";
    } else {
      construct.status = "pass";
    }
  }
}

void WriteOutputs(const std::vector<Construct>& constructs,
                  std::uint64_t mapping_rows,
                  std::uint64_t invalid_mapping_barcodes,
                  const PassStats& first,
                  const PassStats& second, const Arguments& arguments) {
  std::ofstream table(arguments.output_tsv);
  if (!table) {
    throw ReconstructionError("cannot create construct output table");
  }
  table << "construct\tallele\tbarcode_count\tmapped_reads"
        << "\tminimum_position_coverage\tminimum_consensus_fraction"
        << "\texact_consensus_reads\texact_consensus_fraction"
        << "\tconsensus_107bp\tstatus\n";
  table << std::setprecision(12);
  std::uint64_t passing_constructs = 0;
  for (const Construct& construct : constructs) {
    const bool alternative =
        construct.name.size() >= 4 &&
        construct.name.compare(construct.name.size() - 4, 4, "_Mut") == 0;
    table << construct.name << '\t' << (alternative ? "alt" : "ref") << '\t'
          << construct.barcode_count << '\t' << construct.mapped_reads << '\t'
          << construct.minimum_position_coverage << '\t'
          << construct.minimum_consensus_fraction << '\t'
          << construct.exact_consensus_reads << '\t'
          << construct.exact_consensus_fraction << '\t' << construct.consensus
          << '\t' << construct.status << '\n';
    if (construct.status == "pass") {
      ++passing_constructs;
    }
  }
  table.close();
  if (!table) {
    throw ReconstructionError("construct output table write failed");
  }
  std::ofstream stats(arguments.output_stats);
  if (!stats) {
    throw ReconstructionError("cannot create reconstruction statistics");
  }
  stats << "{\n"
        << "  \"schema_version\": \"masld-bench-gse281364-fastq-reconstruction-v1\",\n"
        << "  \"status\": \"pass\",\n"
        << "  \"mapping_rows\": " << mapping_rows << ",\n"
        << "  \"invalid_mapping_barcodes\": "
        << invalid_mapping_barcodes << ",\n"
        << "  \"constructs\": " << constructs.size() << ",\n"
        << "  \"passing_constructs\": " << passing_constructs << ",\n"
        << "  \"first_pass_paired_reads\": " << first.paired_reads << ",\n"
        << "  \"first_pass_exact_adapters\": " << first.exact_adapters << ",\n"
        << "  \"first_pass_mapped_barcodes\": " << first.mapped_barcode << ",\n"
        << "  \"first_pass_unmapped_barcodes\": " << first.unmapped_barcode << ",\n"
        << "  \"first_pass_invalid_barcodes\": " << first.invalid_barcode << ",\n"
        << "  \"first_pass_oligos_with_non_acgt\": "
        << first.oligos_with_non_acgt << ",\n"
        << "  \"second_pass_paired_reads\": " << second.paired_reads << ",\n"
        << "  \"second_pass_exact_adapters\": " << second.exact_adapters << ",\n"
        << "  \"second_pass_mapped_barcodes\": " << second.mapped_barcode << ",\n"
        << "  \"minimum_position_coverage\": "
        << arguments.minimum_position_coverage << ",\n"
        << "  \"minimum_consensus_fraction\": "
        << arguments.minimum_consensus_fraction << ",\n"
        << "  \"minimum_exact_consensus_fraction\": "
        << arguments.minimum_exact_fraction << ",\n"
        << "  \"observed_count_columns_loaded\": false,\n"
        << "  \"sealed_outcomes_loaded\": false\n"
        << "}\n";
  stats.close();
  if (!stats) {
    throw ReconstructionError("reconstruction statistics write failed");
  }
}

}  // namespace

int main(int argc, char** argv) {
  try {
    const Arguments arguments = ParseArguments(argc, argv);
    std::vector<Construct> constructs;
    constructs.reserve(11000);
    std::unordered_map<std::uint64_t, std::uint32_t> barcode_map;
    std::uint64_t mapping_rows = 0;
    std::uint64_t invalid_mapping_barcodes = 0;
    LoadBarcodeMap(arguments.mapping_counts, constructs, barcode_map,
                   mapping_rows, invalid_mapping_barcodes);
    const PassStats first = ScanPairedFastq(
        arguments.fastq_r1, arguments.fastq_r2, barcode_map,
        [&constructs](std::uint32_t construct_id, const std::string& oligo) {
          Construct& construct = constructs.at(construct_id);
          ++construct.mapped_reads;
          for (std::size_t position = 0; position < kOligoLength; ++position) {
            const int base = BaseIndex(oligo[position]);
            if (base >= 0) {
              auto& count = construct.base_counts[position][base];
              if (count == std::numeric_limits<std::uint32_t>::max()) {
                throw ReconstructionError("per-base count overflow");
              }
              ++count;
            }
          }
        });
    if (first.paired_reads != 71771807 || first.mapped_barcode == 0) {
      throw ReconstructionError("raw paired-read or mapped-barcode census differs");
    }
    BuildConsensus(constructs);
    const PassStats second = ScanPairedFastq(
        arguments.fastq_r1, arguments.fastq_r2, barcode_map,
        [&constructs](std::uint32_t construct_id, const std::string& oligo) {
          Construct& construct = constructs.at(construct_id);
          if (oligo == construct.consensus) {
            ++construct.exact_consensus_reads;
          }
        });
    if (first.paired_reads != second.paired_reads ||
        first.exact_adapters != second.exact_adapters ||
        first.mapped_barcode != second.mapped_barcode ||
        first.unmapped_barcode != second.unmapped_barcode ||
        first.invalid_barcode != second.invalid_barcode) {
      throw ReconstructionError("paired FASTQ passes differ");
    }
    FinalizeStatuses(constructs, arguments);
    WriteOutputs(constructs, mapping_rows, invalid_mapping_barcodes,
                 first, second, arguments);
  } catch (const std::exception& error) {
    std::cerr << "GSE281364 reconstruction failed: " << error.what() << '\n';
    return 1;
  }
  return 0;
}
