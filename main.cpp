#include "DSP.h"

#include <cstring>
#include <fstream>
#include <iostream>
#include <string>

std::vector<Point<double>> getData(const std::string& filename) {
  std::vector<Point<double>> ret{};
  std::fstream fs(filename, std::ios::in);
  std::string line;
  while (std::getline(fs, line)) {
    Point<double> p{};
    std::size_t start = 0;
    std::size_t pos = line.find(',', start);
    p.push_back(std::stod(line.substr(start, pos - start)));
    while (pos != std::string::npos) {
      start = pos + 1;
      pos = line.find(',', start);
      p.push_back(std::stod(line.substr(start, pos - start)));
    }
    ret.push_back(std::move(p));
  }
  return ret;
}

void writePartition(const std::vector<std::pair<std::vector<std::pair<double, double>>, double>>& partitions, const std::string& filename) {
  std::fstream fs(filename, std::ios::out);
  for (std::size_t i = 0; i < partitions.size(); ++i) {
    for (std::size_t j = 0; j < partitions[i].first.size(); ++j) {
      fs << partitions[i].first[j].first << "," << partitions[i].first[j].second << ",";
    }
    fs << partitions[i].second << "\n";
  }
  fs.close();
}

int main(int argc, char** argv) {
  if (argc != 8) {
    std::cout << "please provide\n"
                 "  1) input file (string): csv file;\n"
                 "  2) output file (string): csv file;\n"
                 "  3) splits (int): positive integer;\n"
                 "  4) maxDimRatio (double): positive double to control the shape of the partition (maxDimSize / minDimSize < maxDimRatio);\n"
                 "  5) theta (double): refer to paper;\n"
                 "  6) subsample (int): downsampling to compute discrepancy;\n"
                 "  7) minimum sample (int): only split subrectangle with #points larger than minimum sample;" << std::endl;
    exit(1);
  }
  std::cout << "reading " << argv[1] << std::endl;

  const std::vector<Point<double>> points = getData(argv[1]);
  const std::vector<std::pair<double, double>> boundaries{{0, 1}, {0, 1}};
  const std::size_t splits = std::stoul(argv[3]);
  const double maxDimRatio = std::stod(argv[4]);
  double theta = std::stod(argv[5]);
  const int subsample = std::stoi(argv[6]);
  const std::size_t minsample = std::stoul(argv[7]);
  const auto densities = DSP(points, boundaries, splits, maxDimRatio, theta, subsample, minsample);

  std::cout << "writing " << argv[2] << std::endl;
  writePartition(densities, argv[2]);
  return 0;
}

