#include <Rcpp.h>
#include <fstream>
#include <string>
#include <vector>

using namespace Rcpp;

// [[Rcpp::export]]
IntegerMatrix readFImputeCpp(std::string genotype_file, int nrows, int nsnps) {
  std::ifstream infile(genotype_file.c_str());
  if (!infile.is_open()) stop("Could not open genotype file.");

  std::string line;
  std::getline(infile, line); // skip header

  IntegerMatrix mat(nrows, nsnps);
  int row = 0;

  while (std::getline(infile, line) && row < nrows) {
    size_t first_tab = line.find('\t');
    size_t second_tab = line.find('\t', first_tab + 1);
    if (first_tab == std::string::npos || second_tab == std::string::npos) {
      stop("Unexpected format in line ", row + 2);
    }

    std::string geno_str = line.substr(second_tab + 1);
    for (int i = 0; i < nsnps; ++i) {
      int code = geno_str[i] - '0';
      if (code == 2) mat(row, i) = 3;
      else if (code == 1) mat(row, i) = 2;
      else if (code == 0) mat(row, i) = 1;
      else mat(row, i) = NA_INTEGER;
    }
    row++;
  }

  infile.close();
  return mat;
}
