#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
void write_fimpute_cpp(CharacterMatrix geno_mat, CharacterVector ids, std::string file_path) {
  int n_ind = geno_mat.nrow();
  int n_snp = geno_mat.ncol();

  FILE *f = fopen(file_path.c_str(), "wt");
  if (!f) stop("❌ Cannot open file for writing.");

  fprintf(f, "ID    Chip                   Call...\n");

  int id_width = 0;
  for (int i = 0; i < ids.size(); i++) {
    if ((int) Rf_length(ids[i]) > id_width) {
      id_width = Rf_length(ids[i]);
    }
  }

  for (int i = 0; i < n_ind; i++) {
    std::string id = Rcpp::as<std::string>(ids[i]);

    std::string geno_line;
    for (int j = 0; j < n_snp; j++) {
      std::string g = Rcpp::as<std::string>(geno_mat(i, j));
      if (g == "A/A") geno_line += "0";
      else if (g == "A/B") geno_line += "1";
      else if (g == "B/B") geno_line += "2";
      else geno_line += "5";
    }

    fprintf(f, "%-*s 1 %s\n", id_width, id.c_str(), geno_line.c_str());

    if (i % 500 == 0) Rcpp::Rcout << "Written " << i << " individuals...\n";
  }

  fclose(f);
}
