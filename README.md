<p align="center">
  <img src="inst/figures/logo.png" alt="SNPtools logo" width="200"/>
</p>

# SNPtools

`SNPtools` is an R package designed for manipulation, organization, and analysis of genotypic data, with a strong focus on integration with tools such as **FImpute** and **PLINK**.

It provides robust S4-based data structures for storing genotypes and marker maps, along with functions to combine different genotype panels, summarize data, and prepare files for imputation and selection pipelines.

---

## 📦 Installation

To install directly from GitHub, use the `devtools` package:

```r
# Install devtools if not already installed
install.packages("devtools")

# Install SNPtools from GitHub
devtools::install_github("viniciusjunqueira/SNPtools")
```
---

## Usage

```r
library(SNPtools)

# Example: print a formatted header
qc_header("Starting genotype preprocessing")

# Example configuration list
configs <- list(
  list(path = "panel1", fields = list(sample = 2, snp = 1, allele1 = 7, allele2 = 8, confidence = 9)),
  list(path = "panel2", fields = list(sample = 2, snp = 1, allele1 = 7, allele2 = 8, confidence = 9), threshold = 0.10)
)

# Import and combine genotype data
combined_data <- import_geno_list(configs)

# Print a summary
summary(combined_data)
```