\name{PCA_GD}
\alias{PCA_GD}
\title{Principal Component Analysis (PCA) on Genotype Data}
\usage{
PCA_GD(data_path, metadata_path, geno_format = "pfile",
output_filename = "geno_data", pruned_data = "no", maf = 0.01, geno = 0.02,
indep_pairwise = c(50, 5, 0.2), hwe = 1e-8, max_alleles = 2, ...)
}

\description{
Performs principal component analysis (PCA) on genotype data, allowing for various data filtering and manipulation parameters to fine-tune the analysis process.
}

\details{
The \code{PCA_GD} function is designed to handle genotype data for principal component analysis. It supports a range of parameters for data import, filtering based on minor allele frequency, genotype missingness, linkage disequilibrium pruning, and more.
}

\arguments{
  \item{data_path}{Path to the genotype data file.}
  \item{metadata_path}{Path to the metadata file associated with the genotype data.}
  \item{geno_format}{Format of the genotype data. Default is \code{"pfile"}. Other formats must be specified.}
  \item{output_filename}{Base name for output files. Default is \code{"geno_data"}.}
  \item{pruned_data}{Specifies whether to prune the data. Accepts \code{"yes"} or \code{"no"}. Default is \code{"no"}.}
  \item{maf}{Minor allele frequency threshold for filtering variants. Default is \code{0.01}.}
  \item{geno}{Genotype missingness threshold for filtering variants. Default is \code{0.02}.}
  \item{indep_pairwise}{Parameters for the linkage disequilibrium pruning: window size, step size, and r^2 threshold, respectively. Default is \code{c(50, 5, 0.2)}.}
  \item{hwe}{Threshold for Hardy-Weinberg equilibrium test p-value for filtering variants. Default is \code{1e-8}.}
  \item{max_alleles}{Maximum number of alleles for a variant to be included in the analysis. Default is \code{2}.}
  \item{\dots}{Additional arguments passed to internal methods.}
}

\examples{
PCA_GD("path/to/geno_data", "path/to/metadata", geno_format = "pfile", sep = "\t")
}

\description{
This function performs PCA on genotype data, providing insights into the genetic structure of the sample population. It is highly configurable, allowing users to adjust the analysis based on their specific research needs.
}

