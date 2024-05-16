\name{PCA_GD}
\alias{PCA_GD}
\title{Principal Component Analysis (PCA) on Genotype Data}
\usage{
PCA_GD(data_path, metadata_path,geno_format = "pfile",
       output_filename = "geno_data",pruned_data = "no",maf = 0.01,
       geno = 0.02,indep_pairwise = c(50,5,0.2),hwe = 1e-8,max_alleles = 2)
}

\description{
The \code{PCA_GD} function is designed to handle genotype data for principal component analysis. It supports a range of parameters for data import, filtering based on minor allele frequency, genotype missingness, linkage disequilibrium pruning, and other parameters described in GenomeAsia 100K project PCA pipeline. The values of the parameters can be configured based on user requirement.
}

\arguments{
  \item{data_path}{Path to the genotype data file.}
  \item{metadata_path}{Path to the metadata file associated with the genotype data(tab separated file)}
  \item{geno_format}{Format of the genotype data. Default is \code{"pfile"}. Other formats must be specified.}
  \item{output_filename}{Base name for output files. Default is \code{"geno_data"}.}
  \item{pruned_data}{Specifies whether genotype data is pruned. Accepts \code{"yes"} or \code{"no"}. Default is \code{"no"}.}
  \item{maf}{Minor allele frequency threshold for filtering variants. Default is \code{0.01}.}
  \item{geno}{Genotype missingness threshold for filtering variants. Default is \code{0.02}.}
  \item{indep_pairwise}{Parameters for the linkage disequilibrium pruning: window size, step size, and r^2 threshold, respectively. Default is \code{c(50, 5, 0.2)}.}
  \item{hwe}{Threshold for Hardy-Weinberg equilibrium test p-value for filtering variants. Default is \code{1e-8}.}
  \item{max_alleles}{Maximum number of alleles for a variant to be included in the analysis. Default is \code{2}.}
}

\examples{
library(PopulationStructure)
data_path <- "hgdp_chr21/hgdp_chr21"
metadata_path <- "hgdp_chr21/hgdp_chr21.psam"
PCA_GD(data_path = data_path, metadata_path = metadata_path, geno_format = "pfile")
}

