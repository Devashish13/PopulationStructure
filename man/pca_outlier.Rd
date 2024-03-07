\name{PCA_outliers}
\alias{PCA_outliers}
\title{Identify Outlier Individuals in PCA}
\usage{
PCA_outliers(pcadata_path,metadata_path,color_var,sampleid_var,
              p_threshold = 1e-4,ncomp = 3,x = 1,y=2,shape_var = NULL,...)
}

\description{
The \code{PCA_outliers} function analyzes PCA results, using a specified pvalue threshold based on Mahalanobis distance to determine outliers in the context of principal components. It  produces a PDF file containing population specific plots along with outlier individuals colored red based on the user-define pvalue threshold. The function saves the population wise plots in a pdf (\code{pop_specific_plots_D2_pvalue.pdf}) and saves the mahalanobis distance and corresponding pvalues in a text file (\code{PCA_data_outlier_pvalue.txt})
}

\arguments{
  \item{pcadata_path}{Path to the PCA eigenvalue and eigenvector files (without extension).}
  \item{metadata_path}{Path to the metadata file associated with the genotype data(tab separated file)}
  \item{color_var}{Metadata variable to map point colors to}
  \item{sampleid_var}{The variable in the PCA data that matches the sample ID in the metadata.}
  \item{p_threshold}{P-value threshold for outlier detection. Individuals with p-values below this threshold are considered outliers. Default is \code{1e-4}.}
  \item{ncomp}{Number of principal components to consider for outlier detection. Default is 3.}
  \item{x}{The principal component to plot on the x-axis for visualizing outliers. Default is 1.}
  \item{y}{The principal component to plot on the y-axis for visualizing outliers. Default is 2.}
  \item{shape_var}{Optional metadata variable to add population grouping. Default is NULL.}
  \item{\dots}{Additional arguments passed to plotting functions for further customization.}
}

\examples{
pcadata_path <- "geno_data_QC_PCA_RESULT"
metadata_path <- "hgdp_chr21/hgdp_chr21.psam"
PCA_outliers(pcadata_path,metadata_path,color_var = "population",
             sampleid_var = "#IID", p_threshold = 1e-4,ncomp = 3,x = 1,y=2)
}
