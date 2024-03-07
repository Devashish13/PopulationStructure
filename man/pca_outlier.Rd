\name{PCA_outliers}
\alias{PCA_outliers}
\title{Identify Outlier Individuals in PCA}
\usage{
PCA_outliers(pcadata_path, metadata_path, color_var, sampleid_var,
             p_threshold = 1e-4, ncomp = 3, x = 1, y = 2, shape_var = NULL,
             shape_vec = NULL, ...)
}

\description{
Identifies and visualizes outlier individuals for each population based on principal component analysis (PCA) results and a user-defined significance threshold. This function facilitates the exploration of genetic variation within and between populations.
}

\details{
The \code{PCA_outliers} function analyzes PCA results, using a specified p-value threshold to determine outliers in the context of principal components. It optionally produces a plot to visually represent these outliers, highlighting their distribution across specified principal components.
}

\arguments{
  \item{pcadata_path}{Path to the PCA eigenvalue and eigenvector files (without extension).}
  \item{metadata_path}{Path to the metadata CSV file associated with the PCA data.}
  \item{color_var}{Metadata variable to map point colors to, typically representing different populations.}
  \item{sampleid_var}{The variable in the PCA data that matches the sample ID in the metadata.}
  \item{p_threshold}{P-value threshold for outlier detection. Individuals with p-values below this threshold are considered outliers. Default is \code{1e-4}.}
  \item{ncomp}{Number of principal components to consider for outlier detection. Default is 3.}
  \item{x}{The principal component to plot on the x-axis for visualizing outliers. Default is 1.}
  \item{y}{The principal component to plot on the y-axis for visualizing outliers. Default is 2.}
  \item{shape_var}{Optional metadata variable to map point shapes to. Default is NULL.}
  \item{\dots}{Additional arguments passed to plotting functions for further customization.}
}

\examples{
PCA_outliers("path/to/pcadata", "path/to/metadata", "population", "#IID",
             p_threshold = 1e-5, ncomp = 4, x = 1, y = 2)
}

\description{
This function offers a powerful tool for PCA-based genetic studies, enabling researchers to easily identify and investigate outliers within populations. It supports comprehensive analysis and visualization capabilities to enhance the interpretation of PCA results.
}
