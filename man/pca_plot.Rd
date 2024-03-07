\name{PCA_plot}
\alias{PCA_plot}
\title{Plotting Principal Component Analysis (PCA) Results}
\usage{
PCA_plot(pcadata_path, metadata_path, x = 1, y = 2,
         sampleid_var = "#IID", color_var = "population",
         shape_var = "region", color_vec = NULL, shape_vec = NULL,
         figfile = "PCA_plot", ...)
}

\description{
Generates a plot of PCA results, highlighting the distribution of samples based on principal components. Allows customization of aesthetic elements such as color and shape by population or region.
}

\details{
The \code{PCA_plot} function takes PCA results and metadata, and produces a customizable scatter plot of the specified principal components. It supports dynamic aesthetic mappings for color and shape based on metadata variables, allowing for detailed visual analysis of PCA results.
}

\arguments{
  \item{pcadata_path}{Path to the PCA eigenvalue and eigenvector files (without extension).}
  \item{metadata_path}{Path to the metadata CSV file associated with the PCA data.}
  \item{x}{The principal component to plot on the x-axis. Default is 1.}
  \item{y}{The principal component to plot on the y-axis. Default is 2.}
  \item{sampleid_var}{The variable in the PCA data that matches the sample ID in the metadata. Default is \code{"#IID"}.}
  \item{color_var}{Metadata variable to map point colors to. Default is \code{"population"}.}
  \item{shape_var}{Metadata variable to map point shapes to. Default is \code{"region"}.}
  \item{color_vec}{Optional named vector of colors to use for \code{color_var}. If NULL, distinct colors are generated automatically.}
  \item{shape_vec}{Optional named vector of shapes to use for \code{shape_var}. If NULL, shapes are assigned automatically.}
  \item{figfile}{Base filename for saving the generated plot. Default is \code{"PCA_plot"}.}
  \item{\dots}{Additional arguments passed to \code{geom_point} and other ggplot2 functions for further customization.}
}

\examples{
PCA_plot("path/to/pcadata", "path/to/metadata", x = 1, y = 2,
         sampleid_var = "#IID", color_var = "population",
         shape_var = "region", figfile = "my_PCA_plot")
}

\description{
This function creates a scatter plot of PCA results, using metadata to enhance the visualization with colors and shapes representing different groups. It provides an insightful way to explore and present the genetic structure uncovered by PCA.
}
