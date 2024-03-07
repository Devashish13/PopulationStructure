\name{PCA_plot}
\alias{PCA_plot}
\title{Plotting Principal Component Analysis (PCA) Results}
\usage{
PCA_plot(pcadata_path,metadata_path,sampleid_var,color_var,x = 1,y = 2,
        shape_var = NULL,color_vec = NULL,shape_vec= NULL,figfile = "PCA_plot",...)
}

\description{
The \code{PCA_plot} function takes individual PC loadings data, metadata, and produces a scatter plot of the specified principal components. It supports dynamic aesthetic mappings for color and shape based on metadata variables, allowing for detailed visual analysis of PCA results. It also produces a population centroid PC plots for better understanding of populations distribution along with individuals PC plot.It provides publication ready plot representing genetic structure uncovered by PCA.
}


\arguments{
  \item{pcadata_path}{Path to the PCA eigenvalue and eigenvector files (without extension).}
  \item{metadata_path}{Path to the metadata file associated with the genotype data(tab separated file)}
  \item{sampleid_var}{The variable in the PCA data that matches the sample ID in the metadata.}
  \item{color_var}{Metadata variable to map point colors to}

  \item{x}{The principal component to plot on the x-axis. Default is 1.}
  \item{y}{The principal component to plot on the y-axis. Default is 2.}
  \item{shape_var}{Metadata variable to map point shapes to}
  \item{color_vec}{Optional named vector of colors to use for \code{color_var}. If NULL, distinct colors are generated automatically.}
  \item{shape_vec}{Optional named vector of shapes to use for \code{shape_var}. If NULL, shapes are assigned automatically.}
  \item{figfile}{Base filename for saving the generated plot. Default is \code{"PCA_plot"}.}
  \item{\dots}{Additional arguments passed to \code{geom_point} and other ggplot2 functions for further customization.}
}

\examples{
library(PopulationStructure)

library(gridExtra)
pcadata_path <- "geno_data_QC_PCA_RESULT"
metadata_path <- "hgdp_chr21/hgdp_chr21.psam"

plots_set <- PCA_plot(pcadata_path = pcadata_path,metadata_path = metadata_path,
         sampleid_var = "#IID",color_var = "population",shape_var = "region",size = 4)
do.call("grid.arrange", c(plots_set, ncol=2))
}

