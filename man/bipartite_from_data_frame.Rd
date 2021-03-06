% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_graphs.R
\name{bipartite_from_data_frame}
\alias{bipartite_from_data_frame}
\title{two-mode network from a data.frame}
\usage{
bipartite_from_data_frame(df, type1, type2)
}
\arguments{
\item{df}{data.frame}

\item{type1}{column name of mode 1}

\item{type2}{column name of mode 2}
}
\value{
two mode network as igraph object
}
\description{
Create a two-mode network from a data.frame
}
\author{
David Schoch
}
