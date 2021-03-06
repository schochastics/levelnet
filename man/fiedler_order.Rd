% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lazarus_count.R
\name{fiedler_order}
\alias{fiedler_order}
\title{Permutation induced by Fiedler vector}
\usage{
fiedler_order(g, mode = "cols")
}
\arguments{
\item{g}{an igraph object or a (0,1)-matrix}

\item{mode}{one of "mcl" (clique vertex matrix), "cols" (Lazarus count of columns) "rows" (Lazarus count of rows) or "sym" (Lazarus count of both columns and rows).}
}
\value{
numeric vector
}
\description{
Returns the permutation induced by sorting the Fiedler vector of the Laplacian matrix of a graph
}
\author{
David Schoch
}
