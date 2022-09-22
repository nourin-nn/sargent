# Sargent package documentation and import directives

#' The sargent package
#' 
#' \code{sargent} is a member of the Sanofi - Precision Medicine & Computation
#' Biology (PMCB) global research function - Single Cell Biology (SCB) framework
#' of tools and serves as the designated package for annotating scRNA-seq data.
#' For additional details regarding the use of the \code{sargent} package see
#' the vignettes:\cr
#' \code{browseVignettes("sargent")}
#' 
#' @section  Celltype annotating:
#' \itemize{
#'   \item  \link{sargentAnnotation}:     Infer the cell-type of origin of
#'                                        each single cell.
#'   \item  \link{fetchAssignment}:       Retrieves cell-type assignments in
#'                                        a \link{sargentObject}.
#'   \item  \link{fetchCelltype}:         Retrieves cell-type(s) for given
#'                                        cell(s) in a \link{sargentObject}.
#'   \item  \link{fetchCell}:             Retrieves cell-id(s) for given
#'                                        cell-type(s) in a \link{sargentObject}.
#'   \item  \link{fetchAssignmentScore}:  Retrieves calculated score(s) for
#'                                        given cell(s) or cell-type(s) in a 
#'                                        \link{sargentObject}.
#'   \item  \link{fetchGiniIndex}:        Retrieves calculated gini coefficients 
#'                                        for given cell(s) or cell-type(s) in a 
#'                                        \link{sargentObject}.
#'   \item  \link{fetchDensityPlot}:      Generates a density-plot for 
#'                                        calculated gini coefficients in a 
#'                                        \link{sargentObject}.
#'   \item  \link{fetchDotPlot}:          Generates a gene-expression dot-plot 
#'                                        for given gene(s) and cell-type(s) in 
#'                                        a \link{sargentObject}.
#' }
#' 
#' @name     sargent
#' @docType  package
#' @references
#' \enumerate{
#'   \item  Nima Nouri, et al.: A marker gene based method for identifying the 
#'            cell-type of origin from single-cell RNA sequencing data. 
#'            <journal>. <year> <issue>(<volume>):<page-page>.
#' }
#'
#' @import      ggplot2
#' @import      methods
#' @importFrom  ggdist      stat_slab stat_pointinterval
#' @importFrom  dplyr       n %>% group_by ungroup
#'                          filter select mutate summarize rename
#' @importFrom  plyr        ldply
#' @importFrom  igraph      graph_from_adjacency_matrix as_adjacency_matrix
#'                          igraph_opt as_adj_list
#' @importFrom  Matrix      sparseMatrix rowSums colSums
#' @importFrom  rlang       sym syms 
#' @importFrom  stats       setNames sd 
#' @importFrom  tibble      rownames_to_column column_to_rownames
#' @importFrom  tidyr       gather spread unite
#' @importFrom  DescTools   Gini
#' @importFrom  pbapply     pblapply pbsapply pboptions
NULL

# Package loading actions
.onAttach <- function(libname, pkgname) {
    msg <- paste("sargent package, a Sanofi production.",
                 "Confidential: internal use only.",
                 sep="\n")
    packageStartupMessage(msg)
}
