# Classes, Generics and Methods

#### Generics ####

setClassUnion("NumNULL", members=c("numeric", "NULL"))
setClassUnion("CharNULL", members=c("character", "NULL"))
setClassUnion("DfNULL", members=c("data.frame", "NULL"))
setClassUnion("MtxNULL", members=c("matrix", "NULL"))


#### Sargent classes ####

#' S4 class defining sargentObject class
#'
#' \code{sargentObject} defines sargentAnnotation outputs.
#' 
#' @slot  cells             vector of input cell ids.
#' @slot  cells_type        list of annotated cells.
#' @slot  cells_state       vector specifying each cell annotation status:
#'   \itemize{
#'            \item    \code{classified}: annotated cell. 
#'            \item    \code{classified-smooth}: annotated cell post-smoothing. 
#'            \item    \code{unclassified}: unannotated cell. 
#'            \item    \code{unclassified-smooth}: unannotated cell post-smoothing. 
#'            }
#' @slot  cells_score       numeric matrix of scores (gene-sets as rows, cells
#'                          as columns).
#' @slot  cells_gini        vector specifying gini-score of each cell.
#' @slot  gini_min          minimum threshold on the gini index.
#' @slot  threshold         numeric value specifying the threshold used to
#'                          idenfity unannotated cells
#' @slot  celltype_summary  data.frame containing a summary of annotation
#'                          associated statistics: 
#'  \itemize{ 
#'            \item    \code{celltype}: given cell-types for annotation. 
#'            \item    \code{gene}: given markers for annotation. 
#'            \item    \code{cell_count}: number of cells annotated with given 
#'                     cell-types. 
#'            \item    \code{cell_exp_count}: number of cells with detectable 
#'                     (>0) expression for given genes and cell-types. 
#'            \item    \code{pct_exp}: fractione of cells with detectable (>0) 
#'                     expression for given genes and cell-types. 
#'            \item    \code{mean_exp}: average expression for given genes and 
#'                     cell-types. 
#'            \item    \code{z_score}: z-score calculated per gene across 
#'                     cell-types. 
#'            }                     
#' 
#' @name         sargentObject-class
#' @rdname       sargentObject-class
#' @aliases      sargentObject
#' @exportClass  sargentObject
#' @return       sargentObject
setClass("sargentObject", 
         slots=c(cells="character",
                 cells_type="list",
                 cells_state="CharNULL",
                 cells_score="MtxNULL",
                 cells_gini="NumNULL",
                 gini_min="NumNULL",
                 threshold="NumNULL",
                 celltype_summary="DfNULL"))


#' S4 class defining scoreObject class
#'
#' \code{scoreObject} defines sargentAnnotation outputs if `only.score` is TRUE.
#' 
#' @slot  cells_score       numeric matrix of scores (gene-sets as rows, cells
#'                          as columns).
#' 
#' @name         scoreObject-class
#' @rdname       scoreObject-class
#' @aliases      scoreObject
#' @exportClass  scoreObject
#' @return       scoreObject
setClass("scoreObject", 
         slots=c(cells_score="matrix"))


#### Methods ####

#' @param    x    Sargent object
#' 
#' @rdname   sargentObject-class
#' @aliases  sargentObject-method
#' @export
setMethod("print", c(x="sargentObject"), function(x) { lengths(x@cells_type) })


