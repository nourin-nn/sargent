# Main function, helper functions, and summary finction for sargent

#' @include Classes.R
NULL


#### Main function ####

#' @title Annotate scRNA-seq data
#' @description   Returns the cell-type of origin for each cell, given a gene 
#'                expression matrix and sets of gene-markers
#'
#' @param    gex               A numeric matrix of single-cell expression values
#'                             where rows are genes and columns are cells. The 
#'                             matrix can be stored in sparse format (e.g., 
#'                             \code{count} slot from a \code{Seurat} object).
#' @param    gene.sets         A list of gene-markers. Each element of list is
#'                             named by the cell-type and contains a vector of 
#'                             associated gene-markers.
#' @param    gene.sets.neg     A list of gene-markers. Each element of list is
#'                             named by the cell-type and contains a vector of 
#'                             associated gene-markers.
#' @param    gini.min          A numeric scalar, between 0 and 1, specifying the
#'                             minimum threshold on the gini index
#' @param    gini.sigma        A numeric scalar specifying the number of sigmas
#'                             to use for identifying ambiguous calls (outliers) 
#'                             in the distribution of gini indexs.
#' @param    adjacent.mtx      A numeric matrix representing the nearest-neighbor 
#'                             graph. Matrix is binary, 1 for neigbors and 0 for
#'                             not-neighbors, and square with number of rows and 
#'                             coulmns equal to the number of cells.
#' @param    cells             Vector of cells to annotate (default is all
#'                             cells).
#' @param    n.neighbors       An intreger specifying the number of nearest
#'                             neighbors.
#' @param    score.threshold   A numeric scalar specifiying the minmum threshold 
#'                             on scores to annotate cells by a given cell-type. 
#'                             Only used if one set of gene-markers is given.
#' @param    only.score        if TRUE, only calculates assignment scores. 
#' 
#' @return   a \link{sargentObject}.
#'
#' @seealso  \link{sargentObject} for a sargent object.
#' 
#' @examples
#' # getting the data:
#' x <- get(data("sargentDemoData"))
#' 
#' # running cell-type assignment:
#' srgnt <- sargentAnnotation(gex=x$gex,
#'                            gene.sets=x$gene.sets,
#'                            adjacent.mtx=x$adjacent.mtx)
#' print(srgnt)                                 
#' table(predicted=fetchAssignment(srgnt), truth=x$mt.data$celltype)
#' 
#' @export
sargentAnnotation <- function(gex, gene.sets, 
                              gene.sets.neg=NULL, 
                              cells = NULL,
                              gini.min=0.5, gini.sigma=3, 
                              adjacent.mtx=NULL, n.neighbors=10,
                              score.threshold=NULL, only.score=FALSE) {
    start_time <- Sys.time()
  # ===================================
  # checks
  chk <- checks(gex=gex, gsets=gene.sets)
  if (chk != TRUE) stop(chk)
  # ===================================
  message("+++++++++++++++++++++++++++++++++++")
  dims_i <- dim(gex)
  if (!is.null(cells)) {
    # keep cells
    gex <- gex[, colnames(gex) %in% cells]
    dims_f <- dim(gex)
    message(dims_i[2] - dims_f[2], " cell(s) removed")
    # remove any genes with zero expression
    if (any(rowSums(gex) == 0)) {
      gex <- gex[rowSums(gex) != 0, ]
      dims_f <- dim(gex)
      message(dims_i[1] - dims_f[1], 
              " gene(s) removed. No expression across cells.")
    }
  }
  # ===================================
  message("A matrix with ", 
          format(dim(gex)[1], big.mark=","), " genes and ", 
          format(dim(gex)[2], big.mark=","), " cells.")
  cells.state <- setNames(rep("classified", dim(gex)[2]), colnames(gex))
  # ===================================
  # checks gene sets
  n_gsts <- length(gene.sets)
  message("gene markers: ")
  for (i in seq_along(gene.sets)) {
    gns <- gene.sets[[i]]
    message("* ", names(gene.sets)[i], ": ", sum(gns %in% rownames(gex)), "/", 
            length(gns), " gene(s) observed.")
    if (sum(gns %in% rownames(gex)) < length(gns)) {
      msg <- paste(gns[!gns %in% rownames(gex)], collapse=", ")
      message("-- Warning: gene(s) not present in the data: ", msg)
    }
  }
  # ===================================
  # checks
  chk <- checks(gex=gex, gsets=gene.sets)
  if (chk != TRUE) stop(chk)
  # ===================================
  cells.scrs <- scoreCells(mat=gex, 
                           gsets=gene.sets,
                           gsets.neg=gene.sets.neg) 
  # ===================================
  if (only.score) {
    # returns
    score_obj <- new("scoreObject",
                     cells_score=cells.scrs)
    message("+++++++++++++++++++++++++++++++++++")
    return(score_obj)
  }
  # ===================================
  if (n_gsts > 1) {
    helper.res <- annotByMultipleGeneSets(cells.scrs=cells.scrs, 
                                          cells.state=cells.state,
                                          gini.min=gini.min, 
                                          gini.sigma=gini.sigma,
                                          adjacent.mtx=adjacent.mtx, 
                                          n.neighbors=n.neighbors)
    assigns <- helper.res$assigns
    gini.thr <- helper.res$gini.thr
    gini.score <- helper.res$gini.score
    cells.state <- helper.res$cells.state
  } else {
    helper.res <- annotBySingleGeneSet(cells.scrs=cells.scrs, 
                                       score.thr=score.threshold)
    assigns <- helper.res$assigns
  }
  # =================================== 
  # print(table(assigns))
  stopifnot(length(assigns) == dim(gex)[2])
  # ===================================
  helper.res <- packingRes(gex=gex, gsets=gene.sets, assigns.vec=assigns)
  # =================================== 
  if (n_gsts == 1) {
    gini.score <- NULL
    gini.min <- NULL
  }
  # ===================================  
  # returns
  sargent_obj <- new("sargentObject",
                     cells=colnames(gex),
                     cells_type=helper.res$assigns.ls,
                     cells_state=cells.state,
                     cells_score=cells.scrs,
                     cells_gini=gini.score,
                     gini_min=gini.min,
                     threshold=ifelse(n_gsts > 1, gini.thr, score.threshold),
                     celltype_summary=helper.res$summ.df)
  # ===================================  
  end_time <- Sys.time()
  dt <- round(as.numeric(difftime(end_time, start_time, units = "mins")), 2)
  message(paste("time:", dt, "min"))
  message("+++++++++++++++++++++++++++++++++++")
  # ===================================  
  return(sargent_obj)
}


#### Helper functions ####


# Helper function for checkning inputs
checks <- function(gex, gsets) {
    if (is.null(colnames(gex))) {
        msg <- "No cell id. gex column names are empty."
        return(msg)
    }
    # ===================================  
    if (is.null(rownames(gex))) {
        msg <- "No gene name. gex row names are empty."
        return(msg)
    }
    # ===================================  
    genes <- unique(unlist(gsets, use.names = FALSE)) 
    if (sum(genes %in% rownames(gex)) == 0) {
        msg <- "None of the gene markers are present."
        return(msg)
    }
    # ===================================  
    obs.gsets <- lapply(gsets, function(x){ 
        x[x %in% rownames(gex)]
    })
    df <- plyr::ldply(obs.gsets, rbind) %>%
        tibble::column_to_rownames(var = ".id")
    dpcts <- rownames(df)[which(duplicated(df) | duplicated(df, fromLast = TRUE))]
    if (length(dpcts) > 0) {
        msg <- paste(dpcts, collapse=", ")
        msg <- paste("Identical gene sets. Check:", msg, sep=" ")
        return(msg)
    }
    # ===================================  
    # Return TRUE if all checks pass
    return(TRUE)
}


# Helper function for annotating cells using multiple gene sets
annotByMultipleGeneSets <- function(cells.scrs, cells.state,
                                    gini.min=0.5, gini.sigma=3,
                                    adjacent.mtx=NULL, n.neighbors=10) {
  assigns <- apply(cells.scrs, 2, which.max)
  assigns <- setNames(rownames(cells.scrs)[assigns], names(assigns))
  # print(table(assigns))
  stopifnot(length(assigns) == dim(cells.scrs)[2])
  # =================================== 
  trim.res <- trimAssignment(cells.scrs=cells.scrs, 
                             assigns.vec=assigns, 
                             cells.state=cells.state,
                             gini.min=gini.min, 
                             gini.sigma=gini.sigma)
  gini.score <- trim.res$gini.score
  gini.thr <- trim.res$gini.thr
  assigns <- trim.res$assigns.trimed
  cells.state <- trim.res$cells.state
  # ===================================
  stopifnot(all(names(cells.state[cells.state != "classified"]) %in% 
                    names(assigns[assigns == "unclassified"])))
  stopifnot(all(names(assigns[assigns == "unclassified"]) %in% 
                    names(cells.state[cells.state != "classified"])))
  # ===================================
  # smoothing
  if (!is.null(adjacent.mtx)) {
    stopifnot(dim(adjacent.mtx)[2] == length(assigns))
    smooth.res <- smoothAssigns(adjacent.mtx=adjacent.mtx,
                                assigns.vec=assigns,
                                cells.state=cells.state,
                                n.neighbors=n.neighbors)
    assigns <- smooth.res$assigns.smooth
    cells.state <- smooth.res$cells.state
    # ===================================  
    stopifnot(all(names(cells.state[!grepl("^classified", cells.state)]) %in% 
                      names(assigns[assigns == "unclassified"])))
    stopifnot(all(names(assigns[assigns == "unclassified"]) %in% 
                      names(cells.state[!grepl("^classified", cells.state)])))
  }
  # ===================================
  # returns 
  return(list("assigns"=assigns,
              "gini.thr"=gini.thr,
              "gini.score"=gini.score,
              "cells.state"=cells.state))
}


# Helper function for annotating cells using single gene set
annotBySingleGeneSet <- function(cells.scrs, score.thr=NULL) {
  # ===================================
  assigns <- names(which(cells.scrs[rownames(cells.scrs), ] > score.thr))
  assigns <- setNames(rep(rownames(cells.scrs), length(assigns)), assigns)
  # ===================================
  not_assigns <- setdiff(colnames(cells.scrs), names(assigns)) 
  not_assigns <- setNames(rep("unclassified", length(not_assigns)), not_assigns)
  assigns <- c(assigns, not_assigns)
  # print(table(assigns))
  stopifnot(length(assigns) == dim(cells.scrs)[2])
  # ===================================
  # returns 
  return(list("assigns"=assigns))
}


# Helper function to identify ambiguous annotations
trimAssignment <- function(cells.scrs, assigns.vec, cells.state,
                           gini.min=0.5, gini.sigma=3) {
  message("identifying ambiguous calls...")
  # Gini index < 0.2 represents perfect income equality, 
  # 0.2–0.3 relative equality, 0.3–0.4 adequate equality, 
  # 0.4–0.5 large inequality, and above 0.5 represents severe inequality.
  gini.score <- apply(cells.scrs, 2, DescTools::Gini)
  # ===================================
  gini.thr <- mean(gini.score, na.rm=TRUE) - gini.sigma*sd(gini.score, na.rm=TRUE)
  # ===================================
  x <- names(which(gini.score < gini.thr & gini.score < gini.min))
  assigns.vec[x] <- "unclassified"
  cells.state[x] <- "low-gini"
  # ===================================
  x <- names(which(is.na(gini.score)))
  assigns.vec[x] <- "unclassified" 
  cells.state[x] <- "no-marker"
  # ===================================
  x <- names(which(colSums(cells.scrs) == 0))
  assigns.vec[x] <- "unclassified"
  cells.state[x] <- "no-marker"
  # ===================================
  x <- names(which(gini.score == 0))
  assigns.vec[x] <- "unclassified"   
  cells.state[x] <- "multiplet"
  # ===================================
  x <- names(which(apply(cells.scrs, 2, function(x) sum(x == max(x))) > 1))
  x <- x[!is.na(gini.score[x])]
  assigns.vec[x] <- "unclassified"
  cells.state[x] <- "multiplet"
  # ===================================
  # print(table(assigns.vec))
  # print(table(cells.state))
  # ===================================
  return(list("gini.score"=gini.score,
              "gini.thr"=gini.thr,
              "assigns.trimed"=assigns.vec,
              "cells.state"=cells.state))
}


# Helper function to calculate scores for each cell
scoreCells <- function(mat, gsets, gsets.neg=NULL) {
  message("calculating assignment scores...")
  # maxrank <- ceiling(maxrank * nrow(mat))
  # set.seed(12345)
  cols_ls <- split(seq_len(ncol(mat)), ceiling(seq_len(ncol(mat))/1000))
  # ===================================
  op <- pboptions(style=3, char="=")
  # ===================================
  # sum(cumsum(head(y, maxrank) %in% gns)) 
  cells.scrs <- pblapply(cols_ls, function(x){
    apply(mat[, x], 2, function(y){
      y <- names(rev(sort(y[y > 0])))
      vapply(gsets, function(gns) {
        sum(cumsum(y %in% gns)) 
      }, FUN.VALUE=numeric(1))
    })
  })
  # ===================================
  # sum(cumsum(head(y, maxrank) %in% gns)) 
  if (!is.null(gsets.neg)) {
    pnlty.scrs <- pblapply(cols_ls, function(x){
      apply(mat[, x], 2, function(y){
        y <- names(rev(sort(y[y > 0])))
        vapply(gsets.neg, function(gns){
          sum(cumsum(y %in% gns)) 
        }, FUN.VALUE=numeric(1))
      })
    })
    # ===================================
    # sum(head(y, maxrank) %in% gns)
    bonus.scrs <- pblapply(cols_ls, function(x){
      apply(mat[, x], 2, function(y){
        y <- names(y[y == 0])
        vapply(gsets.neg, function(gns){
          sum(y %in% gns)
        }, FUN.VALUE=numeric(1))
      })
    })
  }
  # ===================================
  if (length(gsets) > 1) {
    cells.scrs <- do.call(cbind, cells.scrs)
    if (!is.null(gsets.neg)) {
      if (length(gsets.neg) == 1) {
        pnlty.scrs <- t(as.matrix(Reduce(c, pnlty.scrs)))
        rownames(pnlty.scrs) <- names(gsets.neg)
        bonus.scrs <- t(as.matrix(Reduce(c, bonus.scrs)))
        rownames(bonus.scrs) <- names(gsets.neg)
      } else {
        pnlty.scrs <- do.call(cbind, pnlty.scrs)
        bonus.scrs <- do.call(cbind, bonus.scrs)
      }
      stopifnot(all(colnames(cells.scrs) == colnames(pnlty.scrs)))
      stopifnot(all(colnames(cells.scrs) == colnames(bonus.scrs)))
      cells.scrs[rownames(pnlty.scrs), ] <- cells.scrs[rownames(pnlty.scrs), ] - pnlty.scrs
      cells.scrs[rownames(bonus.scrs), ] <- cells.scrs[rownames(bonus.scrs), ] + bonus.scrs
    }
    stopifnot(dim(cells.scrs)[2] == dim(mat)[2])    
  } else if (length(gsets) == 1) {
    names(cells.scrs) <- NULL
    cells.scrs <- unlist(cells.scrs)
    stopifnot(length(cells.scrs) == dim(mat)[2])    
  } else {
    stop()
  }
  # ===================================
  cells.scrs <- cells.scrs/max(abs(cells.scrs))
  if (any(cells.scrs < 0)) {
    cells.scrs <- (cells.scrs - min(cells.scrs)) / (max(cells.scrs) - min(cells.scrs))
  }
  # ===================================
  if (length(gsets) == 1) {
    cells.scrs <- t(as.matrix(cells.scrs))
    rownames(cells.scrs) <- names(gsets)      
  }
  # ===================================
  return(cells.scrs)
}


# Helper function to smmoth assignments
smoothAssigns <- function(adjacent.mtx, assigns.vec, 
                          cells.state, n.neighbors=10){
    message("smoothing assignments...")
    # make shure adjacent.mtx is symmetric 
    diag(adjacent.mtx) <- 0
    g <- graph_from_adjacency_matrix(adjacent.mtx, 
                                     mode="undirected", 
                                     weighted=NULL)
    adjM_ls <- as_adj_list(g, mode="all")
    # ===================================
    stopifnot(length(adjM_ls) == ncol(adjacent.mtx))
    # do not seek verdicts for cells with less than 10 neighbors
    adjM_ls[lengths(adjM_ls) < n.neighbors] <- NULL 
    # print(table(lengths(adjM_ls)))
    # ===================================
    op <- pboptions(style=3, char="=")
    # ===================================
    # seek above 50% verdicts
    verdicts <- pblapply(adjM_ls, function(x){
        names(which(table(assigns.vec[x])/length(x) > 0.5))
        })
    # print(table(lengths(verdicts)))
    verdicts[lengths(verdicts) == 0] <- NULL
    verdicts <- unlist(verdicts)
    ids <- verdicts == assigns.vec[names(verdicts)]
    # print(table(ids))
    switches <- names(ids[!ids])
    assigns.smooth <- assigns.vec
    assigns.smooth[switches] <- verdicts[switches]
    # ===================================
    logicx <- switches[verdicts[switches] != "unclassified"]
    cells.state[logicx] <- "classified-smooth"
    
    logicx <- switches[verdicts[switches] == "unclassified"]
    cells.state[logicx] <- "unclassified-smooth"
    # ===================================
    # print(table(assigns.smooth))
    # print(table(cells.state))  
    # ===================================
    return(list("assigns.smooth"=assigns.smooth,
                "cells.state"=cells.state))
}


# Helper function to pack results
packingRes <- function(gex, gsets, assigns.vec) {
  message("packing results...")
  
  # Hack for visibility of dplyr variables
  . <- NULL
  cell <- umi <- gene <- celltype <- z_score <- NULL
  pct_exp <- cell_count <- cell_exp_count <- mean_exp <- NULL
  
  ls1 <- split(names(assigns.vec), as.vector(assigns.vec))
  celltypes <- unique(as.vector(assigns.vec))
  genes <- unique(unlist(gsets, use.names = FALSE)) 
  gex <- gex[rownames(gex) %in% genes, , drop = FALSE]
  # gex.cpm <- Log1pNormGex(mat = gex)
  # gex.cpm <- log1p(1e6*t(t(gex)/colSums(gex))) ***
  if (!is.matrix(gex)) {
    gex <- as.matrix(gex)
  }
  suppressMessages(
    df <- data.frame(cell=colnames(gex)[col(gex)],
                     gene=rownames(gex)[row(gex)], 
                     umi=c(gex)) %>%
      dplyr::mutate(celltype=assigns.vec[!!rlang::sym("cell")]) %>%
      dplyr::group_by(!!rlang::sym("celltype")) %>%
        # number of cells in the cluster
      dplyr::mutate(cell_count=length(unique(!!rlang::sym("cell")))) %>% 
      dplyr::ungroup() %>%
      dplyr::group_by(!!!rlang::syms(c("celltype", "gene"))) %>%
      dplyr::summarise(cell_count=unique(!!rlang::sym("cell_count")), 
                       # cells with detectable (>0) expression of that gene
                       cell_exp_count=sum(!!rlang::sym("umi") > 0),      
                       # average experssion of all cells in the group
                       mean_exp=mean(!!rlang::sym("umi"))) %>%           
      dplyr::ungroup() %>%
      dplyr::mutate(pct_exp=!!rlang::sym("cell_exp_count")/!!rlang::sym("cell_count")) %>%
      dplyr::group_by(!!rlang::sym("gene")) %>%
        # scale across celltypes (gene scaled)
      dplyr::mutate(z_score=scale(mean_exp, center=TRUE, scale=TRUE)[,1]) %>% 
      dplyr::ungroup() %>%
      dplyr::mutate(gene=factor(gene, levels = genes, ordered = TRUE),
                    celltype=factor(celltype, levels=celltypes, ordered=TRUE)) %>%
      dplyr::select(celltype, gene, cell_count, cell_exp_count, 
                    pct_exp, mean_exp, z_score)  
  )
  # ===================================
  # returns 
  return(list("assigns.ls"=ls1,
              "summ.df"=df))
}


# Helper function to build base theme
getBaseTheme <- function() {
  # Define universal plot settings
  base_theme <- theme_bw() + 
    theme(text=element_text(size=10),
          plot.margin=grid::unit(c(0.5,0.5,0.5,0.5), units="line")) +
    theme(plot.title=element_text(size=12, hjust=0.5)) +
    theme(panel.border = element_rect(size = 1.25)) +
    theme(strip.background=element_blank(),
          strip.text=element_text(size=11, face = "bold")) +
    theme(axis.title=element_text(size=11, vjust=0.25, face = "bold"),
          axis.text=element_text(size=10, face = "bold")) +
    theme(legend.margin=margin(t=-0.1, r=0, b=0, l=0, unit="line"),
          legend.position='bottom',
          legend.title=element_text(size=10),
          legend.text=element_text(size=9),
          legend.key.height=grid::unit(10, "points"), 
          legend.key.width=grid::unit(15, "points"),
          legend.spacing=grid::unit(2, "points"))
  
  return(base_theme)
}


#### Summary functions ####

#' @title Annotation outcome visualization via dot-plot
#' @description   Displays how gene-markers are expressed across cell-types.
#'
#' @param    data       sargent object.
#' @param    genes      vector of genes to show. Must be chosen from the given
#'                      gene-markers.
#' @param    min.pct    the fraction of cells at which to draw the smallest dot
#'                      (default is 0).
#'   
#' @return   A ggplot object.
#'
#' @seealso  \link{sargentObject} for a sargent object.
#'
#' @examples
#' # getting the data:
#' x <- get(data("sargentDemoData"))
#' 
#' # running cell-type assignment:
#' srgnt <- sargentAnnotation(gex=x$gex,
#'                            gene.sets=x$gene.sets,
#'                            adjacent.mtx=x$adjacent.mtx)
#'
#' fetchDotPlot(srgnt, 
#'              genes = c("CD3D","CD3E","IL7R","ANPEP","SLC13A3","SLC17A3",
#'                        "EMCN","NOTCH4","ADGRL4","CD79A","MS4A1","CD19",
#'                        "CD14","MS4A7","FCGR3A"), 
#'              min.pct=0.1)
#' 
#' @export
fetchDotPlot <- function(data, genes=NULL, min.pct=0) {
  
  # Hack for visibility of dplyr variables
  . <- NULL
  gene <- pct_exp <- z_score <- celltype <- NULL
  
  x <- data@celltype_summary %>%
    dplyr::filter(pct_exp > min.pct) 
  
  if (!is.null(genes)) {
    x <- x %>%
      dplyr::filter(gene %in% genes) %>%
      dplyr::mutate(gene = factor(gene, levels=genes, ordered=TRUE))
  }
  
  p <- ggplot(x,
              aes(x=celltype, y=gene, color=z_score, size=pct_exp*100)) +
    getBaseTheme() + 
    theme(axis.title=element_blank(),
          legend.position="right",
          legend.box="vertical",
          legend.key.width=grid::unit(20, "points"),
          legend.key.height=grid::unit(20, "points"),
          legend.spacing=grid::unit(15, "points"),
          panel.grid.major=element_line(size=0.1),
          axis.text.x=element_text(angle=90, hjust=1, vjust=0.5)
    ) +
    geom_point() + 
    scale_size(name="Percent Expressing",
               range=c(0,6), breaks=seq(0,100,25)) +
    scale_color_gradient(low="lightgrey", high="firebrick", 
                         # limits=c(-0.5,2.5), oob=scales::squish,
                         name='Average Expression') +
    coord_flip()
  return(p)
}


#' @title Annotation confidency visualization via density-plot
#' @description   Displays distribution of gini coefficients calculated for each
#'                cell across cell-types.
#'
#' @param    data       sargent object.
#'
#' @return   A ggplot object.
#'
#' @details A gini coefficient above 0.5 represents good confidency in
#' assignments. The vertical dashed line represents the threshold used to
#' identfiy assignemnts with low confidency. If one gene-set is given for
#' annotation, the gini coefficients will not be calculated and the distribution
#' of cells' scores are shown instead. In this scenario, the vertical dashed
#' line represents the threshold used to identfiy cells that belong to the given
#' gene-set.
#' 
#' @seealso  \link{sargentObject} for a sargent object.
#' 
#' @examples
#' # getting the data:
#' x <- get(data("sargentDemoData"))
#' 
#' # running cell-type assignment:
#' srgnt <- sargentAnnotation(gex=x$gex,
#'                            gene.sets=x$gene.sets,
#'                            adjacent.mtx=x$adjacent.mtx)
#'                            
#' fetchDensityPlot(srgnt)
#'                            
#' @export
fetchDensityPlot <- function(data) {
  
  if (inherits(data, "sargentObject")) {
    
    if (nrow(data@cells_score) > 1) {
      x <- data@cells_gini
      xlabel <- "Gini coefficient"
    } else {
      x <- data@cells_score
      xlabel <- "Cell score"
    }
    
    xntrcpt <- data@threshold
    p <- ggplot(data=data.frame(x=x[!is.na(x)]), 
                aes(x=x)) +
      getBaseTheme() +
      theme(axis.text=element_text(size=10, face="bold"),
            axis.title=element_text(size=12, face="bold"),
            axis.ticks.y=element_blank(),
            panel.grid.major=element_blank(),
            panel.grid.minor=element_blank()) +
      xlab(xlabel) + ylab("Density") + 
      coord_cartesian(xlim=c(0, 1), ylim=c(0, 1)) +
      ggdist::stat_slab(alpha=.55, fill= "steelblue") +
      ggdist::stat_pointinterval(stroke=1) + 
      geom_vline(xintercept=xntrcpt, linetype=2, size=0.5, color="firebrick")
    
    if (!is.null(data@gini_min)) {
      p <- p + geom_vline(xintercept=data@gini_min, linetype=2, size=0.5, color="black")
    }    
  }

  if (inherits(data, "scoreObject")) {
    
    # Hack for visibility of dplyr variables
    cell <- celltype <- score <- NULL
    
    x <- data@cells_score
    df <- data.frame(cell=colnames(x)[col(x)],
                     celltype=rownames(x)[row(x)], 
                     score=c(x))
    
    xlabel <- "Cell score"
    p <- ggplot(data=df, 
                aes(x=score)) +
      getBaseTheme() +
      theme(axis.text=element_text(size=10, face="bold"),
            axis.title=element_text(size=12, face="bold"),
            axis.ticks.y=element_blank(),
            panel.grid.major=element_blank(),
            panel.grid.minor=element_blank()) +
      xlab(xlabel) + ylab("Density") + 
      coord_cartesian(xlim=c(0, 1), ylim=c(0, 1)) +
      ggdist::stat_slab(alpha=.55, fill= "steelblue") +
      ggdist::stat_pointinterval(stroke=1)
  }
  
  return(p)
}


#' @title Access gini coefficients
#' @description   Retrieves gini coefficients in a Sargent Object.
#'
#' @param    data          sargent object.
#' @param    cells         a vector of cells to collect data. Default is all
#'                         cells.
#' @param    cell.types    a vector of cell-types to collect data. Default is
#'                         all cell-types.
#' 
#' @return   A \code{data.frame}. 
#' 
#' @seealso  \link{sargentObject} for a sargent object.
#' 
#' @examples
#' # getting the data:
#' x <- get(data("sargentDemoData"))
#' 
#' # running cell-type assignment:
#' srgnt <- sargentAnnotation(gex=x$gex,
#'                            gene.sets=x$gene.sets,
#'                            adjacent.mtx=x$adjacent.mtx)
#' 
#' df <- fetchGiniIndex(srgnt)
#' head(df)
#'                             
#' @export
fetchGiniIndex <- function(data, cells=NULL, cell.types=NULL) {
  
  # Hack for visibility of dplyr variables
  . <- NULL
  
  x <- fetchAssignment(data)
  kp.cells <- names(x)
  if (any(c(!is.null(cells), !is.null(cell.types)))) {
    kp.cells <- kp.cells[kp.cells %in% unique(c(cells, names(x[as.vector(x) %in% 
                                                                   cell.types])))]  
  }
  
  x <- data@cells_gini
  x <- x[names(x) %in% kp.cells]
  cstate <- data@cells_state
  x <- as.data.frame(x) %>%
    tibble::rownames_to_column(var="cell") %>%
    dplyr::rename(gini=x) %>%
    dplyr::mutate(state=cstate[.$cell])
  return(x)
}


#' @title Access assignment scores
#' @description   Retrieves assignment scores in a Sargent Object.
#'
#' @param    data          sargent object.
#' @param    cells         a vector of cells to collect data. Default is all
#'                         cells.
#' @param    cell.types    a vector of cell-types to collect data. Default is
#'                         all cell-types.
#' 
#' @return   A \code{data.frame}. 
#' 
#' @seealso  \link{sargentObject} for a sargent object.
#' 
#' @examples
#' # getting the data:
#' x <- get(data("sargentDemoData"))
#' 
#' # running cell-type assignment:
#' srgnt <- sargentAnnotation(gex=x$gex,
#'                            gene.sets=x$gene.sets,
#'                            adjacent.mtx=x$adjacent.mtx)
#' 
#' df <- fetchAssignmentScore(srgnt)
#' head(df)
#' 
#' @export
fetchAssignmentScore <- function(data, cells=NULL, cell.types=NULL) {
  
  # Hack for visibility of dplyr variables
  . <- NULL
  
  x <- fetchAssignment(data)
  kp.cells <- names(x)
  if (any(c(!is.null(cells), !is.null(cell.types)))) {
    kp.cells <- kp.cells[kp.cells %in% 
                             unique(c(cells, names(x[as.vector(x) %in% 
                                                         cell.types])))]  
  }
  
  x <- data@cells_score
  x <- x[, colnames(x) %in% kp.cells, drop=FALSE]
  cstate <- data@cells_state
  x <- as.data.frame(t(x)) %>%
    tibble::rownames_to_column(var="cell") %>%
    dplyr::mutate(state=cstate[.$cell])
  return(x)
}


#' @title Access cell-type assignments
#' @description    Retrieves cell-type assignments in a Sargent Object.
#'
#' @param    data          sargent object.
#' 
#' @return   A named vector of assignments. 
#' 
#' @seealso  \link{sargentObject} for a sargent object.
#' 
#' @examples
#' # getting the data:
#' x <- get(data("sargentDemoData"))
#' 
#' # running cell-type assignment:
#' srgnt <- sargentAnnotation(gex=x$gex,
#'                            gene.sets=x$gene.sets,
#'                            adjacent.mtx=x$adjacent.mtx)
#' 
#' vec <- fetchAssignment(srgnt)
#' table(vec)
#' 
#' @export
fetchAssignment <- function(data) {
  x <- data@cells_type
  x <- setNames(unlist(x, use.names=FALSE), rep(names(x), lengths(x)))
  x <- setNames(names(x), x)
  x <- x[data@cells]
  return(x)
}


#' @title Access cell-type assignments
#' @description   Retrieves cell-type assignments for specific cells in a
#'   Sargent Object.
#'
#' @param    data          sargent object.
#' @param    cells         a vector of cells to collect data. Default is all
#'                         cells.
#' @param    form          indicates the shape of out puts. See \code{details}
#'
#' @return   A \code{data.frame}.
#'
#' @details If \code{form} = \code{gather} returns a \code{data.frame} with two
#'   columns indicating cell ids and associated assignments, respectively. If
#'   \code{form} = \code{spread} returns a \code{data.frame} in which the first
#'   column contains cell ids, next columns contain the associated assignments
#'   of each cell in a binary format (0/1) for given cell-types, and the last
#'   column contains information about the state of the assignemts.
#'
#' @seealso  \link{sargentObject} for a sargent object.
#'
#' @examples
#' # getting the data:
#' x <- get(data("sargentDemoData"))
#' 
#' # running cell-type assignment:
#' srgnt <- sargentAnnotation(gex=x$gex,
#'                            gene.sets=x$gene.sets,
#'                            adjacent.mtx=x$adjacent.mtx)
#'
#' cells <- colnames(x$gex)[seq_len(10)]
#' df <- fetchCelltype(srgnt, cells=cells)
#' head(df)
#'
#' @export
fetchCelltype <- function(data, cells=NULL, form=c("gather", "spread")) {
  
  # Hack for visibility of dplyr variables
  . <- NULL
  cell <- yesno <- celltype <- NULL

  # Check arguments
  form <- match.arg(form)
  
  x <- fetchAssignment(data)
  if (!is.null(cells)) {
    x <- x[names(x) %in% cells]
  }
  
  x <- as.data.frame(x) %>%
    tibble::rownames_to_column(var="cell") %>%
    dplyr::rename(celltype=x)
  
  if (form == "spread") {
    cstate <- data@cells_state
    x <- x %>%
      dplyr::mutate(yesno=1) %>%
      dplyr::distinct() %>%
      tidyr::spread(key=celltype, value=yesno, fill=0) %>%
      dplyr::mutate(state=cstate[.$cell])
  }
  return(x)
}


#' @title Access cells
#' @description   Retrieves cells assigned by specific cell-types in a Sargent
#'                Object.
#'
#' @param    data          sargent object.
#' @param    cell.types    a vector of cell-types to collect data. Default is
#'                         all cell-types.
#' @param    form          indicates the shape of out puts. See \code{details}
#' 
#' @return   A \code{data.frame}. 
#' 
#' @details  
#' If \code{form} = \code{gather} returns a \code{data.frame} with two columns
#' indicating cell ids and associated assignments, respectively. If \code{form}
#' = \code{spread} returns a \code{data.frame} in which the first column
#' contains cell ids, next columns contain the associated assignments of each
#' cell in a binary format (0/1) for given cell-types, and the last column
#' contains information about the state of the assignemts.
#' 
#' @seealso  \link{sargentObject} for a sargent object.
#' 
#' @examples
#' # getting the data:
#' x <- get(data("sargentDemoData"))
#' 
#' # running cell-type assignment:
#' srgnt <- sargentAnnotation(gex=x$gex,
#'                            gene.sets=x$gene.sets,
#'                            adjacent.mtx=x$adjacent.mtx)
#' 
#' df <- fetchCell(srgnt, cell.types="TNK")
#' head(df)
#' 
#' @export
fetchCell <- function(data, cell.types=NULL, form=c("gather", "spread")) {
  
  # Hack for visibility of dplyr variables
  . <- NULL
  cell <- yesno <- celltype <- NULL
  
  # Check arguments
  form <- match.arg(form)
  
  x <- fetchAssignment(data)
  if (!is.null(cell.types)) {
    x <- x[as.vector(x) %in% cell.types]
  }
  
  cstate <- data@cells_state
  if (form == "gather") {
    x <- as.data.frame(x) %>%
      tibble::rownames_to_column(var="cell") %>%
      dplyr::rename(celltype=x) %>%
      dplyr::mutate(state=cstate[.$cell])        
  } else if (form == "spread") {
    x <- as.data.frame(x) %>%
      tibble::rownames_to_column(var="cell") %>%
      dplyr::rename(celltype=x) %>%
      dplyr::mutate(yesno=1) %>%
      dplyr::distinct() %>%
      tidyr::spread(key=celltype, value=yesno, fill=0) %>%
      dplyr::mutate(state=cstate[.$cell])        
  }
  return(x)
}


#' @title Merge sargent objects
#' @description   Merge sargent objects from mutiple level of annotating.
#'
#' @param    ...             sargent objects to merge.
#' @param    srgnt.list      a list of sargent objects to merge.
#' 
#' @return   A \code{data.frame} with two columns: 1) cell ids (\code{cell}) and
#'           cell-type (\code{celltype}). 
#' 
#' @details  
#' A common strategy to annotate single-cell RNA-seq data is to annotate cells
#' according to hierarchy of known cell-types. Then, given multiple sargent
#' objects made in each level of hierarchy, \code{fetchMergeSargents} merges
#' cell-type assignments based on the cell ids into one column
#' (\code{celltype}). column \code{celltype} contains assignmets from given
#' objects separated by ".".
#'  
#' @seealso  \link{sargentObject} for a sargent object.
#' 
#' @examples 
#' # getting the data:
#' x <- get(data("sargentDemoData")) 
#' # annotating single cell RNA-seq example data 
#' 
#' # running cell-type assignment:
#' srgnt <- sargentAnnotation(gex=x$gex,
#'                            gene.sets=x$gene.sets,
#'                            adjacent.mtx=x$adjacent.mtx)                                   
#' 
#' # get cells annotated as T or Natural Killer cells (TNKs)
#' cells <- fetchCell(srgnt, cell.types="TNK")$cell
#'
#' # annotate TNK subset  
#' srgnt_tnk <- sargentAnnotation(gex=x$gex, cells=cells,
#'                                gene.sets=x$gene.tnk.sets,
#'                                adjacent.mtx=x$adjacent.tnk.mtx)
#' print(srgnt_tnk)
#' 
#' @export
fetchMergeSargents <- function(..., srgnt.list=NULL){
  # Hack for visibility of dplyr variables
  . <- NULL
  cell <- celltype <- NULL
  
  srgnts <- c(list(...), srgnt.list)
  if(length(srgnts) == 0) return(NULL)
  
  # fetch cell-types
  srgnts <- lapply(srgnts, function(x){
    fetchCelltype(x)
  }) 
  
  # Create own merging function
  my_merge <- function(df1, df2){                                
    merge(df1, df2, by="cell", all=TRUE)
  }
  
  # Combine celltype columns
  srgnts <- Reduce(my_merge, srgnts) %>%
    tidyr::unite(., col="celltype", sep=".", colnames(.)[colnames(.) != "cell"], 
                 na.rm=TRUE)
  
  # return
  return(srgnts)
}
