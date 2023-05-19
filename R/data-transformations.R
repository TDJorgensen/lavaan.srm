### Terrence D. Jorgensen
### Last updated: 19 May 2023
### (currently hidden) functions to format round-robin data
### between wide, long, and matrix



#TODO: allow passing covariates at different levels

## transform from long-long (univariate outcome) format to square format,
## like an adjacency matrix / sociometric data
## - y = vector of round-robin variable names/indices in "data"
## - IDout, IDin = actor & partner (perceiver & target) ID names/indices
## - group = optional group ID name/index
## - returnList indicates whether multiple groups should be treated as 1 group (FALSE)
longUni2square <- function(y, data, IDout, IDin, group = NULL,
                           returnList = TRUE) {
  ## all variable names in the data?
  if (is.numeric(y)) {
    y <- colnames(data)[y]
  } else if (is.character(y)) {
    y.miss <- setdiff(y, colnames(data))
    if (length(y.miss)) stop('The following variable names in "y" are missing ',
                             'from the data set:\n', paste(y.miss, collapse = ","))
  } else stop('"y" must be a character vector indicating variable name(s) or',
              ' a numeric vector of indices in "data"')

  if (is.numeric(IDout)) IDout <- colnames(data)[IDout]
  if (is.numeric(IDin)) IDin <- colnames(data)[IDin]
  if (is.numeric(group)) group <- colnames(data)[group]
  id.miss <- setdiff(c(IDout, IDin, group), colnames(data))
  if (length(id.miss)) stop('The following ID variable names are missing from ',
                            'the data set:\n', paste(id.miss, collapse = ","))
  data <- as.data.frame(data[ , c(IDout, IDin, group, y)], stringsAsFactors = FALSE)

  ## unique IDs across groups
  #FIXME: necessary?  mvsrm() checks for duplicate case IDs.
  #       This overwrites user's IDs before mvsrm() saves them in @ID slot
  if (!is.null(group)) {
    gFac <- factor(data[ , group])

    if (!returnList) {
      ## unique IDs across groups
      data[ , group] <- as.integer(gFac)
      data[ , IDout] <- paste0("g", data[ , group], "_", data[ , IDout])
      data[ , IDin ] <- paste0("g", data[ , group], "_", data[ , IDin ])
    }
  }

  ## list of data frames (one per group)
  dataList <- if (is.null(group) || !returnList) list(data) else {
    lapply(levels(gFac), function(g) {
      data[data[,group] == g, ]
    })
  }

  ## loop over groups to store y in square matrix
  outList <- lapply(dataList, function(gdat) {
    ## all person/case-level IDs
    gdat[ , IDout] <- as.character(gdat[ , IDout])
    gdat[ , IDin] <- as.character(gdat[ , IDin])
    allIDs <- unique(c(gdat[ , IDout], gdat[ , IDin]))
    ## empty array of matrices
    gmat <- array(NA, dim = c(length(allIDs), length(allIDs), length(y)),
                  dimnames = list(ego = allIDs, alter = allIDs, variable = y))
    attr(gmat, "warn") <- 0L
    ## loop over rows of data (within each y) to copy values to matrix
    warn <- FALSE
    for (yy in y) for (n in 1:nrow(gdat)) {
      ## check if there was already a value
      if (yy == y[1] && !is.na(gmat[gdat[n, IDout], gdat[n, IDin], yy]))
        attr(gmat, "warn") <- attr(gmat, "warn") + 1L
      gmat[gdat[n, IDout], gdat[n, IDin], yy] <- gdat[n, yy]
    }
    gmat
  })
  #TODO: add diagnostics, or suggest using checkDupUni()

  ## check list for multiple values in a cell
  multiVals <- sapply(outList, attr, which = "warn")
  if (!all(multiVals == 0L)) {
    for (G in which(multiVals > 0L)) {
      attr(outList[[G]], "warn") <- paste('  Group', G, 'had', multiVals[G],
                                          'cell(s) with > 1 elligible value.')
    }
    warning('\n\n"data" contained multiple values for some cells in a matrix:\n',
            paste(sapply(outList[multiVals > 0L], attr, which = "warn"),
                  collapse = "\n"),
            '\nIn each instance, the last value was kept. Please check that ',
            'each subject only provides one datum about each other subject in ',
            'data[ , y].  Store repeated measures as additional variables.\n\n')
  }
  for (G in which(multiVals == 0L)) attr(outList[[G]], "warn") <- NULL

  ## convert arrays to 2-D matrices if univariate
  if (length(y) == 1L) for (G in seq_along(outList)) outList[[G]] <- outList[[G]][ , , 1]

  ## return result (remove list if 1 group)
  if (is.null(group) || !returnList) {
    return(outList[[1]])
  } else {
    names(outList) <- levels(gFac)
  }
  outList
}


## transform from long format (bivariate outcome) to square format,
## like an adjacency matrix / sociometric data
# longBi2square <- function(Ys, data, ID1, ID2, IDdyad = NULL, group = NULL) {}


## transform from square (adjacency matrix) to long (bivariate outcome) format
## - mat = a 2-D matrix or 3-D array (third dim = variable)
## - dropNA indicates BOTH dyadic observations must be NA to drop that row
##          MUST BE IGNORED when "mat" is a 3-D array (patterns can differ)
## - label = the round-robin variable name, used as a prefix for both columns
##           when mat is a 2-D matrix (otherwise, use names of third dimension)
## - suffix = labels for first and second cases in dyad, used for ID names
##            and to distinguish between columns of round-robin variable
## - group = label used for the grouping variable when is.list(mat)
## - groupMC = whether to group mean center each RR variable
##' @importFrom utils combn
square2longBi <- function(mat, dropNA = FALSE, label = "y", group = "group",
                          suffix = c("i","j"), groupMC = FALSE) {
  if (is.list(mat)) {
    isMat <- sapply(mat, function(m) is.array(m) && is.numeric(m))
    eqDim <- sapply(mat, function(m) diff(dim(m)[1:2]) == 0)
    eqLab <- sapply(mat, function(m) all.equal(dimnames(m)[[1]], dimnames(m)[[2]]))
    if (!all(isMat)) stop('Every element in the list "mat" must be a numeric matrix/array')
    if (!all(eqDim)) stop('Every element in the list "mat" must be a square matrix,',
                          ' or the first 2 dimensions of each 3-D array must match')
    if (!all(eqLab)) stop('row and column names must match in each square matrix')
    multigroup <- TRUE
  } else {
    if (!(is.array(mat) && is.numeric(mat)))
      stop('"mat" must be a numeric matrix/array (or list of them)')
    if (diff(dim(mat)[1:2]) != 0)
      stop('"mat" must be a square matrix (or list of them), or the first 2',
           ' dimensions of a 3-D array must match')
    if ( !all.equal(dimnames(mat)[[1]], dimnames(mat)[[2]]) )
      stop('rownames(mat) must == colnames(mat)')
    mat <- list(mat)
    multigroup <- FALSE
  }
  if (!is.character(suffix) || length(suffix) != 2L)
    stop('"suffix" must be a length-2 character vector')

  outList <- sapply(mat, function(M) {
    ## save each pair of IDs
    RN <- rownames(M)
    if (is.null(RN)) RN <- 1:nrow(M)
    IDs <- combn(RN, m = 2, simplify = FALSE,
                 FUN = function(x) c(ID_A = x[1], ID_B = x[2]))

    ## save each dyadic (i.e., bivariate) outcome for each pair of IDs
    if (length(dim(M)) == 2L) {
      nYs <- 1L # 1 round-robin variable
      label <- as.character(label[1]) # use "label=" argument for variable names
      yLabels <- paste0(label, "_", suffix[1:2], suffix[2:1])

      ## (group-)mean center?
      if (groupMC) {
        Mcopy <- M
        diag(Mcopy) <- NA # ignore any self-reports in diagonal
        M <- M - mean(Mcopy, na.rm = TRUE)
      }

      ## for each dyad (i != j)
      yList <- utils::combn(RN, m = 2, simplify = FALSE,
                            FUN = function(x) c(AB = M[ x[1] , x[2] ],
                                                BA = M[ x[2] , x[1] ]))
      ## check for missing?
      if (dropNA) {
        bothNA <- sapply(yList, function(x) all(is.na(x)))
        IDs <- IDs[which(!bothNA)]
        yList <- yList[which(!bothNA)]
      }
      Ys <- do.call(rbind, yList)


      ## if multivariate, loop over third dimension
    } else if (length(dim(M)) == 3L) {
      nYs <- dim(M)[[3]] # how many round-robin variables?
      if (nYs == length(label)) {
        ## use "label=" argument for variable names
        yLabels <- paste0(rep(as.character(label), each = 2),
                          "_", suffix[1:2], suffix[2:1])

      } else if (!is.null(dimnames(M)[[3]])) {
        ## use list names to make variable names
        yLabels <- paste0(rep(dimnames(M)[[3]], each = 2),
                          "_", suffix[1:2], suffix[2:1])

      } else {
        ## use the first "label=" argument and add numbers
        yLabels <- paste0(as.character(label[1]), rep(1:nYs, each = 2),
                          "_", suffix[1:2], suffix[2:1])

      }

      yList <- list()
      if (is.null(dimnames(M)[[3]])) dimnames(M)[[3]] <- 1:(dim(M)[3])
      for (y in dimnames(M)[[3]]) {
        M2 <- M[ , , y, drop = TRUE] # extract slice for this RR variable

        ## (group-)mean center?
        if (groupMC) {
          Mcopy <- M2
          diag(Mcopy) <- NA # ignore any self-reports in diagonal
          M2 <- M2 - mean(Mcopy, na.rm = TRUE)
        }

        ## for each dyad (i != j)
        tempY <- combn(RN, m = 2, simplify = FALSE,
                       FUN = function(x) c(AB = M2[ x[1] , x[2] ],
                                           BA = M2[ x[2] , x[1] ]))
        yList[[y]] <- do.call(rbind, tempY)
      }
      Ys <- do.call(cbind, yList)

    } else stop('mat= argument must be a (list of) 2- or 3-dimensional array(s).')

    ## combine each dyad's information into a data.frame
    out <- data.frame(do.call(rbind, IDs), Ys, stringsAsFactors = FALSE)
    ## assemble requested variable names
    colnames(out)[1:2] <- paste("ID", suffix, sep = "_")
    colnames(out)[1:(2*nYs) + 2L] <- yLabels
    ## return result
    out
  }, simplify = FALSE)

  if (multigroup) {
    ## add group labels
    gNames <- names(mat)
    if (is.null(gNames)) gNames <- seq_along(mat)
    for (g in gNames) outList[[g]][ , group] <- g
  }

  result <- do.call(rbind, outList)
  if (multigroup) result <- result[ , c(group, setdiff(colnames(result), group))]
  rownames(result) <- NULL
  result
}


## TODO: extract case-level data from diagonal of each RR-variable matrix
# square2case <- function(mat) {}

## transform from square (adjacency matrix) to long-long (univariate) format
# square2longUni <- function(mat, dropNA = TRUE) {}


## transform from bivariate (long) to univariate (long-long) format
# longBi2longUni <- function(Ys, data, ID1, ID2, IDdyad = NULL, group = NULL) {}


## transform from univariate (long-long) to bivariate (long) format in 2 steps:
## longUni2square() then square2longBi()
longUni2longBi <- function(y, data, IDout, IDin, group = NULL, groupMC = FALSE,
                           returnList = TRUE, dropNA = FALSE,
                           label = NULL, suffix = c("i","j")) {
  ## convert to squares
  SQ <- longUni2square(y = y, data = data, IDout = IDout, IDin = IDin,
                       group = group, returnList = returnList)
  ## group ID lost, so recover the name from original data
  if (is.numeric(group)) group <- colnames(data)[group]
  ## convert to long (bivariate) format
  square2longBi(mat = SQ, dropNA = dropNA, suffix = suffix,
                group = group, groupMC = groupMC,
                label = if (is.null(label)) y else label)
}


