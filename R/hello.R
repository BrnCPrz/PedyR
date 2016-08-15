#' Checks for errors and conformities in a pedigree file
#'
#' This function checks for common errors in genealogical structured data
#'
#' @param pedigree a data.frame with three columns: id, sire and dam
#' @param id integer containning the column where animal ids are located in the pedigree
#' @param sire integer containning the column where sire ids are located in the pedigree
#' @param dam integer containning the column where dam ids are located in the pedigree
#' @examples
#' id = c(1,2,3,4,5,6,7,8,9,10)
#' sire = c(0,0,1,1,3,5,5,5,7,7)
#' dam = c(0,0,2,2,2,4,4,6,6,8)
#' ped = as.data.frame(cbind(id,sire,dam))
#' checkPed(ped)
#'
#' @export
checkPed <- function (pedigree, id=1, sire=2, dam=3, rmsingle = FALSE, verbose = FALSE) {
  options(warn=-1)

  if (missing(pedigree))
    stop("  Pedigree file not specified!  ")

  if (!is.data.frame(pedigree))
    stop (" Pedigree file must be a data.frame object! ")

  anim = unique(as.character(pedigree[,id]))
  sires = unique(as.character(pedigree[,sire]))
  dams = unique(as.character(pedigree[,dam]))

  #check for duplicated ids
  if (length(anim) != length(pedigree[,id]))
  {
    dupid = pedigree[,1][duplicated(pedigree[,1])]
    cat(" Duplicated ids detected ! \n")
    cat(" -> ", (dupid[!duplicated(dupid)]), "\n")
    cat(" Duplicated ids position are represented by -1! \n")
    pedigree[,1][duplicated(pedigree[,1])] = -1
    return(pedigree[,1])
  }

  sireclean = sires[sires!="0"]
  damclean = dams[dams!="0"]
  wrongsex = sireclean %in% damclean

  if(TRUE %in% wrongsex) {
    stop (" Sire = Dam  ")
  }

  # remove singleton animals (no parents or progeny in the pedigree)
  if(rmsingle == TRUE){
      singleton.vec <- pedigree[,1][pedigree[,2] == 0 & pedigree[,3]== 0 & !(pedigree[,1] %in% c(pedigree[,2], pedigree[,3]))]

    if(length(singleton.vec) > 0) {
        pedigree = pedigree[-(pedigree[,1] %in% singleton.vec),]
    if (verbose == TRUE) cat("A total of", length(singleton.vec), "found and removed. \n")
      }
  }

  # check for individuals in sire/dam that are not in id and add them as founders
  if(any(!pedigree[,2] %in% pedigree[,1]))
    sirenew <- data.frame(id = pedigree[which(!pedigree[,2] %in% pedigree[,1]), 2],
                                                                sire = 0, dam = 0)
  if(any(!pedigree[,3] %in% pedigree[,1]))
    damnew <- data.frame(id = pedigree[which(!pedigree[,3] %in% pedigree[,1]), 3],
                                                                sire = 0, dam = 0)

  sirenew <- unique(subset(sirenew, id != 0))
  damnew <- unique(subset(damnew, id != 0))

  correctPed <- as.data.frame(rbind(sirenew, damnew, pedigree))

  addFounders = FALSE

  if(length(correctPed) > length(pedigree)) {
    addFounders == TRUE
    cat("Individuals appearing as sire/dam but not as individuals were added as founders. \n")
  }

  # wrapps up pedigree stats and output information
  cat(paste('Total number of animals = ', length(anim),sep=''),'\n')
  #cat('  > Animal ids: ', anim, '\n')
  cat(paste('Total number of unique sires = ', length(sireclean),sep=''),'\n')
  #cat('  > Sire ids: ', sireclean, '\n')
  cat(paste('Total number of unique dams = ', length(damclean),sep=''),'\n')
  #cat('  > Dam ids: ', damclean, '\n')

  cat('Pedigree is OK! \n')
  if (addFounders == TRUE) return(correctPed)
  else return(pedigree)
}


#' Calculates inbreeding coefficients of individuals in the pedigree.
#'
#' @param pedigree data.frame with three columns: id, sire and dam
#'
#' @details Uses Rcpp package for calculations
#' @return Numeric vector
#' @examples
#' id = c(1,2,3,4,5,6,7,8,9,10)
#' sire = c(0,0,1,1,3,5,5,5,7,7)
#' dam = c(0,0,2,2,2,4,4,6,6,8)
#' ped = as.data.frame(cbind(id,sire,dam))
#' ped$F = getF(ped)
#'
#' @export
getF <- function(pedigree){
  if (nargs()==0){
    stop(" Pedigree file not specified! ")
  }

  if(!is.data.frame(pedigree)){
    stop("Pedigree must be a data.frame object!")
  }

  s = pedigree[,2]
  d = pedigree[,3]

  n = length(s)
  N = n + 1
  A = matrix(0, ncol=N, nrow=N)

  s = (s == 0)*(n) + (s-1)
  d = (d == 0)*n + (d-1)

  aMatrix(A, s, d)
  #return(A[1:n, 1:n])
}


#' Calculates generations for each animal in the pedigree
#'
#' @param pedigree data.frame with three columns: id, sire and dam
#'
#' @details Uses Rcpp package for calculations
#'
#' @return Integer vector
#' @examples
#' id = c(1,2,3,4,5,6,7,8,9,10)
#' sire = c(0,0,1,1,3,5,5,5,7,7)
#' dam = c(0,0,2,2,2,4,4,6,6,8)
#' ped = as.data.frame(cbind(id,sire,dam))
#' ped$Gen = getGen(ped)
#'
#' @export
getGen <- function(pedigree){
  if (nargs()==0){
    stop(" Pedigree file not specified! ")
  }

  if(!is.data.frame(pedigree)){
    stop("Pedigree must be a data.frame object!")
  }

  nrow = length(pedigree[,1])
  s = pedigree[,2]
  d = pedigree[,3]

  genFinal = rep.int(0, nrow)
  calcGen(s, d, nrow, genVec=genFinal)
}


#' Orders a pedigree (parents before offspring) using Kahn's Algorithm.
#'
#' @param pedigree data.frame with three columns: id, sire and dam
#' @details Kahn's algorithm is a topological sorting method commonly applied to networks. It's only feasible for acyclic directed graphs, so any cycles (ex. an individual that has itself as a parent) will return an error and must be corrected.
#' @references Kahn, Arthur B. (1962), "Topological sorting of large networks", Communications of the ACM
#'
#' @return Integer vector
#' @examples
#' id = c("A","B","C","D","E","F","G","H","I","J","K")
#' sire = c("0","0","0","0","A","A","C","A","G","G","A")
#' dam = c("0","0","0","0","B","H","D","D","B","E","H")
#' ped = as.data.frame(cbind(id,sire,dam))
#' newOrd = getOrdPed(ped)
#' pedReord = ped[ped$id[c(newOrd)],]
#'
#' @export
getOrdPed <- function(pedigree, verbose = FALSE){

  if (missing(pedigree))
    stop("  Pedigree file not specified!  ")

  if (!is.data.frame(pedigree))
    stop (" Pedigree file must be a data.frame object! ")

  id = pedigree[,1]
  s = pedigree[,2]
  d = pedigree[,3]

  nrow=length(id)
  in_degree= rep.int(0, nrow)

  for (i in 1:nrow) {
    if (pedigree[i,2] != 0) {
      in_degree[i] = in_degree[i] + 1
    if (pedigree[i,3] != 0) {
      in_degree[i] = in_degree[i] + 1
      }
    }
  }
  if (verbose == TRUE) cat(in_degree, "\n")
  Q = c()
  for (i in 1:nrow){
    if (in_degree[i] == 0){
      Q = append(Q, pedigree[i,1], 0)
    }
  }

  order_list = c()
  while (length(Q)>0){
    u = Q[length(Q)]
    if (verbose == TRUE) cat("Q= ", Q, "u= ", u, "\n")
    Q = Q[-length(Q)]
    order_list = append(order_list, u, length(order_list))

    progeny_s = which(s %in% id[u])
    progeny_d = which(d %in% id[u])

    if (verbose == TRUE) cat("pr_s", progeny_s, "\n")
    if (verbose == TRUE) cat("pr_d", progeny_d, "\n")

    if (id[u] %in% s){
      for (i in progeny_s){
        in_degree[i] = in_degree[i] - 1
        if (in_degree[i] == 0){
          Q = append(Q, id[i], 0)
    }
    }
    }
    if (id[u] %in% d){
      for (i in progeny_d){
        in_degree[i] = in_degree[i] - 1
        if (in_degree[i] == 0){
          Q = append(Q, id[i], 0)
    }
    }
    }
  }

  if (length(order_list) == length(id)){
    return(order_list)
  }
  else {
    cat("At least one cycle detected \n")
    cat("Individuals involved in cycles: ", Q)
    stop("Correct errors and try again.")
  }}


