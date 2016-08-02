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
checkPed <- function (pedigree, id=1, sire=2, dam=3) {
  options(warn=-1)

  if (missing(pedigree))
    stop("  Pedigree file not specified!  ")

  if (!is.data.frame(pedigree))
    stop (" Pedigree file must be a data.frame object! ")
  #if (missing(ids))
  #stop("  Need to specify the file with the interested id animals   ")

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
    #stop (" Duplicated ids detected")
  }

  sireclean = sires[sires!="0"]
  damclean = dams[dams!="0"]
  wrongsex = sireclean %in% damclean
  #wrongid = anim %in% sireclean || anim %in% damclean

  if(TRUE %in% wrongsex) {
    stop (" Sire = Dam  ")
  }
  #if(wrongid == "TRUE")
  #	stop (" Anim = Sire or Anim = Dam")

  cat(paste('Total number of animals = ',length(anim),sep=''),'\n')
  #cat('  > Animal ids: ', anim, '\n')
  cat(paste('Total number of unique sires = ',length(sireclean),sep=''),'\n')
  #cat('  > Sire ids: ', sireclean, '\n')
  cat(paste('Total number of unique dams = ',length(damclean),sep=''),'\n')
  #cat('  > Dam ids: ', damclean, '\n')

  cat('Pedigree is OK!')
  # return(pedigree)
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
#' @return Character vector
#' @examples
#' id = c(1,2,3,4,5,6,7,8,9,10)
#' sire = c(0,0,1,1,3,5,5,5,7,7)
#' dam = c(0,0,2,2,2,4,4,6,6,8)
#' ped = as.data.frame(cbind(id,sire,dam))
#' newOrd = getOrdPed(ped)
#'
#' @export
getOrdPed <- function(pedigree){

  id = pedigree[,1]
  s = pedigree[,2]
  d = pedigree[,3]

  nrow=length(id)
  in_degree= rep.int(0, nrow)

  for (i in 1:nrow) {
    if (pedigree[i,2] != 0) {
      in_degree[i] = in_degree[i] + 1
    }
  }
  for (i in 1:nrow) {
    if (pedigree[i,3] != 0) {
      in_degree[i] = in_degree[i] + 1
    }
  }

  Q = c()
  for (i in 1:nrow){
    if (in_degree[i] == 0){
      Q = append(Q,id[i],0)
    }
  }

  order_list = c()
  while (length(Q)!=0){
    u = Q[length(Q)]
    order_list = append(order_list, u, length(order_list))
    Q = Q[-length(Q)]
    for (i in 1:nrow){
      in_degree[i] = in_degree[i] -1
      if (in_degree[i] == 0){
        Q = append(Q, i, 0)
      }
    }
  }
  if (length(order_list) == length(id)){
    return(order_list)
  }
  else {
    cat("At least one cycle detected")
    cat("Individuals involved in cycles: ", Q)
    stop("Correct errors and try again.")
  }
  }

