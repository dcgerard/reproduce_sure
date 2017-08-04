## Miscilaneous functions that are used in big_sims.R and changing_core.R

#' Generates a scaled all-orthonormal array of dimension p with Frobenius norm equal to fnorm_mean.
#'
#' @param p A vector positive integers. The dimension of the array.
#' @param r A vector of positive integers. The multilinear rank of the array.
#' @param fnorm_mean The Frobenius norm of the tensor.
#'
#' @author David Gerard
#'
#' @examples
#' p <- c(5, 7, 11)
#' r <- c(2, 3, 5)
#' fnorm_mean <- 1
#' fout <- generate_core_same(p, r, fnorm_mean)
#' tensr::fnorm(fout)
generate_core_same <- function(p, r, fnorm_mean)
{
  stopifnot(length(p) == length(r))
  stopifnot(all(p >= r))
  x_temp <- array(stats::rnorm(prod(r)), dim = r)
  x_holq <- tensr::holq(x_temp)
  x_isvd <- tensr::get_isvd(x_holq)
  s_temp <- array(0, dim = p)
  indices_s_temp <- as.matrix(expand.grid(lapply(r, seq, from = 1)))
  s_temp[indices_s_temp] <- x_isvd$V
  U_mult <- tensr::hosvd(array(stats::rnorm(prod(p)), dim = p))$U
  theta <- tensr::atrans(s_temp / fnorm(s_temp) * fnorm_mean, U_mult)
  return(theta)
}


#' A rehashing of \code{hose::sure_given_c} but wiht slightly different inputs.
#'
#' NOTE: This only works when tau2 = 1.
#'
#' @param c_x This is the same as \code{obj} in \code{hose::sure_given_c}.
#' @param func A list of functions.
#' @param dfunc A list of derivatives of \code{func}.
#' @param lambda This is the biggest difference between this function and \code{hose::sure_given_c}.
#'     Instead of being a list, this is a matrix, where the first element is hte lambda of the first
#'     mode and the second element is the scaling multiplier.
#'
#' @author David Gerard
sure_optim_wrapper <- function(c_x, func, dfunc, lambda_vec) {
  ## assumes first arg in func and dfunc have two parameters and rest have one
  n <- length(func)
  lambda <- list()
  lambda[[1]] <- lambda_vec[1:2]
  for(mode_index in 2:n)
  {
    lambda[[mode_index]] <- lambda_vec[mode_index + 1]
  }
  hose::sure_given_c(obj = c_x, func = func, dfunc = dfunc, lambda = lambda)$sure_val
}


#' Creates the appropriate parameter list for use in \code{hose::sure} from a vector of thresholding
#' parameters and a scaling parameter.
#'
#' @param lambda A vector of the mode specific soft-thresholding parameters.
#' @param c_current The scaling parameter.
#'
#' @author David Gerard
lambda_list_wrapper <- function(lambda, c_current) {
  lambda_list <- list()
  lambda_list[[1]] <- c(lambda[1], c_current)
  for (mode_index in 2:length(lambda)) {
    lambda_list[[mode_index]] <- lambda[mode_index]
  }
  return(lambda_list)
}

#' Returns a vector whose elements are the positive parts of the elements of x.
#'
#' @param x A numeric.
#'
#' @author David Gerard
pos_part <- function (x) {
  x[x < 0] <- 0
  return(x)
}
