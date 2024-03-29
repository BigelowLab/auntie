#' Count the length of filepairs
#'
#' @export
#' @param x a list with file pairings as character vectors
#' @return named numeric vector with counts for each element
count_filepairs <- function(x){
  lengths(x)
}

#' Verify that a list of file pairs is suitable for processing.
#'
#' @export
#' @param x a list with file pairings as character vectors
#' @param min_size numeric, files with fewer bytes than this number are considered empty. Set to NA or NULL to ignore file size.
#' @param elements character, the names of the file pair elements to test
#' @param require_reverse logical, if TRUE then insist that reverse files must be present (ala Illumina). If FALSE (the default) then allow reverse files to be absent (ala PacBio)
#' @param action character, either stop or warn upon discovery of mismatched files or all file sizes
#'   fall below the minimum size (if tested) 
#' @return the input list with possibly all elements removed if all files fail the size test
verify_filepairs <- function(x, 
  min_size = NULL, 
  elements = c("forward", "reverse"), 
  require_reverse = FALSE,
  action = c("stop", "warn")[2]){
  
    # a silly convenience function
  stop_or_warn <- function(action, ...){
    switch(tolower(action[1]),
      'stop' = stop(...),
      warning(...))
  }
  
  if (!all((elements %in% names(x)))) {
    stop("input is missing one or more required elements - are they named according to 'element' argument?")
  }
  
  ll <- lengths(x)
  if ((length(ll) < 2) || (ll[2] == 0)){
    if (require_reverse) stop_or_warn("stop", "reverse elements are required")
  } else {
    if (!all(ll %in% ll[1])) stop_or_warn("stop","elements of input must be the same length")
  }
  
  # test each file pair against min_size
  # for any pair with at least one is small, drop it
  # and provide warning
  if (!charlier::is.nullna(min_size[1])){
    
    y <- size_filepairs(x, min_size = min_size[1], verify = FALSE)
    small_ix <- is.na(y$size_pass) | !y$size_pass
    if (all(small_ix)){
      stop_or_warn(action, "All files are smaller than min_size: ", min_size[1])
      # empty the list
      x <- sapply(x, function(x) {x <- character(0); x}, simplify = FALSE)
    } else if (any(small_ix)){
      x <- y %>% 
        dplyr::filter(!small_ix) %>%
        dplyr::select(dplyr::all_of(elements)) %>%
        as.list()
      y <- y %>% 
        dplyr::filter(small_ix) %>%
        `[[`(1)
      msg <- sprintf("%i small filepairs encountered, dropping: %s", 
                     length(y), paste(basename(y), collapse = ", "))
      stop_or_warn("warn", msg)
    }
  }
  x
}


#' Given a filepair list, make a vector of alternating filenames
#'
#' @export
#' @param x a filepairs list (two element list for forard and reverse filenames)
#' @return a vector of interleaved filenames
interleave_filepairs <- function(x){
  nf <- length(x[[1]])
  nr <- length(x[[2]])
  if (nf != nr) stop("input must two element listing of files of equal length")
  n <- nf + nr
  fi <- seq(from = 1, to = n, by = 2)
  ri <- seq(from = 2, to = n, by = 2)
  r <- rep("", n)
  r[fi] <- x[[1]]
  r[ri] <- x[[2]]
  r
}

#' Given a even numbered vector of filenames, deinterleaf into a filepairs list
#'
#' @export
#' @param x a vector of items
#' @param elements character,the names of the elements
#' @return return a filepairs list
deinterleave_filepairs <- function(x, elements = c("forward", "reverse")){
  n <- length(x)
  if (as.logical(n %% 2)) stop("input must have even number of elements")
  fi <- seq(from = 1, to = n, by = 2)
  ri <- seq(from = 2, to = n, by = 2)
  list(x[fi], x[ri]) %>%
    rlang::set_names(elements)
}


#' Compute and test the size of filepairs
#' 
#' @export
#' @param x named list of filepairs, typically 'forward' and 'reverse'
#' @param min_size numeric or NA, the mininum size of a file required
#' @param verify logical, if TRUE pass the input to \code{verify_filepairs}
#' @param ... other arguments passed to \code{verify_filepairs} except \code{min_size}
#' @return tibble of 
#' \itemize{
#'   \item forward, charcater vector of foreward filenames
#'   \item reverse, charcater vector of reverse filenames
#'   \item forward_size, numeric file sizes in bytes
#'   \item reverse_size, numeric file sizes in bytes
#'   \item size_pass logical, TRUE/FALSE if the size passes and NA when min_size is NA
#' }
size_filepairs <- function(x, min_size = NA, verify = FALSE, ...){
   if (verify) x <- verify_filepairs(x) 
   nm <- names(x)
   sz <- sapply(x, function(x) file.info(x)$size, simplify = FALSE)
   nmsz <- paste0(nm, "_size")
   if (length(x) == 2){
     size_pass <- sz[[1]] > min_size[1] & sz[[2]] > min_size[1]
     r <- dplyr::as_tibble(x) %>%
       dplyr::mutate(
         !!nmsz[1] := sz[[1]],
         !!nmsz[2] := sz[[2]],
         size_pass = size_pass)
   } else {
     size_pass <- sz[[1]] > min_size[1] 
     r <- dplyr::as_tibble(x) %>%
       dplyr::mutate(
         !!nmsz[1] := sz[[1]],
         size_pass = size_pass)
   }
   r
}


#' List files and separate into forward and reverse file lists based upon filenaming patterns.
#'
#' If reverse files are not present (ala \code{PacBio}) then that element is set to a zero-length
#' character element.
#'
#' @export
#' @param path character, the input path
#' @param pattern_forward file pattern, the default is to match "^.*R1_001"
#' @param pattern_reverse file pattern, the default is to match either "^.*R2_001"
#' @param patterns_exclude one or more file patterns to exclude, default is "^.*\\.cutadapt_output\\.txt$"
#'        Unlike \code{pattern_forward} and \code{pattern_reverse} this can have multiple elements.  Set to 
#'        NULL of NA to skip.
#' @param glob logical, if \code{TRUE} the input patterns are considered file globs like "*_R1_001.fastq"
#'        and will be converted to regex patterns using \code{\link[utils]{glob2rx}}.
#'        If glob is \code{FALSE} then the the patterns are passed to \code{\link[base]{list.files}}
#'        as provided by the user.
#' @param verify logical, if TRUE test that the filepairs are the same length
#' @param ... further arguments for \code{\link[base]{list.files}}
#' @return named list of sorted forward and reverse filenames (reverse possibly zero-length character)
list_filepairs <- function(path,
                       pattern_forward = "^.*R1_001",
                       pattern_reverse = "^.*R2_001",
                       patterns_exclude = "^.*\\.cutadapt_output\\.txt$",
                       glob = FALSE,
                       verify = TRUE,
                       ...){

  if (glob){
    pattern_forward = utils::glob2rx(pattern_forward)
    pattern_reverse = utils::glob2rx(pattern_reverse)
  }

  x <- list(
    forward = sort(list.files(path, pattern = pattern_forward, full.names = TRUE, ...)),
    reverse = sort(list.files(path, pattern = pattern_reverse, full.names = TRUE, ...)) 
  )
  
  if (!charlier::is.nullna(patterns_exclude[1])){
    if (glob) patterns_exclude = utils::glob2rx(patterns_exclude)
    x <- sapply(x,
      function(x){
        ix <- charlier::mgrepl(patterns_exclude, x, fixed = glob)
        x[!ix]
      }, simplify = FALSE)
  }
  

  if (verify) x <- verify_filepairs(x)

  x
}

