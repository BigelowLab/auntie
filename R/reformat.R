#' Retrieve IDs for FastQ pairs 
#'
#' @export
#' @param x list of filepairs (forward and reverse elements) 
#' @return two element list of IDs
#' @examples
#' \dontrun{
#' library(dplyr)
#'  filepairs = dadautils::list_filepairs("/mnt/storage/data/edna/mothur/projects/subsample/sub")
#'  z = dplyr::as_tibble(fastq_ids(filepairs)[[1]])
#'  head(z)
#'  # # A tibble: 6 × 2
#'  #   forward                                                       reverse                                                      
#'  #   <chr>                                                         <chr>                                                        
#'  # 1 A01346:18:H7KY5DRXY:2:2101:25229:2002 1:N:0:TGGTTCGA+AACCGTGT A01346:18:H7KY5DRXY:2:2101:25229:2002 2:N:0:TGGTTCGA+AACCGTGT
#'  # 2 A01346:18:H7KY5DRXY:2:2101:10673:3004 1:N:0:TGGTTCGA+AACCGTGT A01346:18:H7KY5DRXY:2:2101:10673:3004 2:N:0:TGGTTCGA+AACCGTGT
#'  # 3 A01346:18:H7KY5DRXY:2:2101:2211:3223 1:N:0:TGGTTCGA+AACCGTGT  A01346:18:H7KY5DRXY:2:2101:2211:3223 2:N:0:TGGTTCGA+AACCGTGT 
#'  # 4 A01346:18:H7KY5DRXY:2:2101:7907:3724 1:N:0:TGGTTCAA+AACCGTGT  A01346:18:H7KY5DRXY:2:2101:7907:3724 2:N:0:TGGTTCAA+AACCGTGT 
#'  # 5 A01346:18:H7KY5DRXY:2:2101:15248:4820 1:N:0:TGGTTCGA+AACCGTGT A01346:18:H7KY5DRXY:2:2101:15248:4820 2:N:0:TGGTTCGA+AACCGTGT
#'  # 6 A01346:18:H7KY5DRXY:2:2101:28004:5368 1:N:0:TGGTTCGA+AACCGTGT A01346:18:H7KY5DRXY:2:2101:28004:5368 2:N:0:TGGTTCGA+AACCGTGT
#'  tail(z)
#'  # # A tibble: 6 × 2
#'  #   forward                                                        reverse                                                       
#'  #   <chr>                                                          <chr>                                                         
#'  # 1 A01346:18:H7KY5DRXY:2:2278:20826:35305 1:N:0:TGGTTCGA+TACCGTGT A01346:18:H7KY5DRXY:2:2278:20826:35305 2:N:0:TGGTTCGA+TACCGTGT
#'  # 2 A01346:18:H7KY5DRXY:2:2278:5222:35869 1:N:0:TGGTTCGA+TACCGTGT  A01346:18:H7KY5DRXY:2:2278:5222:35869 2:N:0:TGGTTCGA+TACCGTGT 
#'  # 3 A01346:18:H7KY5DRXY:2:2278:12780:36151 1:N:0:TGGTTCGA+TACCGTGT A01346:18:H7KY5DRXY:2:2278:12780:36151 2:N:0:TGGTTCGA+TACCGTGT
#'  # 4 A01346:18:H7KY5DRXY:2:2278:26964:36385 1:N:0:TGGTTCGA+TACCGTGT A01346:18:H7KY5DRXY:2:2278:26964:36385 2:N:0:TGGTTCGA+TACCGTGT
#'  # 5 A01346:18:H7KY5DRXY:2:2278:2103:36417 1:N:0:TGGTTCGA+TACCGTGT  A01346:18:H7KY5DRXY:2:2278:2103:36417 2:N:0:TGGTTCGA+TACCGTGT 
#'  # 6 A01346:18:H7KY5DRXY:2:2278:8956:36667 1:N:0:TGGTTCGA+TACCGTGT  A01346:18:H7KY5DRXY:2:2278:8956:36667 2:N:0:TGGTTCGA+TACCGTGT 
#' }
fastq_ids <- function(x){
  
  nn <- lengths(x)
  if (nn[1] != nn[2]) stop("x must have the same number of files listed in each element")
    
  get_ids <- function(filename){
    x <- ShortRead::id(ShortRead::readFastq(filename))
    as.character(x)
  }
    
  ids <- lapply(seq_along(x[[1]]),
    function(i){
      list(forward = get_ids(x[[1]][i]),
           reverse = get_ids(x[[2]][i]))
    })
  ids
}

#' Retrieve a path to the bbmap reformat.sh executable
#'
#' @export
#' @param default character, the default to use if \code{Sys.which()} doesn't find the application
#' @return character, path to the executable
which_reformat.sh <- function(default = '/mnt/scgc_nfs/opt/common/bbmap/35.85/reformat.sh'){
  app <- Sys.which("reformat.sh")
  if (nchar(app) == 0) {
    java <- Sys.which("java")
    if (nchar(java) == 0) stop("java isn't installed - did you load the 'bbmap' module")
    app <- default[1]
  }
  app
}


#' A wrapper around bbmap's \href{https://github.com/BioInfoTools/BBMap/blob/a9ceda047a7c918dc090de0fdbf6f924292d4a1f/sh/reformat.sh}{reformat.sh} to subsample.
#'
#' @export
#' @param x list or character vector.  If a list then 2 elements which are forward and reverse filenames
#' @param reformat_args character, default for subsampling, see \href{https://github.com/BioInfoTools/BBMap/blob/a9ceda047a7c918dc090de0fdbf6f924292d4a1f/sh/reformat.sh}{reformat.sh}
#' @param path character, the default path to write the subsampled files to
#' @param app character, the path tot he reformat.sh executable
#' @param stdout see \code{system2} - but really for future use if we can figure out how to capture java output - ignore for now
#' @return numeric, 0 for success for all inputs
#' @examples
#' \dontrun{
#' $ module use /mod/scgc
#' $ module load anaconda3/2019.07
#' $ module load bbmap
#' $ module load dada2
#' $ R
#' > filepairs = auntie::list_filepairs("/mnt/storage/data/edna/mothur/projects/subsample/raw")
#' > sub_path <- "/mnt/storage/data/edna/mothur/projects/subsample/sub"
#' > charlier::make_path("/mnt/storage/data/edna/mothur/projects/subsample/sub")
#' > ok <- auntie::bbmap_subsample(filepairs, path = sub_path, reformat_args = "ow=t")
#' }
bbmap_subsample <- function(x,
  path = ".",
  app = which_reformat.sh(),
  reformat_args = "ow=f samplereadstarget=15000 sampleseed=7",
  overwrite = FALSE,
  stdout = ""){
  
  if (is.list(x)){
    # here we handle filepairs
    x <- verify_filepairs(x, action = "stop")
        
    # generate output paths
    outfiles <- sapply(x,
      function(x, path = "."){
        file.path(path[1], basename(x))
      }, path = path, simplify = FALSE)
    
    
    # if no overwrite then compare in/out pairs
    if (!overwrite){
      same <- x[[1]] == outfiles[[1]]
      if (any(same)){
        msg <- sprintf("%i output forward file(s) will overwrite inputs", sum(same))
        warning(msg)
        return(1)
      }
      same <- x[[2]] == outfiles[[2]]
      if (any(same)){
        msg <- sprintf("%i output reverse file(s) will overwrite inputs", sum(same))
        warning(msg)
        return(1)
      }
    }
    # now construct the comand and run  
    ok <- sapply(seq_along(x[[1]]),
      function(i){
        ok <- charlier::make_path(dirname(outfiles[[1]][i]))
        args <- sprintf("in=%s in2=%s out=%s out2=%s %s",
          x[[1]][i], 
          x[[2]][i],
          outfiles[[1]][i],
          outfiles[[2]][i],
          reformat_args)
        
        cat(sprintf("CMD: %s %s", app[1], args), "\n")
        ok <- system2(app[1], args = args, stdout = stdout)
        if (ok != 0) warning("reformat failed for:", x[[1]][i]) 
        ok
      })
  } else {
    #here we have a single vector of files
    outfiles = file.path(path[1], basename(x))
    
    # check for overwrites
    if (!overwrite){
      same <- x == outfiles
      if (any(same)){
        msg <- sprintf("%i output file(s) will overwrite inputs", sum(same))
        warning(msg)
        return(1)
      }
    }
    
    # build commands and run
    ok <- sapply(seq_along(x),
      function(i){
        ok <- charlier::make_path(dirname(outfiles[i]))
        args <- sprintf("in=%s out=%s %s",
          x[i], 
          outfiles[i],
          reformat_args)
          
        cat(sprintf("CMD: %s %s", app[1], args), "\n")
        ok <- system2(app[1], args = args, stdout = stdout)
        if (ok != 0) warning("reformat failed for:", x[i]) 
        ok
      })
    
  }

  # return 0 for all OK and != 0 for all others
  sum(ok)
}