#### My packages

if(R.version$major == 3){
  library(readr)#, lib.loc = "/gpfs/commons/home/jeinson/R/x86_64-redhat-linux-gnu-library/3.4")
  library(magrittr)#, lib.loc = "/gpfs/commons/home/jeinson/R/x86_64-redhat-linux-gnu-library/3.4")
  library(dplyr)#, lib.loc = "/gpfs/commons/home/jeinson/R/x86_64-redhat-linux-gnu-library/3.4")
  library(tibble)#, lib.loc = "/gpfs/commons/home/jeinson/R/x86_64-redhat-linux-gnu-library/3.4")
  library(stringr)#, lib.loc = "/gpfs/commons/home/jeinson/R/x86_64-redhat-linux-gnu-library/3.4")
  library(forcats)#, lib.loc = "/gpfs/commons/home/jeinson/R/x86_64-redhat-linux-gnu-library/3.4")
  library(tidyr)#, lib.loc = "/gpfs/commons/home/jeinson/R/x86_64-redhat-linux-gnu-library/3.4")
  library(purrr)#, lib.loc = "/gpfs/commons/home/jeinson/R/x86_64-redhat-linux-gnu-library/3.4")
} else if(R.version$major == 4){
  library(readr)
  library(magrittr)
  library(dplyr)
  library(tibble)
  library(stringr)
  library(forcats)
  library(tidyr)
  library(purrr)
}

# Custom defined functions
stfu <- function(x) suppressMessages(suppressWarnings(x))
sniffer <- function(x) x[1:10,1:10]

flag_outliers <- function(x){
  quants <- quantile(x)[c(2,4)]
  iqr <- IQR(x)
  return((x > quants[2] + 1.5 * iqr) | (x < quants[1] - 1.5 * iqr))
}

remove_trailing_digit <- function(x) str_remove(x, "\\..+")

exon_to_coord <- function(x) {
  y <- unlist(str_split(x, "_"))
  paste0(y[1], ":", y[2], "-", y[3])
}

snp_id_to_tabix <- function(snpid){
  x <- unlist(str_split(snpid, "_"))
  paste0(x[1], ":", x[2], "-", x[2])
}

gtex_v8_figure_theme <- function() {
	  return(theme(plot.title = element_text(face="plain", size = 8), 
	               text = element_text(size = 8),
	               axis.text=element_text(size = 7),
	               panel.grid.major = element_blank(),
	               panel.grid.minor = element_blank(),
	               panel.background = element_rect(fill = "transparent"),
	               # legend.box.background = element_rect(fill = "transparent"),
	               axis.line = element_line(colour = "black"),
	               legend.text = element_text(size = 7),
	               legend.title = element_text(size = 8)
	  )
	  )
}

# Thanks Paul for figuring out how to do this. 
collapse <- function(x) {
  # Perform unions if there's an intersection
  y <- lapply(
    X = seq_along(along.with = x),
    FUN = function(i) {
      return(Reduce(
        f = function(a, b) {
          if (length(x = intersect(x = a, y = b))) {
            return(union(x = a, y = b))
          }
          return(a)
        }, 
        x = x[seq.int(from = i, to = length(x = x))]
      ))
    }
  )
  # Remove values that had been added previously
  ni <- vector()
  for (i in seq_along(along.with = y)) {
    if (i == 1L) {
      next
    }
    for (j in seq.int(from = 1L, to = i - 1L)) {
      if (all(x[[i]] %in% x[[j]])) {
        ni <- c(ni, i)
        break
      }
    }
  }
  if (length(x = ni)) {
    y <- y[-ni]
  }
  # Run collapse again
  if (!identical(x = y, y = x)) {
    return(collapse(x = y))
  }
  return(y)
}

# Inverse normal transformation function
inv_norm <- function(x) qnorm((rank(x, na.last='keep') - 0.5)/sum(!is.na(x)))

# A function that gets the names elements that are TRUE (or FALSE)
true_names <- function(x, invert = F) {
  if(!invert){
    names(x[x])
  } else {
    names(x[!x])
  }
}
