#!/usr/bin/env Rscript

#' bigdawg_object
#'
#' Create big dog object:
#' @param input The object input
#' @note This function is for internal BIGDAWG use only.
bigdawg_object <- function(input){
  setClass("bigdawg",
           slots=list(name="character",
                      id="numeric",
                      contents="character"))
  obj <- new("bigdawg",
             name="bigdawg",
             id=as.integer(runif(1, 1,10000)),
             contents="content")
  is.object(obj)
  isS4(obj)
  obj@name
  obj@id
  #Use generic function show() below
  setMethod("show",
            "bigdawg",
            function(object) {
              cat("Name:",object@name, "\n")
              cat("Id:",object@id, "\n")
              cat("Content:", object@contents, "\n")
            }
  )
  obj
  return(obj)
}
