#' @title customfun objects
#' @name customfun
#' @rdname customfun
#' @description Include a custom function that can be used to manipulate the projection matrix (A).
NULL

#e.g nick's nifty example, I'll need to chat to nick about this.... I'm not quite sure how he's done this.

# ddfun <- function (habitat,params) {
#   adult_density <- population(habitat)$adult / area(habitat)$area
#   param$p * exp(-adult_density / param$area)
# }
# 
# # turn it into a transfun object
# dd <- pop::as.customfun(ddfun,
#                        param = list(p = 0.9,
#                                area = 100),
#                         type = 'probability')

# other possible functions - recruitment, dispersal, stochastic, catastrophies, 