#' @title module objects
#' @name module
#' @rdname module
#' @description Include a module that drives changes in a dynamic landscape.
NULL

#' as.module <- function(){
#'   
#' }
#' 
#' # things related to the transfun class
#' 
#' transfunClasses <- function () {
#'   # list all available classes of transfun
#'   c('probability', 'rate', 'dispersal', 'compound')
#' }
#' 
#' transfunType <- function (x) {
#'   # get the type of the transfun object
#'   stopifnot(is.transfun(x))
#'   classes <- class(x)
#'   matches <- na.omit(match(transfunClasses(), classes))
#'   if (length(matches) == 0) {
#'     stop ('this transfun object does not correspond to any known transfun types')
#'   } else if (length(matches) > 1) {
#'     stop ('this transfun object correspond to multiple transfun types')
#'   } else {
#'     type <- classes[matches]
#'   }
#'   return (type)
#' }
#' 
#' #' @title transfun objects
#' #' @name transfun
#' #' @rdname transfun
#' #' @param x a transfun object to print or an object to test as a transfun object
#' #' @description utility functions for the \code{transfun} class. \code{transfun}
#' #'   objects are created by functions such as \code{\link{probability}}.
#' #' @export
#' #' @examples
#' #' prob <- p(0.3)
#' #' is.transfun(prob)
#' #'
#' is.transfun <- function (x) inherits(x, 'transfun')
#' 
#' #' @rdname transfun
#' #' @param \dots further arguments passed to or from other methods.
#' #' @export
#' #' @examples
#' #' prob
#' print.transfun <- function(x, ...) {
#'   if (containsUserTransfun(x)) {
#'     text <- sprintf('user-specified %s transfun',
#'                     transfunType(x))
#'   } else {
#'     landscape <- as.landscape(NULL)
#'     text <- sprintf('%s transfun with expectation %s\n',
#'                     transfunType(x),
#'                     x(landscape))
#'   }
#'   
#'   cat(text)
#' }
#' 
#' is.compound <- function (x) inherits(x, 'compound')
#' 
#' as.compound <- function (x) {
#'   # define a compound transfun class
#'   if (!is.compound(x)) {
#'     class(x) <- c('compound', 'transfun', class(x))
#'   }
#'   return (x)
#' }
#' 
#' #' @rdname transfun
#' #' @param y a transfun object to be multiplied with another with the same
#' #'   pathway
#' #' @details multiplication of transfun objects with the same pathway results in
#' #'   a compound transfun object (also of class \code{transfun}). When used in a
#' #'   stochastic model, the two stochastic transitions are evaluated one after
#' #'   another. When analysed deterministically, the expectation of the compound
#' #'   transition function is taken as the product of the expectations of the two
#' #'   basis transfuns.
#' #' @export
#' #' @examples
#' #' (compound <- prob * r(4.3))
#' #'
#' `*.transfun` <- function (x, y) {
#'   # given two transfun objects, combine them into a compound transfun
#'   stopifnot(is.transfun(x))
#'   stopifnot(is.transfun(y))
#'   
#'   # make sure a dispersal is not combined with a rate or another dispersal
#'   combineDispersalCheck(x, y)
#'   
#'   # if it contains a dispersal, the other must be only a probability
#'   if (contains(x, 'dispersal') |
#'       contains(y, 'dispersal')) {
#'     
#'     # if its a dispersal and probability
#'     z <- function (landscape) {
#'       probdisp(x, y, landscape)
#'     }
#'     
#'   } else {
#'     
#'     # if none are dispersals, just the product
#'     z <- function (landscape) {
#'       x(landscape) * y(landscape)
#'     }
#'     
#'   }
#'   
#'   # coerce class & return
#'   z <- as.compound(z)
#'   return (z)
#'   
#' }
#' 
#' #' @title create a transition function
#' #' @name as.transfun
#' #' @description A utility function to enable users to create bespoke transition
#' #'   functions (\code{transfun} objects) for use in \code{transition}s.
#' #' @param fun an R function describing the transition. This must take only one
#' #'   argument: \code{landscape} and return a numeric vector (see
#' #'   \code{details}).
#' #' @param param a named list of the parameters of \code{fun} (see
#' #'   \code{details}).
#' #' @param type what type of transition this function represents, a probability
#' #'   or a rate
#' #' @details \code{fun} must take only one argument, \code{landscape}, an object
#' #'   of class \code{\link{landscape}}. \code{landscape} objects contain three
#' #'   elements which may be used in the function: \code{population}, a dataframe
#' #'   giving the number of individuals of each stage (columns) in each patch
#' #'   (rows); \code{area}; a numeric vector giving the area of each patch in
#' #'   square kilometres; and \code{features}, a dataframe containing
#' #'   miscellaneous features (columns) of each habitat patch (rows), such as
#' #'   measures of patch quality or environmental variables. See examples for an
#' #'   illustration of how to these objects. Parameters of the transfun should be
#' #'   passed to \code{as.transfun} as a named list. These can then be used in
#' #'   \code{fun} by accessing them from this list. Note that \code{param} isn't
#' #'   an argument to \code{fun}, instead it's modified directly in the function's
#' #'   envirnment (because \emph{reasons}).
#' #' @export
#' #' @examples
#' #' # a very simple (and unnecessary, see ?p) transfun
#' #' fun <- function(landscape) param$prob
#' #' prob <- as.transfun(fun, param = c(prob = 0.3), type = 'probability')
#' #'
#' #' # a density-dependent probability
#' #' dd_fun <- function (landscape) {
#' #'     adult_density <- population(landscape, 'adult') / area(landscape)
#' #'     param$p * exp(- adult_density/param$range)
#' #' }
#' #'
#' #' dd_prob <- as.transfun(dd_fun,
#' #'                        param = list(p = 0.8,
#' #'                                     range = 10),
#' #'                        type = 'probability')
#' #'
#' as.transfun <- function (fun,
#'                          param,
#'                          type = c('probability', 'rate', 'dispersal')) {
#'   
#'   # line up the transfun type
#'   type <- match.arg(type)
#'   
#'   # check it's a function
#'   stopifnot(is.function(fun))
#'   
#'   # check landscape is the only argument
#'   args <- names(formals(fun))
#'   if (length(args) != 1 || args != 'landscape') {
#'     stop ("transfun objects must only take the argument 'landscape'
#'           see ?as.transfun for details and examples")
#'   }
#'   
#'   # define parameters and fun here so fun can see param
#'   param <- param
#'   environment(fun) <- environment()
#'   
#'   # assign type and return
#'   fun <- switch(type,
#'                 probability = as.probability(fun),
#'                 rate = as.rate(fun),
#'                 dispersal = as.dispersal(fun))
#'   
#'   attr(fun, 'user-defined') <- TRUE
#'   
#'   return (fun)
#'   }
#' 
#' containsUserTransfun <- function (transfun) {
#'   # test whether a transfun object contains a user-defined transfun
#'   
#'   # get transfun type, if it's a compound, call this function recursively
#'   type <- transfunType(transfun)
#'   
#'   if (type == 'compound') {
#'     # expand and test components
#'     tf_x <- environment(transfun)$x
#'     tf_y <- environment(transfun)$y
#'     ans <- containsUserTransfun(tf_x) |
#'       containsUserTransfun(tf_y)
#'   } else {
#'     # otherwise test this
#'     ans <- attr(transfun, 'user-defined')
#'     if (is.null(ans)){
#'       ans <- FALSE
#'     }
#'   }
#'   
#'   return (ans)
#'   
#' }
#' 
#' #' @rdname transfun
#' #' @export
#' #' @examples
#' #' # extract the transfun parameters
#' #' (param_prob <- parameters(prob))
#' #' (param_compound <- parameters(compound))
#' #'
#' parameters.transfun <- function (x) {
#'   if (is.compound(x)) {
#'     tf_x <- environment(x)$x
#'     tf_y <- environment(x)$y
#'     param <- c(parameters(tf_x), parameters(tf_y))
#'   } else {
#'     param <- environment(x)$param
#'   }
#'   return (param)
#' }
#' 
#' #' @rdname transfun
#' #' @export
#' #' @param value a named list of parameters matching those currently defined for \code{x}
#' #' @examples
#' #' # update the parameters of these transfuns
#' #' param_prob$p <- 0.6
#' #' parameters(prob) <- param_prob
#' #' parameters(prob)
#' #'
#' #' param_compound$r <- 15
#' #' parameters(compound) <- param_compound
#' #' parameters(compound)
#' `parameters<-.transfun` <- function (x, value) {
#'   
#'   if (is.compound(x)) {
#'     
#'     # get components
#'     components <- list(environment(x)$x, environment(x)$y)
#'     
#'     # do components in turn
#'     for (i in 1:2) {
#'       new_param <- old_param <- parameters(components[[i]])
#'       
#'       # loop through parameters in this component
#'       for (j in 1:length(new_param)) {
#'         
#'         # get first match
#'         which_value <- which(names(value) == names(old_param)[j])[1]
#'         
#'         # update in new_param
#'         new_param[j] <- value[which_value]
#'         
#'         # remove from value
#'         value <- value[-which_value]
#'         
#'       }
#'       
#'       # update transfun
#'       parameters(components[[i]]) <- new_param
#'       
#'     }
#'     
#'     # recombine the components
#'     x <- components[[1]] * components[[2]]
#'     
#'   } else {
#'     # otherwise update basis transfun
#'     
#'     # check new parameters
#'     parametersCheck(value, x)
#'     
#'     # if that worked, define param and fun here
#'     param <- value
#'     environment(x) <- environment()
#'     
#'   }
#'   
#'   return (x)
#'   
#' }
#' 
#' parametersCheck <- function (param, transfun = NULL) {
#'   
#'   # check incoming parameters make sense
#'   stopifnot(is.list(param))
#'   stopifnot(all(sapply(param, is.finite)))
#'   stopifnot(all(sapply(param, is.numeric)))
#'   
#'   if (!is.null(transfun)) {
#'     # check they match the transfun
#'     old_param <- parameters(transfun)
#'     stopifnot(all.equal(names(param), names(old_param)))
#'     stopifnot(length(param) == length(old_param))
#'   }
#'   
#' }
#' 
#' probdisp <- function (x, y, landscape) {
#'   # get expected dispersal fraction from a probability and a dispersal transfun.
#'   # dispersal should have diagonal giving probability of staying, off-diagonals
#'   # giving probability of moving to each other patch, w/ all rows summing to 1.
#'   # probability is probability of leaving (1-probability of staying).
#'   
#'   # work out which way round
#'   if (is.probability(x) & is.dispersal(y)) {
#'     prob <- x(landscape)
#'     disp <- y(landscape)
#'   } else {
#'     prob <- y(landscape)
#'     disp <- x(landscape)
#'   }
#'   
#'   # multiply each row by the dispersal probability
#'   disp <- sweep(disp, 1, prob, '*')
#'   
#'   # add fraction not attempting dispersal back onto diagonal
#'   diag(disp) <- diag(disp) + 1 - prob
#'   
#'   return (disp)
#'   
#' }
#' 
#' combineDispersalCheck <- function (x, y) {
#'   bad_thing1 <- contains(x, 'dispersal') &
#'     (contains(y, 'dispersal') | contains(y, 'rate'))
#'   bad_thing2 <- contains(y, 'dispersal') &
#'     (contains(x, 'dispersal') | contains(x, 'rate'))
#'   if (bad_thing1 | bad_thing2) {
#'     stop ('dispersal transfuns can only be combined with probability transfuns')
#'   }
#' }