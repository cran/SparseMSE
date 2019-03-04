#' New Orleans data
#'
#' Victims related to modern slavery and human trafficking in New Orleans
#'
#' These data are collected into 8 lists.  For reasons of confidentiality the labels of the lists are anonymised.
#'  Fuller details are given in Bales, Murphy and Silverman (2018).
#'
#' @references K. Bales, L. Murphy and B. W. Silverman (2018). How many trafficked and enslaved people are there in New Orleans? Available from \url{https://tinyurl.com/ybfb9tg6}.
"NewOrl"
#' New Orleans data five list version
#'
#' New Orleans data consolidated into five lists
#'
#' This reduces the New Orleans data \code{\link{NewOrl}} into five lists, constructed by combining
#'  the four smallest lists B, E, F and G into a single list.
#'
#'
"NewOrl_5"
#' Artificial data set to demonstrate possible instabilities
#'
#' This is a simple data set based on three lists, which shows that there is not necessarily any clear hierarchical
#'  relationship between models that fail on one or the other of the criteria tested by \code{\link{checkident}}.
#'
#' If all three interactions are included in the fitted model then then the linear program
#' in \code{\link{checkident}} yields a strictly positive value but the matrix A is not of full column rank, so the parameters are not identifiable.
#' If the model contains AB either alone or in conjunction with one of AC and BC, then the linear program result is zero, so the MLE does not exist.
#' If only main effects are considered, or if either or both of AC and BC, but not AB are included,
#' then the model passes both tests.
#'
"Artificial_3"
#'Victims related to sex trafficking in a U.S. Western site
#'
#'These data are collected into 5 lists. For reasons of confidentiality the labels of the lists are anonymised. Fuller details are
#'  given in Farrell, Dank, Kfafian, Lockwood, Pfeffer, Hughes, Vincent (2019).
#'
#' @references Farrell, A., Dank, M., Kfafian, M., Lockwood, S., Pfeffer, R., Hughes, A., and Vincent, K. (2019).
#' Capturing human trafficking victimization through crime reporting. Technical Report 2015-VF-GX-0105, National Institute of Justice. Available from \url{https://www.ncjrs.gov/pdffiles1/nij/grants/252520.pdf}.
#'
"Western"
