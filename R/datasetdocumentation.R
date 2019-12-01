#' New Orleans data
#'
#' Victims related to human trafficking in Greater New Orleans
#'
#' These data are collected into 8 lists.  For reasons of confidentiality the lists are only labelled as A, B, ..., H.
#'  Full details are given in Bales, Murphy and Silverman (2019).
#'
#' @references K. Bales, L. Murphy and B. W. Silverman (2019). How many trafficked people are there in Greater New Orleans?
#' \emph{Journal of Human Trafficking, to appear}.
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
#' This is a simple data set based on three lists, which gives examples of models that
#'  fail on one or the other of the criteria tested by \code{\link{checkident}}. This is Table 2 in Chan, Silverman and Vincent (2019).
#'
#' If all three two-list effects are included in the fitted model then the linear program
#' in \code{\link{checkident}} yields a strictly positive value but the matrix A is not of full column rank, so the parameters are not identifiable.
#' If the model contains AB either alone or in conjunction with one of AC and BC, then the linear program result is zero, so the MLE does not exist.
#' If only main effects are considered, or if either or both of AC and BC, but not AB are included,
#' then the model passes both tests.
#'
#' @references
#' Chan, L., Silverman, B. W., and Vincent, K. (2019).
#'  Multiple Systems Estimation for Sparse Capture Data: Inferential Challenges when there are Non-Overlapping Lists. Available from \url{https://arxiv.org/abs/1902.05156}.
#'
"Artificial_3"
#'Victims related to sex trafficking in a U.S. Western site
#'
#'These data are collected into 5 lists. For reasons of confidentiality the lists are only labelled as A, B, C, D and E. Full details are
#'  given in Farrell, Dank, Kfafian, Lockwood, Pfeffer, Hughes and Vincent (2019).
#'
#' @references Farrell, A., Dank, M., Kfafian, M., Lockwood, S., Pfeffer, R., Hughes, A., and Vincent, K. (2019).
#' Capturing human trafficking victimization through crime reporting. Technical Report 2015-VF-GX-0105, National Institute of Justice. Available from \url{https://www.ncjrs.gov/pdffiles1/nij/grants/252520.pdf}.
#'
"Western"
#' The Netherlands data
#'
#' Victims related to human trafficking in the Netherlands
#'
#'These data are collected into six lists. Full details are given in Table 2. of Silverman (2019).
#'
#'@references Silverman, B. W. (2019). Model fitting in Multiple Systems Analysis for the quantification of Modern Slavery: Classical and Bayesian approaches
#'\emph{Journal of Royal Statistical Society: Series A to appear}.
#'
"Ned"
#' Netherlands data five list version
#'
#' Netherlands data consolidated into five lists
#'
#' This reduces the Netherlands data \code{\link{Ned}} into five lists, constructed by combining
#'  the two smallest lists I and O into a single list.
#'
"Ned_5"
#' UK data
#'
#'
#'Data from the UK 2013 strategic assessment
#'
#' This is a table of six lists used in the research \href{https://tinyurl.com/ydegfjaw}{published by the Home Office} as part of the strategy leading
#'  to the Modern Slavery Act 2015. The data are considered in six lists, labelled as follows:  LA--Local authorities; NG--Non-government organisations;
#' PF--Police forces; GO--Government organisations; GP--General public; NCA--National Crime Agency.  Each of the first six columns in the data frame
#'  corresponds to one of these lists.
#'   Each of the rows of the data frame corresponds to a possible combination of lists, with value 1 in the relevant
#'  column if the list is in that particular combination.
#'  The last column of the data frame states the number of cases observed in that particular combination of lists.
#'  Combinations of lists for which zero cases are observed are omitted.
#'
#' @references \url{https://www.gov.uk/government/publications/modern-slavery-an-application-of-multiple-systems-estimation}
#'
"UKdat"
#' UK data five list version
#'
#' UK data consolidated into five lists
#'
#' This reduces the UK data \code{\link{UKdat}} into five lists, constructed by combining
#'  the PF and NCA lists into a single PFNCA list
#'
#'
"UKdat_5"

