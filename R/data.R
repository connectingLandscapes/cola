#' Coordinate reference system list
#'
#' A short description of what the dataset contains and where it comes from.
#'
#' @format A data frame with 3 rows and 2 variables:
#' \describe{
#'   \item{auth_name}{Authority name}
#'   \item{code}{CRS code}
#'   \item{name}{CRS short name}
#'   \item{table_name}{Type of CRS}
#'   \item{crs_code}{CRS auth_name:code}
#'   \item{label}{Label to show combining crs_code, name, table_name}
#' }
#' @source \url{https://example.com/where-data-came-from}
"crs_df"

# crs_df$label <- paste(crs_df$crs_code, crs_df$name, gsub('_crs', '', crs_df$table_name) )
# usethis::use_data(crs_df, overwrite = TRUE)
# usethis::use_data(per, overwrite = TRUE)
# usethis::use_data(uper, overwrite = TRUE)

#' Performance times 
#'
#' Time, RAM and CPU consumption from CoLa and UNICOR when running kernels and corridors
#'
#' @format A data frame with 504 rows and 10 variables:
#' \describe{
#'   \item{soft}{Software}
#'   \item{scen}{Scenario name}
#'   \item{size}{Size label, 1-1 match to pixel size}
#'   \item{spix}{Pixels size}
#'   \item{npix}{Number of pixels}
#'   \item{npts}{Number of points}
#'   \item{func}{Function tested}
#'   \item{variable}{Variable measured}
#'   \item{var}{Pretty variable}
#'   \item{value}{Measured value}
#' }
#' @source \url{https://example.com/where-data-came-from}
"per"


#' Performance scenarios details 
#'
#' Scenario name, number and pixel sizes, number of points 
#'
#' @format A data frame with 42 rows and 5 variables:
#' \describe{
#'   \item{scen}{Scenario name}
#'   \item{npix}{Number of points}
#'   \item{spix}{Size pixels}
#'   \item{size}{Size label, 1-1 match to pixel size}
#'   \item{npts}{Number of points}
#' }
"uper"