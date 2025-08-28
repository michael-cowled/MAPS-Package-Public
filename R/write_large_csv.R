#' Write a Large Data Frame to CSV in Chunks
#'
#' Writes a data frame to a CSV file chunk-by-chunk to handle large datasets
#' that may exceed memory limits when written all at once.
#'
#' @param df A data frame to write.
#' @param path A string specifying the output file path.
#' @param chunk_size An integer for the number of rows to write per chunk.
#'
#' @return `TRUE` invisibly on success.
#' @export
#' @importFrom utils write.table
#'
#' @examples
#' \dontrun{
#'   # Create a large dummy data frame
#'   large_df <- data.frame(a = 1:100000, b = rnorm(100000))
#'   # Write it to disk in chunks of 10000 rows
#'   write_large_csv(large_df, "large_data.csv", chunk_size = 10000)
#' }
write_large_csv <- function(df, path, chunk_size = 50000) {
  total_rows <- nrow(df)
  if (total_rows == 0) {
    write.csv(df, path, row.names = FALSE)
    return(invisible(TRUE))
  }
  num_chunks <- ceiling(total_rows / chunk_size)
  write.table(df[0, , drop = FALSE], file = path, sep = ",", col.names = TRUE,
              row.names = FALSE, append = FALSE, qmethod = "double")
  for (i in seq_len(num_chunks)) {
    start_row <- (i - 1) * chunk_size + 1
    end_row <- min(i * chunk_size, total_rows)
    chunk <- df[start_row:end_row, , drop = FALSE]
    write.table(chunk, file = path, sep = ",", col.names = FALSE,
                row.names = FALSE, append = TRUE, qmethod = "double")
  }
  invisible(TRUE)
}
