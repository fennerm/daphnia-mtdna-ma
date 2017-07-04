
#' @export
replace_rows <- function(df, old_index, new_row) {
    # exclude indel rows from table
    new_df <- df[-old_index,]

    pre_half <- df[1:(old_index[1] - 1),]
    post_half <- df[(tail(old_index, n = 1) + 1):nrow(df),]
    # replace with the new row

    new_df <- rbind(pre_half, new_row, post_half)
    row.names(new_df) <- NULL
    new_df
}

#' @export
ulapply <- function(...) {
    unlist(lapply(...))
}
