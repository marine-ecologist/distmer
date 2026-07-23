pri_to_df <- function(path, long = FALSE) {
  if (length(path) != 1L || !file.exists(path)) {
    stop("`path` must name one existing .pri file.", call. = FALSE)
  }
  if (length(long) != 1L || is.na(long) || !is.logical(long)) {
    stop("`long` must be TRUE or FALSE.", call. = FALSE)
  }

  bytes <- readBin(path, what = "raw", n = file.info(path)$size)
  n_bytes <- length(bytes)

  u8 <- function(i) as.integer(bytes[i])
  u32 <- function(i) {
    if (i < 1L || i + 3L > n_bytes) return(NA_real_)
    z <- as.integer(bytes[i:(i + 3L)])
    sum(z * 256^(0:3))
  }
  text_at <- function(i, n) {
    if (n < 0L || i < 1L || i + n - 1L > n_bytes) return(NA_character_)
    tryCatch(suppressWarnings(rawToChar(bytes[i:(i + n - 1L)])),
             error = function(e) NA_character_)
  }
  lp_string <- function(i) {
    if (i > n_bytes) return(NULL)
    len <- u8(i)
    if (len > 127L || i + len > n_bytes) return(NULL)
    value <- if (len == 0L) "" else text_at(i + 1L, len)
    non_printing <- suppressWarnings(grepl("[^[:print:]\\t ]", value))
    if (is.na(value) || is.na(non_printing) || non_printing) return(NULL)
    list(value = value, next_i = i + len + 1L)
  }
  string_sequence <- function(i, number) {
    out <- character(number)
    for (j in seq_len(number)) {
      item <- lp_string(i)
      if (is.null(item)) return(NULL)
      out[j] <- item$value
      i <- item$next_i
    }
    list(values = out, next_i = i)
  }
  doubles_at <- function(i, number) {
    last <- i + number * 8L - 1L
    if (i < 1L || last > n_bytes) return(NULL)
    con <- rawConnection(bytes[i:last], open = "rb")
    on.exit(close(con))
    readBin(con, what = "double", n = number, size = 8L, endian = "little")
  }

  if (n_bytes < 64L || u32(1L) != 6 || text_at(5L, 6L) != "PRIMER") {
    stop("Not a recognised classic PRIMER .pri file.", call. = FALSE)
  }

  # Worksheet header: numbers of variables and samples, followed by a repeated
  # variable count and the length-prefixed variable names.
  header_candidates <- integer()
  for (i in seq.int(18L, max(18L, n_bytes - 24L))) {
    p0 <- u32(i)
    n0 <- u32(i + 4L)
    if (!is.na(p0) && p0 >= 1 && p0 <= 10000 &&
        n0 >= 1 && n0 <= 1e7 && u32(i + 17L) == p0) {
      candidate <- string_sequence(i + 21L, p0)
      if (!is.null(candidate) && all(nzchar(candidate$values))) {
        header_candidates <- c(header_candidates, i)
      }
    }
  }
  if (length(header_candidates) != 1L) {
    stop("Could not uniquely identify the .pri worksheet header.", call. = FALSE)
  }

  header <- header_candidates[[1L]]
  p <- as.integer(u32(header))
  n <- as.integer(u32(header + 4L))
  variable_info <- string_sequence(header + 21L, p)
  variable_names <- make.unique(variable_info$values)
  after_variables <- variable_info$next_i

  # Find the sample-name sequence: a uint32 count followed by n strings.
  sample_hits <- list()
  upper <- min(n_bytes - 4L, after_variables + 4096L + n * 64L)
  for (i in seq.int(after_variables, upper)) {
    if (u32(i) == n) {
      z <- string_sequence(i + 4L, n)
      if (!is.null(z) && all(nzchar(z$values))) {
        sample_hits[[length(sample_hits) + 1L]] <- z
      }
    }
  }
  if (!length(sample_hits)) {
    sample_names <- as.character(seq_len(n))
    after_samples <- after_variables
  } else {
    sample_names <- sample_hits[[1L]]$values
    after_samples <- sample_hits[[1L]]$next_i
  }

  # Factor blocks can describe either rows (n values) or matrix columns
  # (p values). Scan for both instead of trusting the legacy header count,
  # which does not consistently describe all stored factor worksheets.
  find_factor_blocks <- function(number) {
    values <- list()
    ends <- integer()
    if (n_bytes < 6L) return(list(values = values, ends = ends))
    for (i in seq.int(1L, n_bytes - 5L)) {
      nm <- lp_string(i)
      if (is.null(nm) || !nzchar(nm$value) || u32(nm$next_i) != number) next
      z <- string_sequence(nm$next_i + 4L, number)
      if (is.null(z) || any(!nzchar(z$values))) next
      if (is.null(values[[nm$value]])) {
        values[[nm$value]] <- z$values
        ends <- c(ends, z$next_i)
      }
    }
    list(values = values, ends = ends)
  }

  row_blocks <- find_factor_blocks(n)
  column_blocks <- find_factor_blocks(p)

  # When the matrix is square, avoid representing the same block twice.
  if (n == p && length(column_blocks$values)) {
    row_blocks$values[names(column_blocks$values)] <- NULL
  }

  factor_ends <- c(row_blocks$ends, column_blocks$ends)
  count <- n * p
  search_from <- max(c(after_samples, factor_ends), na.rm = TRUE)
  search_to <- n_bytes - count * 8L + 1L
  if (search_to < search_from) stop("The .pri file is truncated.", call. = FALSE)

  best <- list(score = -Inf, values = NULL, at = NA_integer_)
  for (i in seq.int(search_from, search_to)) {
    x <- doubles_at(i, count)
    if (is.null(x) || any(!is.finite(x))) next
    normal <- x == 0 | (abs(x) >= 1e-100 & abs(x) <= 1e100)
    if (mean(normal) < 0.98 || max(abs(x)) < 1e-12) next
    rounded <- abs(x - round(x, 8L)) <= 1e-10 * pmax(1, abs(x))
    score <- mean(normal) + mean(rounded) + 0.1 * mean(abs(x) > 1e-12)
    if (score > best$score) best <- list(score = score, values = x, at = i)
  }
  if (is.null(best$values) || best$score < 1.75) {
    stop("Could not validate a numerical matrix in this .pri layout.", call. = FALSE)
  }

  data_matrix <- matrix(
    best$values, nrow = n, ncol = p, byrow = TRUE,
    dimnames = list(NULL, variable_names)
  )

  row_metadata <- data.frame(Sample = sample_names, stringsAsFactors = FALSE,
                             check.names = FALSE)
  for (nm in names(row_blocks$values)) {
    row_metadata[[nm]] <- factor(row_blocks$values[[nm]])
  }

  variable_metadata <- data.frame(Site = variable_names,
                                  stringsAsFactors = FALSE,
                                  check.names = FALSE)
  for (nm in names(column_blocks$values)) {
    variable_metadata[[nm]] <- factor(column_blocks$values[[nm]])
  }

  if (long) {
    df <- row_metadata[rep(seq_len(n), each = p), , drop = FALSE]
    df$Site <- rep(variable_names, times = n)
    if (ncol(variable_metadata) > 1L) {
      site_rows <- match(df$Site, variable_metadata$Site)
      df <- cbind(
        df,
        variable_metadata[site_rows, -1L, drop = FALSE]
      )
    }
    df$Abundance <- as.vector(t(data_matrix))
  } else {
    # Community matrix for site-level analyses: sites are rows and species are
    # columns. Site identifiers are row names so every data column is numeric.
    df <- as.data.frame(t(data_matrix), check.names = FALSE)
    names(df) <- make.unique(sample_names)
    #rownames(df) <- variable_names
  }

  if (long) rownames(df) <- NULL
  #rownames(variable_metadata) <- variable_metadata$Site
  #rownames(row_metadata) <- row_metadata$Sample

  result <- list(
    df = df,
    metadata = list(
      sites = variable_metadata,
      species = row_metadata
    )
  )
  attr(result, "pri_file") <- normalizePath(path, mustWork = TRUE)
  attr(result, "pri_matrix_byte_offset") <- best$at - 1L
  result
}
