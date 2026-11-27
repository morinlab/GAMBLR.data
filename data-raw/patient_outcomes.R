# patient_outcomes.R
#
# Build-time helper for the public patient_outcomes table of gambl_mutations.db: one row per patient in the public
# sample_meta whose study may release outcomes, with only the outcome groups that study released. `source()` this
# file (it defines functions; it does not run anything). Used by assemble_bundled_data.R.
#
#   data-raw/public_outcome_columns.txt  every column of get_gambl_outcomes(), included or excluded ("-name"):
#                                         the build stops on a column nobody has decided about
#   data-raw/public_outcome_studies.tsv  which groups (OS, PFS, IPI, ...) each study may release

# Columns of each outcome group
outcome_groups <- list(
    OS = c("CODE_OS", "OS_YEARS", "alive"),
    PFS = c("CODE_PFS", "PFS_YEARS"),
    DSS = c("CODE_DSS", "DSS_YEARS"),
    TTP = c("CODE_TTP", "TTP_YEARS"),
    TTT = c("CODE_TTT", "TTT_YEARS"),
    IPI = "IPI",
    FLIPI = "FLIPI",
    age = c("age", "is_adult"),
    treatment = c("tx_primary_class", "tx_ritu_maintenance")
)

# The include/exclude lists of data-raw/public_outcome_columns.txt: bare names are included, "-name" excluded
read_outcome_column_lists <- function(path = "data-raw/public_outcome_columns.txt") {
    lines <- trimws(readLines(path))
    lines <- lines[nzchar(lines) & !startsWith(lines, "#")]
    list(include = lines[!startsWith(lines, "-")],
         exclude = sub("^-", "", lines[startsWith(lines, "-")]))
}

# data-raw/public_outcome_studies.tsv as a named list: study -> outcome groups
read_outcome_studies <- function(path = "data-raw/public_outcome_studies.tsv") {
    s <- utils::read.delim(path, comment.char = "#", stringsAsFactors = FALSE)
    groups <- lapply(strsplit(s$groups, ","), trimws)
    unknown <- setdiff(unlist(groups), names(outcome_groups))
    if (length(unknown) > 0) {
        stop("data-raw/public_outcome_studies.tsv names unknown outcome group(s): ", paste(unknown, collapse = ", "),
             ". Known: ", paste(names(outcome_groups), collapse = ", "), call. = FALSE)
    }
    stats::setNames(groups, s$study)
}

# A time of exactly 0 with a censored event code is a placeholder for "no follow-up" (e.g. Schmitz patients not in
# the survival analysis): both become NA, endpoint by endpoint
drop_zero_follow_up <- function(outcomes) {
    for (ep in c("OS", "PFS", "DSS", "TTP", "TTT")) {
        time <- paste0(ep, "_YEARS")
        code <- paste0("CODE_", ep)
        if (!all(c(time, code) %in% colnames(outcomes))) next
        none <- outcomes[[time]] %in% 0 & outcomes[[code]] %in% 0
        outcomes[[time]][none] <- NA
        outcomes[[code]][none] <- NA
        if (ep == "OS" && "alive" %in% colnames(outcomes)) outcomes$alive[none] <- NA
    }
    outcomes
}

# outcomes: GAMBLR.results::get_gambl_outcomes() (years, cleaned); sample_meta: the final public sample_meta, with
# patient_id and study. Stops if outcomes has a column the column lists don't mention, if an included or released
# column is missing, or on duplicate patients. Each patient keeps only the groups their studies may release;
# patients with nothing left get no row.
assemble_patient_outcomes <- function(outcomes, sample_meta,
                                      column_lists = read_outcome_column_lists(),
                                      studies = read_outcome_studies()) {
    cols <- setdiff(colnames(outcomes), "patient_id")
    unlisted <- setdiff(cols, c(column_lists$include, column_lists$exclude))
    if (length(unlisted) > 0) {
        stop("get_gambl_outcomes() has column(s) not listed in data-raw/public_outcome_columns.txt: ",
             paste(unlisted, collapse = ", "), ". Add each one there, as included or excluded (\"-name\"), ",
             "before building.", call. = FALSE)
    }
    absent <- setdiff(column_lists$include, cols)
    if (length(absent) > 0) {
        stop("data-raw/public_outcome_columns.txt includes column(s) get_gambl_outcomes() doesn't return: ",
             paste(absent, collapse = ", "), call. = FALSE)
    }
    released_cols <- unique(unlist(outcome_groups[unique(unlist(studies))]))
    not_included <- setdiff(released_cols, column_lists$include)
    if (length(not_included) > 0) {
        stop("data-raw/public_outcome_studies.tsv releases column(s) that public_outcome_columns.txt doesn't ",
             "include: ", paste(not_included, collapse = ", "), call. = FALSE)
    }
    if (any(duplicated(outcomes$patient_id))) {
        stop("get_gambl_outcomes() has more than one row for some patients: ",
             paste(utils::head(unique(outcomes$patient_id[duplicated(outcomes$patient_id)]), 5), collapse = ", "),
             call. = FALSE)
    }
    outcomes <- drop_zero_follow_up(as.data.frame(outcomes))

    # the groups each patient may release: the union over the patient's studies
    patient_studies <- unique(sample_meta[sample_meta$study %in% names(studies), c("patient_id", "study")])
    out <- outcomes[outcomes$patient_id %in% patient_studies$patient_id, c("patient_id", released_cols), drop = FALSE]
    for (i in seq_len(nrow(out))) {
        groups <- unique(unlist(studies[patient_studies$study[patient_studies$patient_id == out$patient_id[i]]]))
        hidden <- setdiff(released_cols, unlist(outcome_groups[groups]))
        out[i, hidden] <- NA
    }
    has_value <- rowSums(!is.na(out[, released_cols, drop = FALSE])) > 0
    out <- out[has_value, , drop = FALSE]
    rownames(out) <- NULL
    out
}
