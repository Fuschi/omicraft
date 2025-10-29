#' Filter `omic` or `omics` Objects by Sample Metadata and Abundance Fields
#'
#' Filters the samples in an `omic` or `omics` object based on user-specified 
#' conditions provided through `...`. Each condition can reference columns in the 
#' sample metadata (`meta`) or abundance-related fields (`abun`, `rela`, `norm`). 
#' When multiple conditions are provided, they are combined using a logical AND, 
#' meaning a sample must satisfy **all** conditions to be retained.
#'
#' @param object An `omic` or `omics` object. For `omics`, each contained 
#'   `omic` can be filtered according to the combined conditions, with a column 
#'   `omic` distinguishing them in any joined data frames.
#' @param ... Unquoted expressions for filtering. These expressions can reference 
#'   columns in the sample metadata (`meta`) or any of the abundance-related fields 
#'   (`abun`, `rela`, `norm`). Each expression is applied sequentially, and only 
#'   samples that satisfy **all** expressions are retained.
#' @param .by A character vector of grouping columns. The default differs by object 
#'   type:
#'   \itemize{
#'     \item For `omic`, the default is `.by = "sample_id"`.
#'     \item For `omics`, the default is `.by = c("omic", "sample_id")`.
#'   }
#'   If `.by` is `NULL`, no grouping is applied, and conditions are evaluated 
#'   across all samples globally. If `.by` is missing, any existing meta grouping 
#'   (from `group_meta()`) is used instead.
#'
#' @details
#' This function integrates closely with **tidyverse** conventions. Each expression 
#' in `...` is captured and evaluated using tidy evaluation, allowing you to write 
#' conditions similarly to `dplyr::filter()`. When an expression references 
#' abundance variables (`abun`, `rela`, or `norm`), the function automatically 
#' joins the sample metadata with the relevant abundance data (`long_abun`) before 
#' filtering.
#'
#' **Grouping Logic**:
#' 1. **If you provide an explicit `.by`** (a non-empty character vector), those 
#'    columns define the grouping.
#' 2. **If `.by` is `NULL`**, no grouping is applied.
#' 3. **If `.by` is missing**:
#'    - **If the object has existing meta grouping** (set via `group_meta()`), that 
#'      grouping is used.
#'    - **Else**, if any filter expression references abundance variables, 
#'      grouping defaults to:
#'      \itemize{
#'        \item `c("omic", "sample_id")` for `omics` objects.
#'        \item `"sample_id"` for `omic` objects.
#'      }
#'    - **Otherwise**, no grouping is applied.
#'
#' **Multiple Conditions**:
#' Each condition in `...` is applied sequentially, narrowing down the set of 
#' remaining samples. After processing all conditions, only those samples that 
#' satisfied **all** conditions remain.
#'
#' **Behavior with `omics`**:
#' - The `omic` column in the joined data frames identifies which `omic` each sample 
#'   belongs to.
#' - Filtering expressions can apply to metadata or abundance columns, referencing 
#'   each `omic`'s data.
#' - Samples are filtered **independently** for each `omic` based on the conditions.
#'
#' **Return Value**:
#' - For an `omic` object, an updated `omic` containing only the samples that 
#'   passed all filters.
#' - For a `omics`, an updated list of `omic` objects. Each `omic` is 
#'   subset to the remaining samples.
#'
#' This approach provides a flexible, `dplyr`-like interface for complex filtering 
#' of both sample metadata and abundance data, respecting the specified grouping context.
#'
#' @return An `omic` or `omics` object containing only the samples (and their 
#'   corresponding abundance data, if any) that satisfy all provided filter conditions. 
#'   The rest are removed. If no samples match, an object with zero samples is returned 
#'   (but the same structure otherwise).
#'
#' @export
#' @aliases filter_meta,omic-method filter_meta,omics-method
setGeneric("filter_meta", function(object, ...) {
  standardGeneric("filter_meta")
})

#' @rdname filter_meta
#' @export
setMethod("filter_meta", "omic", function(object, ..., .by) {
  
  # 1) Check for samples
  if (miss_sample(object)) {
    stop("Error: No samples available in the 'omic' object.")
  }
  
  # 2) Capture filter expressions
  exprs <- rlang::enquos(...)
  check_reserved_keywords(exprs)
  check_forbidden_expressions(exprs)
  
  # 3) Create long abundance once, if necessary
  meta_df <- meta(object, .fmt = "tbl")
  long_abun <- long_abundance_join(object, get_abundance_keys(exprs))
  if(length(long_abun) && length(meta_df)){
    long_abun <- dplyr::left_join(long_abun, meta_df, by = "sample_id") 
  }
  
  # 4) Apply filter expressions
  filtered_lists <- vector("list", length(exprs))
  
  for (i in seq_along(exprs)) {
    # Identify variables in the current expression
    evars <- all.vars(exprs[[i]])
    is_abun_expr <- any(evars %in% c("abun", "rela", "norm"))
    taxa_evars <- evars[evars %in% taxa_vars(object)]
    
    if (length(taxa_evars) > 0L && !is_abun_expr) {
      cli::cli_abort(
        c(
          "x" = "You referenced {length(taxa_evars)} taxa variable{?s} ({.val {paste(taxa_evars, collapse = ', ')}}) without using any abundance.",
          "i" = "When taxa variables appear in an expression in mutate_meta, you must also use at least one of {.field abun}, {.field rela}, or {.field norm} to define the taxa-sample context.",
          "*" = "Expression: {.code {rlang::as_label(exprs[[i]])}}"
        ),
        class = "omic_mutate_meta_taxa_without_abundance"
      )
    } else if(length(taxa_evars) > 0L && is_abun_expr){
      long_abun <- taxa(object, .fmt = "tbl") %>%
        dplyr::select(tidyselect::all_of(c("taxa_id", taxa_evars))) %>%
        dplyr::right_join(long_abun, by = "taxa_id")
    }
    
    # Determine grouping columns
    meta_groups <- setdiff(get_group_omic(object), taxa_vars(object))
    if(length(meta_groups) == 0 && is_abun_expr) meta_groups <- "sample_id"
    meta_groups <- rlang::syms(meta_groups)
    
    # Decide if expression references abun/rela/norm
    if (is_abun_expr) {
      # Join abundance
      filtered_lists[[i]] <- long_abun %>%
        dplyr::group_by(!!!meta_groups) %>%
        dplyr::filter(!!exprs[[i]]) %>%
        dplyr::ungroup() %>%
        dplyr::select(-tidyselect::any_of(c("taxa_id","abun","rela","norm", taxa_evars))) %>%
        dplyr::distinct() %>%
        dplyr::pull("sample_id")
    } else {
      filtered_lists[[i]] <- meta_df %>%
        dplyr::group_by(!!!meta_groups) %>%
        dplyr::filter(!!exprs[[i]]) %>%
        dplyr::ungroup() %>%
        dplyr::pull("sample_id")
    }
  }
  
  # 6) Intersect sample sets from all expressions (AND logic)
  final_samples <- purrr::reduce(filtered_lists, intersect)
  
  # Reorder final_samples to match the original object’s sample order
  sample_ids_in_object <- sample_id(object)
  final_samples <- final_samples[order(match(final_samples, sample_ids_in_object))]
  
  # 7) Subset the object by the filtered samples and return
  return(object[final_samples, ])
})

#' @rdname filter_meta
#' @export
setMethod("filter_meta", "omics", function(object, ..., .by) {
  
  # 1) Check for samples
  if (miss_sample(object, "any")) {
    stop("Error: No samples available in at least one 'omic' object.")
  }
  
  # 2) Capture filter expressions
  exprs <- rlang::enquos(...)
  check_reserved_keywords(exprs)
  check_forbidden_expressions(exprs)
  
  # 3) Gather data frames
  meta_df   <- meta(object, .collapse = TRUE)
  long_abun <- long_abundance_join(object, get_abundance_keys(exprs))
  if(length(long_abun) && length(get_abundance_keys(exprs))){
    long_abun <- dplyr::left_join(long_abun, meta_df, by = c("omic", "sample_id")) 
  }
  
  # 4) Apply filter expressions
  filtered_lists <- vector("list", length(exprs))
  
  for (i in seq_along(exprs)) {
    
    # Identify variables in the current expression
    evars <- all.vars(exprs[[i]])
    is_abun_expr <- any(evars %in% c("abun", "rela", "norm"))
    taxa_evars <- evars[evars %in% taxa_vars(object, .fmt = "unique")]
    
    if (length(taxa_evars) > 0L && !is_abun_expr) {
      cli::cli_abort(
        c(
          "x" = "You referenced {length(taxa_evars)} taxa variable{?s} ({.val {paste(taxa_evars, collapse = ', ')}}) without using any abundance.",
          "i" = "When taxa variables appear in an expression in mutate_meta, you must also use at least one of {.field abun}, {.field rela}, or {.field norm} to define the taxa-sample context.",
          "*" = "Expression: {.code {rlang::as_label(exprs[[i]])}}"
        ),
        class = "omic_mutate_meta_taxa_without_abundance"
      )
    } else if(length(taxa_evars) > 0L && is_abun_expr){
      long_abun <- taxa(object, .collapse = TRUE) %>%
        dplyr::select(tidyselect::all_of(c("omic", "taxa_id", taxa_evars))) %>%
        dplyr::right_join(long_abun, by = dplyr::join_by(omic, taxa_id))
    }
    
    # Determine grouping columns
    meta_groups <- setdiff(get_group_omic(object), taxa_vars(object, .fmt = "unique"))
    if(length(meta_groups) == 0 && is_abun_expr) meta_groups <- c("omic", "sample_id")
    meta_groups <- rlang::syms(meta_groups)
    
    # Decide if expression references abun/rela/norm
    if (is_abun_expr) {
      # Join abundance data with metadata
      filtered_lists[[i]] <- long_abun %>%
        dplyr::left_join(meta_df, by = c("omic","sample_id")) %>%
        dplyr::group_by(!!!meta_groups) %>%
        dplyr::filter(!!exprs[[i]]) %>%
        dplyr::ungroup() %>%
        dplyr::select(-tidyselect::any_of(c("taxa_id","abun","rela","norm", taxa_evars))) %>%
        dplyr::distinct() %>%
        dplyr::select(tidyselect::all_of(c("omic", "sample_id")))
    } else {
      # Filter metadata alone
      filtered_lists[[i]] <- meta_df %>%
        dplyr::group_by(!!!meta_groups) %>%
        dplyr::filter(!!exprs[[i]]) %>%
        dplyr::ungroup() %>%
        dplyr::select(tidyselect::all_of(c("omic", "sample_id")))
    }
  }
  
  # 6) Intersect sample sets from all expressions (AND logic)
  final_samples <- purrr::reduce(filtered_lists, dplyr::semi_join, by = c("omic", "sample_id")) 
  
  # 7) Subset each omic in the omics by the final filtered samples
  for(i in names(object)){
    # Extract the sample_ids for this omic that passed all filters
    filtered_i <- dplyr::filter(final_samples, omic == i) %>% dplyr::pull("sample_id")
    
    # Reorder survived_ids to match the original sample order
    object_samples <- sample_id(object[[i]])
    matched_indices <- match(filtered_i, object_samples)
    matched_indices <- matched_indices[!is.na(matched_indices)]
    filtered_i <- object_samples[matched_indices]
    object[[i]] <- object[[i]][filtered_i, ]
  }
  
  # 8) Validate and return the updated omics
  validObject(object)
  return(object)
})






#' Filter `omic` or `omics` Objects by Taxa Metadata and Abundance Fields
#'
#' Filters the taxa in an `omic` or `omics` object based on user-specified 
#' conditions provided through `...`. Each condition can reference columns in the 
#' taxa metadata (`taxa`) or abundance-related fields (`abun`, `rela`, `norm`). 
#' When multiple conditions are provided, they are combined using a logical AND, 
#' meaning a taxon must satisfy **all** conditions to be retained.
#'
#' @param object An `omic` or `omics` object. For `omics`, each contained 
#'   `omic` can be filtered according to the combined conditions, with a column 
#'   `omic` distinguishing them in any joined data frames.
#' @param ... Unquoted expressions for filtering. These expressions can reference 
#'   columns in the taxa metadata (`meta`) or any of the abundance-related fields 
#'   (`abun`, `rela`, `norm`). Each expression is applied sequentially, and only 
#'   taxa that satisfy **all** expressions are retained.
#' @param .by A character vector of grouping columns. The default differs by object 
#'   type:
#'   \itemize{
#'     \item For `omic`, the default is `.by = "taxa_id"`.
#'     \item For `omics`, the default is `.by = c("omic", "taxa_id")`.
#'   }
#'   If `.by` is `NULL`, no grouping is applied, and conditions are evaluated 
#'   across all taxa globally. If `.by` is missing, any existing taxa grouping 
#'   (from `group_taxa()`) is used instead.
#'
#' @details
#' This function integrates closely with **tidyverse** conventions. Each expression 
#' in `...` is captured and evaluated using tidy evaluation, allowing you to write 
#' conditions similarly to `dplyr::filter()`. When an expression references 
#' abundance variables (`abun`, `rela`, or `norm`), the function automatically 
#' joins the taxa metadata with the relevant abundance data (`long_abun`) before 
#' filtering.
#'
#' **Grouping Logic**:
#' 1. **If you provide an explicit `.by`** (a non-empty character vector), those 
#'    columns define the grouping.
#' 2. **If `.by` is `NULL`**, no grouping is applied.
#' 3. **If `.by` is missing**:
#'    - **If the object has existing taxa grouping** (set via `group_taxa()`), that 
#'      grouping is used.
#'    - **Else**, if any filter expression references abundance variables, 
#'      grouping defaults to:
#'      \itemize{
#'        \item `c("omic", "taxa_id")` for `omics` objects.
#'        \item `"taxa_id"` for `omic` objects.
#'      }
#'    - **Otherwise**, no grouping is applied.
#'
#' **Multiple Conditions**:
#' Each condition in `...` is applied sequentially, narrowing down the set of 
#' remaining taxa. After processing all conditions, only those taxa that 
#' satisfied **all** conditions remain.
#'
#' **Behavior with `omics`**:
#' - The `omic` column in the joined data frames identifies which `omic` each taxon 
#'   belongs to.
#' - Filtering expressions can apply to metadata or abundance columns, referencing 
#'   each `omic`'s data.
#' - Taxa are filtered **independently** for each `omic` based on the conditions.
#'
#' **Return Value**:
#' - For an `omic` object, an updated `omic` containing only the taxa that 
#'   passed all filters.
#' - For a `omics`, an updated list of `omic` objects. Each `omic` is 
#'   subset to the remaining taxa.
#'
#' This approach provides a flexible, `dplyr`-like interface for complex filtering 
#' of both taxa metadata and abundance data, respecting the specified grouping context.
#'
#' @return An `omic` or `omics` object containing only the taxa (and their 
#'   corresponding abundance data, if any) that satisfy all provided filter conditions. 
#'   The rest are removed. If no taxa match, an object with zero taxa is returned 
#'   (but the same structure otherwise).
#'
#' @export
#' @aliases filter_taxa,omic-method filter_taxa,omics-method
setGeneric("filter_taxa", function(object, ..., .by) {
  standardGeneric("filter_taxa")
})

#' @rdname filter_taxa
#' @export
setMethod("filter_taxa", "omic", function(object, ..., .by) {
  
  # 1) Check for taxa
  if (miss_taxa(object)) {
    stop("Error: No taxa available in the 'omic' object.")
  }
  
  # 2) Capture filter expressions
  exprs <- rlang::enquos(...)
  check_reserved_keywords(exprs)
  check_forbidden_expressions(exprs)
  
  # 3) Gather data frames
  taxa_df   <- taxa(object, .fmt = "tbl")
  long_abun <- long_abundance_join(object, get_abundance_keys(exprs))
  if (length(long_abun) && length(get_abundance_keys(exprs))) {
    long_abun <- dplyr::left_join(long_abun, taxa_df, by = "taxa_id") 
  }
  
  # 4) Apply filter expressions
  filtered_lists <- vector("list", length(exprs))
  
  for (i in seq_along(exprs)) {
    # Identify variables in the current expression
    evars <- all.vars(exprs[[i]])
    is_abun_expr <- any(evars %in% c("abun", "rela", "norm"))
    meta_evars <- evars[evars %in% meta_vars(object)]
    
    if (length(meta_evars) > 0L && !is_abun_expr) {
      cli::cli_abort(
        c(
          "x" = "You referenced {length(meta_evars)} meta variable{?s} ({.val {paste(meta_evars, collapse = ', ')}}) without using any abundance.",
          "i" = "When meta variables appear in an expression in mutate_taxa, you must also use at least one of {.field abun}, {.field rela}, or {.field norm} to define the meta-sample context.",
          "*" = "Expression: {.code {rlang::as_label(exprs[[i]])}}"
        ),
        class = "omic_mutate_meta_taxa_without_abundance"
      )
    } else if(length(meta_evars) > 0L && is_abun_expr){
      long_abun <- meta(object, .fmt = "tbl") %>%
        dplyr::select(tidyselect::all_of(c("sample_id", meta_evars))) %>%
        dplyr::right_join(long_abun, by = "sample_id")
    }
    
    # Determine grouping columns
    taxa_groups <- setdiff(get_group_omic(object), meta_vars(object))
    if(length(taxa_groups) == 0 && is_abun_expr) taxa_groups <- "sample_id"
    taxa_groups <- rlang::syms(taxa_groups)
    
    # Decide if expression references abun/rela/norm
    if (any(evars %in% c("abun", "rela", "norm"))) {
      # Filter using abundance data
      filtered_lists[[i]] <- long_abun %>%
        dplyr::group_by(!!!taxa_groups) %>%
        dplyr::filter(!!exprs[[i]]) %>%
        dplyr::ungroup() %>%
        # Remove abundance columns to avoid duplication
        dplyr::select(-tidyselect::any_of(c("sample_id", "abun", "rela", "norm", meta_evars))) %>%
        dplyr::distinct() %>%
        dplyr::pull("taxa_id")
    } else {
      # Filter using taxa taxadata alone
      filtered_lists[[i]] <- taxa_df %>%
        dplyr::group_by(!!!taxa_groups) %>%
        dplyr::filter(!!exprs[[i]]) %>%
        dplyr::ungroup() %>%
        dplyr::pull("taxa_id")
    }
  }
  
  # 6) Intersect taxa sets from all expressions (AND logic)
  final_taxa <- purrr::reduce(filtered_lists, intersect)
  
  # Reorder final_taxa to match the original object’s taxa order
  taxa_ids_in_object <- taxa_id(object)
  final_taxa <- final_taxa[order(match(final_taxa, taxa_ids_in_object))]
  
  # 7) Subset the object by the filtered taxa and return
  return(object[, final_taxa])
})

#' @rdname filter_taxa
#' @export
setMethod("filter_taxa", "omics", function(object, ..., .by) {
  
  # 1) Check for taxa
  if (miss_taxa(object, "any")) {
    stop("Error: No taxa available in at least one 'omic' object.")
  }
  
  # 2) Capture filter expressions
  exprs <- rlang::enquos(...)
  check_reserved_keywords(exprs)
  check_forbidden_expressions(exprs)
  
  # 3) Gather data frames
  taxa_df   <- taxa(object, .collapse = TRUE)
  long_abun <- long_abundance_join(object, get_abundance_keys(exprs))
  if (length(long_abun) && length(get_abundance_keys(exprs))) {
    long_abun <- dplyr::left_join(long_abun, taxa_df, by = c("omic","taxa_id")) 
  }
  
  # 4) Apply filter expressions
  filtered_lists <- vector("list", length(exprs))
  
  for (i in seq_along(exprs)) {
    # Identify variables in the current expression
    evars <- all.vars(exprs[[i]])
    is_abun_expr <- any(evars %in% c("abun", "rela", "norm"))
    meta_evars <- evars[evars %in% meta_vars(object, .fmt = "unique")]
    
    if (length(meta_evars) > 0L && !is_abun_expr) {
      cli::cli_abort(
        c(
          "x" = "You referenced {length(meta_evars)} meta variable{?s} ({.val {paste(meta_evars, collapse = ', ')}}) without using any abundance.",
          "i" = "When meta variables appear in an expression in mutate_taxa, you must also use at least one of {.field abun}, {.field rela}, or {.field norm} to define the meta-sample context.",
          "*" = "Expression: {.code {rlang::as_label(exprs[[i]])}}"
        ),
        class = "omic_mutate_meta_taxa_without_abundance"
      )
    } else if(length(meta_evars) > 0L && is_abun_expr){
      long_abun <- meta(object, .collapse = TRUE) %>%
        dplyr::select(tidyselect::all_of(c("omic", "sample_id", meta_evars))) %>%
        dplyr::right_join(long_abun, by  = dplyr::join_by(omic, sample_id))
    }
    
    # Determine grouping columns
    taxa_groups <- setdiff(get_group_omic(object), meta_vars(object, .fmt = "unique"))
    if(length(taxa_groups) == 0 && is_abun_expr) taxa_groups <- c("omic", "sample_id")
    taxa_groups <- rlang::syms(taxa_groups)
    
    # Decide if expression references abun/rela/norm
    if(is_abun_expr) {
      # Join abundance data with taxadata
      filtered_lists[[i]] <- long_abun %>%
        dplyr::left_join(taxa_df, by = c("omic","taxa_id")) %>%
        dplyr::group_by(!!!taxa_groups) %>%
        dplyr::filter(!!exprs[[i]]) %>%
        dplyr::ungroup() %>%
        dplyr::select(-tidyselect::any_of(c("sample_id","abun","rela","norm",meta_evars))) %>%
        dplyr::distinct() %>%
        dplyr::select(tidyselect::all_of(c("omic", "taxa_id")))
    } else {
      # Filter taxadata alone
      filtered_lists[[i]] <- taxa_df %>%
        dplyr::group_by(!!!taxa_groups) %>%
        dplyr::filter(!!exprs[[i]]) %>%
        dplyr::ungroup() %>%
        dplyr::select(tidyselect::all_of(c("omic", "taxa_id")))
    }
  }
  
  # 6) Intersect sample sets from all expressions (AND logic)
  final_taxa <- purrr::reduce(filtered_lists, dplyr::semi_join, by = c("omic", "taxa_id")) 
  
  # 7) Subset each omic in the omics by the final filtered samples
  for(i in names(object)){
    # Extract the sample_ids for this omic that passed all filters
    filtered_i <- dplyr::filter(final_taxa, omic == i) %>% dplyr::pull("taxa_id")
    
    # Reorder survived_ids to match the original sample order
    object_taxa <- taxa_id(object[[i]])
    matched_indices <- match(filtered_i, object_taxa)
    matched_indices <- matched_indices[!is.na(matched_indices)]
    filtered_i <- object_taxa[matched_indices]
    object[[i]] <- object[[i]][, filtered_i]
  }
  
  # 8) Validate and return the updated omics
  validObject(object)
  return(object)
})