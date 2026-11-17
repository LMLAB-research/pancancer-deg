need <- function(packages) {
  absent <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(absent)) stop("Missing packages: ", paste(absent, collapse = ", "), ". See README setup.")
}
bind <- function(x) {
  x <- Filter(function(z) is.data.frame(z) && nrow(z), x)
  if (!length(x)) return(data.frame())
  cols <- unique(unlist(lapply(x, names)))
  x <- lapply(x, function(z) { for (nm in setdiff(cols, names(z))) z[[nm]] <- NA; z[cols] })
  rownames_out <- do.call(rbind, x); rownames(rownames_out) <- NULL; rownames_out
}
write_csv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write.csv(x, path, row.names = FALSE, na = "", fileEncoding = "UTF-8")
}
symbols <- function(path) {
  if (!file.exists(path)) stop("MISSING_INPUT: ", path)
  x <- read.csv(path, check.names = FALSE)
  if (!"gene_symbol" %in% names(x)) stop("Missing gene_symbol: ", path)
  z <- toupper(trimws(as.character(x$gene_symbol)))
  sort(unique(z[!is.na(z) & nzchar(z)]))
}
input_fingerprint <- function() {
  dirs <- c(cancer_input_dir, senescence_input_dir, overlap_root, legacy_ois_root)
  files <- sort(unique(unlist(lapply(dirs[nzchar(dirs)], function(d) list.files(d, "\\.csv$", recursive = TRUE, full.names = TRUE)))))
  list(files = setNames(as.character(tools::md5sum(files)), files),
       code = tools::md5sum(file.path(script_dir, c("config.R", "helpers.R", "prepare_inputs.R", "01_GO.R", "02_MSigDB.R", "03_Reactome.R", "04_Compile.R"))),
       config = list(signatures, padj_threshold, lfc_threshold, minimum_input_genes,
                     minimum_gene_set_size, maximum_gene_set_size, bh_cutoff,
                     go_ontology, legacy_ois_full_table, dirs))
}
prepare_all <- function() {
  files <- sort(list.files(cancer_input_dir, cancer_file_pattern, full.names = TRUE))
  cancer_names <- sub("_DGE_Raw_Results_Tumor_vs_Normal_Adjusted\\.csv$", "", basename(files))
  from_overlaps <- unlist(lapply(c(file.path(overlap_root, signatures), legacy_ois_root), function(d) {
    z <- list.dirs(d, recursive = FALSE, full.names = FALSE); z[grepl("^TCGA-", z)]
  }))
  cancers <- sort(unique(c(cancer_names, from_overlaps)))
  if (!length(cancers)) stop("No TCGA inputs or prepared overlap folders found. Check config.local.R. Searched: ", paste(c(cancer_input_dir, overlap_root, legacy_ois_root), collapse = "; "))
  load_one <- function(file, cancer = FALSE, background = NULL) {
    tryCatch({
      if (!file.exists(file)) stop("MISSING_INPUT: ", file)
      z <- if (cancer) prepare_cancer(file) else if (!is.null(background)) {
        if (!file.exists(background)) stop("MISSING_INPUT: ", background)
        prepare_signature(file, background)
      } else prepare_genes(file)
      z$source <- normalizePath(file)
      if (is.null(z$background)) z$background <- z$data$gene_symbol
      z$data$source_file <- z$source
      z$audit$source_file <- z$source
      z
    }, error = function(e) list(error = conditionMessage(e), source = file))
  }
  cancer <- setNames(lapply(cancers, function(n) {
    f <- files[match(n, cancer_names)]
    if (is.na(f)) return(list(error = paste("MISSING_INPUT: full cancer DEG table", n)))
    load_one(f, TRUE)
  }), cancers)
  sf <- list.files(senescence_input_dir, "^(senescence|quiescence)_[A-Za-z0-9_-]+_genes\\.csv$", full.names = TRUE, ignore.case = TRUE)
  sn <- toupper(sub("_genes\\.csv$", "", sub("^(senescence|quiescence)_", "", basename(sf), ignore.case = TRUE), ignore.case = TRUE))
  if (anyDuplicated(sn)) stop("Duplicate signature filenames")
  sen <- setNames(lapply(signatures, function(s) {
    i <- match(s, sn)
    if (is.na(i)) return(list(error = paste("MISSING_INPUT: signature", s)))
    bg <- sub("_genes\\.csv$", "_background.csv", sf[i], ignore.case = TRUE)
    if (s == "OIS" && legacy_ois_full_table && !file.exists(bg)) bg <- NULL
    load_one(sf[i], background = bg)
  }), signatures)
  list(cancers = cancers, cancer = cancer, senescence = sen)
}
metrics <- function(z, prefix) {
  if (is.null(z) || !is.null(z$error)) return(data.frame(gene_symbol = character()))
  x <- z$data
  cols <- intersect(c("gene_symbol", "gene_id", "ensembl_id", "gene", "gene_name", "log2FoldChange", "pvalue", "p_value", "padj", "duplicate_n", "duplicate_status", "source_file"), names(x))
  x <- x[cols]
  if (!"ensembl_id" %in% names(x)) x$ensembl_id <- rep(NA_character_, nrow(x))
  if (!"pvalue" %in% names(x)) x$pvalue <- if ("p_value" %in% names(x)) x$p_value else rep(NA_real_, nrow(x))
  names(x)[names(x) != "gene_symbol"] <- paste0(prefix, "_", names(x)[names(x) != "gene_symbol"])
  x
}
make_job <- function(kind, cancer, signature, direction, inputs) {
  cdat <- inputs$cancer[[cancer]]; sdat <- inputs$senescence[[signature]]
  tryCatch({
    if (kind == "overlap") {
      d <- strsplit(direction, "_", fixed = TRUE)[[1]]
      folder <- file.path(overlap_root, signature, cancer)
      if (!dir.exists(folder) && signature == "OIS" && nzchar(legacy_ois_root)) folder <- file.path(legacy_ois_root, cancer)
      # Prepared overlap files use cancer first; workbook directions use senescence first.
      file <- file.path(folder, paste0("overlap_", d[2], "-", d[1], ".csv"))
      bgfile <- file.path(folder, "common_tested_genes.csv")
      genes <- symbols(file); universe <- symbols(bgfile)
      if (length(setdiff(genes, universe))) stop("Overlap contains genes outside prepared background: ", file)
      if (is.null(cdat$error) && is.null(sdat$error) && !is.null(cdat$data) && !is.null(sdat$data)) {
        expected_bg <- intersect(cdat$data$gene_symbol, sdat$background)
        if (!setequal(universe, expected_bg)) stop("Prepared background differs from supplied raw inputs: ", bgfile)
      }
      raw <- read.csv(file, stringsAsFactors = FALSE, check.names = FALSE)
      raw$gene_symbol <- toupper(trimws(raw$gene_symbol))
      if (anyDuplicated(raw$gene_symbol)) stop("Duplicate symbols in prepared overlap: ", file)
      detail <- data.frame(gene_symbol = genes)
      for (system in c("cancer", "senescence")) {
        suffix <- if (system == "cancer") "CANCER" else if ("log2FoldChange_SIGNATURE" %in% names(raw)) "SIGNATURE" else "OIS"
        for (m in c("log2FoldChange", "pvalue", "padj")) {
          col <- paste0(m, "_", suffix)
          detail[[paste0(system, "_", m)]] <- if (col %in% names(raw)) raw[[col]][match(genes, raw$gene_symbol)] else rep(NA_real_, length(genes))
        }
      }
      for (nm in setdiff(names(raw), names(detail))) detail[[nm]] <- raw[[nm]][match(genes, raw$gene_symbol)]
      detail$overlap_source <- rep(normalizePath(file), nrow(detail))
      # Raw tables supply p-values and source IDs that upstream CSVs did not retain.
      for (pair in list(list(cdat, "cancer"), list(sdat, "senescence"))) {
        tab <- metrics(pair[[1]], pair[[2]])
        if (nrow(tab)) {
          ix <- match(detail$gene_symbol, tab$gene_symbol)
          for (nm in setdiff(names(tab), "gene_symbol")) {
            v <- tab[[nm]][ix]
            if (nm %in% names(detail)) {
              ok <- !is.na(detail[[nm]]) & !is.na(v)
              if (is.numeric(v) && any(abs(detail[[nm]][ok] - v[ok]) > 1e-8)) stop("Stale overlap metrics differ from raw input: ", nm)
              fill <- is.na(detail[[nm]]); detail[[nm]][fill] <- v[fill]
            } else detail[[nm]] <- v
          }
          z <- pair[[1]]$data
          if (anyNA(ix)) stop("Prepared overlap gene absent from supplied raw table")
          expected <- if (pair[[2]] == "cancer") d[2] else d[1]
          if (any(z$regulation[match(genes, z$gene_symbol)] != expected)) stop("Prepared overlap direction differs from supplied raw table")
        }
      }
      return(list(kind = kind, genes = genes, universe = universe, details = detail, source = c(file, bgfile)))
    }
    z <- if (kind == "cancer") cdat else sdat
    if (!is.null(z$error)) stop(z$error)
    genes <- z$data$gene_symbol[z$data$regulation == direction]
    list(kind = kind, genes = sort(unique(genes)), universe = sort(unique(z$background)),
         details = metrics(z, if (kind == "cancer") "cancer" else "senescence"), source = z$source)
  }, error = function(e) list(error = conditionMessage(e)))
}
job_grid <- function(inputs) {
  bind(list(
    transform(expand.grid(cancer = inputs$cancers, signature = signatures,
                          direction = c("UP_UP", "UP_DOWN", "DOWN_UP", "DOWN_DOWN"), stringsAsFactors = FALSE), kind = "overlap"),
    transform(expand.grid(cancer = inputs$cancers, direction = c("UP", "DOWN"), stringsAsFactors = FALSE), signature = "", kind = "cancer"),
    transform(expand.grid(signature = signatures, direction = c("UP", "DOWN"), stringsAsFactors = FALSE), cancer = "", kind = "senescence")
  ))
}
annotation <- function(db) {
  need(c("AnnotationDbi", "org.Hs.eg.db"))
  dir.create(cache_root, recursive = TRUE, showWarnings = FALSE)
  cache <- file.path(cache_root, paste0(db, "_annotations.rds"))
  if (file.exists(cache)) return(readRDS(cache))
  org <- org.Hs.eg.db::org.Hs.eg.db
  if (db == "GO") {
    need("GO.db")
    a <- suppressMessages(AnnotationDbi::select(org, keys = AnnotationDbi::keys(org, "SYMBOL"),
                 columns = c("GOALL", "ONTOLOGYALL"), keytype = "SYMBOL"))
    a <- a[!is.na(a$GOALL) & a$ONTOLOGYALL == go_ontology, ]
    terms <- AnnotationDbi::select(GO.db::GO.db, unique(a$GOALL), "TERM", "GOID")
    a <- data.frame(pathway_id = a$GOALL, gene_symbol = a$SYMBOL, collection = "GO", subcollection = go_ontology)
    a$pathway <- terms$TERM[match(a$pathway_id, terms$GOID)]
    version <- c(org = as.character(packageVersion("org.Hs.eg.db")), GO = as.character(packageVersion("GO.db")))
  } else if (db == "MSIGDB") {
    need("msigdbr")
    old_cache <- Sys.getenv("R_USER_CACHE_DIR", unset = NA_character_)
    on.exit(if (is.na(old_cache)) Sys.unsetenv("R_USER_CACHE_DIR") else Sys.setenv(R_USER_CACHE_DIR = old_cache), add = TRUE)
    Sys.setenv(R_USER_CACHE_DIR = file.path(cache_root, "downloads"))
    m <- rbind(msigdbr::msigdbr(species = "Homo sapiens", collection = "H"),
               msigdbr::msigdbr(species = "Homo sapiens", collection = "C2"))
    keep <- m$gs_collection == "H" | (m$gs_collection == "C2" & (m$gs_subcollection == "CGP" | grepl("^CP", m$gs_subcollection)))
    m <- m[keep & !grepl("REACTOME", m$gs_subcollection, ignore.case = TRUE) & !grepl("^REACTOME_", m$gs_name), ]
    a <- data.frame(pathway_id = m$gs_id, pathway = m$gs_name, gene_symbol = m$gene_symbol,
                    collection = m$gs_collection, subcollection = m$gs_subcollection)
    version <- unique(m$db_version)
  } else {
    need("reactome.db")
    r <- as.list(reactome.db::reactomeEXTID2PATHID)
    a <- data.frame(entrez_id = rep(names(r), lengths(r)), pathway_id = unlist(r, use.names = FALSE))
    a <- a[grepl("^R-HSA-", a$pathway_id), ]
    ids <- intersect(unique(a$entrez_id), AnnotationDbi::keys(org, "ENTREZID"))
    m <- suppressMessages(AnnotationDbi::select(org, ids, "SYMBOL", "ENTREZID"))
    a <- merge(a, m, by.x = "entrez_id", by.y = "ENTREZID")
    a$gene_symbol <- a$SYMBOL; a$SYMBOL <- NULL
    nm <- as.list(reactome.db::reactomePATHID2NAME)
    a$pathway <- vapply(a$pathway_id, function(id) paste(nm[[id]], collapse = "; "), character(1))
    a$collection <- "REACTOME"; a$subcollection <- "all_levels"
    version <- as.character(packageVersion("reactome.db"))
  }
  a$annotation_symbol <- a$gene_symbol
  a$gene_symbol <- toupper(trimws(a$gene_symbol))
  a <- unique(a[!is.na(a$gene_symbol) & nzchar(a$gene_symbol), ])
  a$signature_direction <- ifelse(grepl("_UP$", a$pathway), "UP", ifelse(grepl("_(DN|DOWN)$", a$pathway), "DN", NA_character_))
  if (db != "MSIGDB") a$signature_direction <- NA_character_
  result <- list(associations = a, version = version, retrieved = as.character(Sys.time()),
                 ontology = if (db == "GO") go_ontology else NA_character_,
                 org_metadata = AnnotationDbi::metadata(org),
                 source_metadata = if (db == "REACTOME") reactome.db::reactome_dbInfo() else NULL)
  saveRDS(result, cache); result
}
ora <- function(job, associations) {
  bg <- intersect(job$universe, unique(associations$gene_symbol))
  genes <- intersect(job$genes, bg)
  audit <- data.frame(gene_symbol = job$universe,
                      in_input = job$universe %in% job$genes,
                      annotated = job$universe %in% bg)
  if (length(genes) < minimum_input_genes) return(list(table = data.frame(), details = data.frame(), audit = audit,
    status = "SKIPPED_TOO_FEW_ELIGIBLE_GENES", message = paste(length(genes), "eligible genes; minimum", minimum_input_genes)))
  a <- associations[associations$gene_symbol %in% bg, ]
  sets <- lapply(split(a$gene_symbol, a$pathway_id), unique)
  sets <- sets[lengths(sets) >= minimum_gene_set_size & lengths(sets) <= maximum_gene_set_size]
  if (!length(sets)) return(list(table = data.frame(), details = data.frame(), audit = audit, status = "NO_ELIGIBLE_PATHWAYS", message = "No pathways pass eligible size limits"))
  hits <- lapply(sets, intersect, y = genes)
  ids <- names(sets); k <- lengths(hits); K <- lengths(sets); n <- length(genes); N <- length(bg)
  meta <- associations[match(ids, associations$pathway_id), c("pathway_id", "pathway", "collection", "subcollection", "signature_direction")]
  tab <- cbind(meta, data.frame(Overlap_count = k, Pathway_size = K, Input_size = n, Universe_size = N,
    Input_size_raw = length(unique(job$genes)), Universe_size_raw = length(unique(job$universe)),
    GeneRatio = paste0(k, "/", n), BgRatio = paste0(K, "/", N), Pathway_coverage = k/K,
    FoldEnrichment = (k/n)/(K/N), p_value = phyper(k - 1, K, N-K, n, lower.tail = FALSE), row.names = NULL))
  # Include zero-hit pathways in the BH family; they have p = 1.
  tab$p_adjust <- p.adjust(tab$p_value, "BH")
  tab$significant <- tab$p_adjust < bh_cutoff
  tab$matched_genes <- vapply(hits, paste, collapse = "; ", character(1))
  detail <- bind(lapply(seq_along(ids), function(i) {
    if (!length(hits[[i]])) return(data.frame())
    z <- data.frame(pathway_id = ids[i], gene_symbol = hits[[i]], stringsAsFactors = FALSE)
    z <- merge(z, job$details, by = "gene_symbol", all.x = TRUE, sort = FALSE)
    z$collection <- meta$collection[i]; z$subcollection <- meta$subcollection[i]; z
  }))
  if (nrow(detail)) {
    ix <- match(detail$pathway_id, tab$pathway_id)
    for (nm in c("pathway", "signature_direction", "Overlap_count", "p_value", "p_adjust")) detail[[paste0("ORA_", nm)]] <- tab[[nm]][ix]
    systems <- if (identical(job$kind, "cancer")) "cancer" else if (identical(job$kind, "senescence")) "senescence" else c("senescence", "cancer")
    text <- detail$gene_symbol
    for (system in systems) {
      nm <- paste0(system, "_log2FoldChange")
      if (!nm %in% names(detail)) detail[[nm]] <- NA_real_
      text <- paste0(text, " [", system, " log2FC=", detail[[nm]], "]")
    }
    text <- split(text, detail$pathway_id)
    tab$gene_details <- vapply(ids, function(id) paste(text[[id]], collapse = "\n"), character(1))
  } else tab$gene_details <- ""
  tab <- tab[order(tab$p_adjust, tab$p_value, tab$pathway_id), ]
  list(table = tab, details = detail, audit = audit,
       status = if (any(tab$significant)) "COMPLETED" else "NO_SIGNIFICANT_HITS", message = "")
}
run_database <- function(db) {
  dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
  inputs <- prepare_all(); grid <- job_grid(inputs); fingerprint <- input_fingerprint()
  ann <- tryCatch(annotation(db), error = function(e) list(error = conditionMessage(e)))
  if (is.null(ann$error) && db == "GO" && !identical(ann$ontology, go_ontology)) stop("Cached GO ontology differs; remove GO annotation cache")
  out <- vector("list", nrow(grid))
  for (i in seq_len(nrow(grid))) {
    g <- grid[i, ]; job <- make_job(g$kind, g$cancer, g$signature, g$direction, inputs)
    r <- tryCatch({
      if (!is.null(job$error)) stop(job$error)
      if (!is.null(ann$error)) stop("ANNOTATION_FAILED: ", ann$error)
      ora(job, ann$associations)
    }, error = function(e) list(table = data.frame(), details = data.frame(), audit = data.frame(),
         status = if (grepl("MISSING_INPUT", conditionMessage(e))) "MISSING_INPUT" else "FAILED", message = conditionMessage(e)))
    r$job <- g; r$source <- job$source; out[[i]] <- r
    message(db, " ", i, "/", nrow(grid), " ", paste(g, collapse = " / "), ": ", r$status, " ", r$message)
  }
  audit_dir <- file.path(output_root, "audit")
  for (system in c("cancer", "senescence")) for (nm in names(inputs[[system]])) {
    z <- inputs[[system]][[nm]]
    if (!is.null(z$audit)) write_csv(z$audit, file.path(audit_dir, paste0(system, "_", nm, "_rows.csv")))
  }
  if (is.null(ann$error)) {
    need(c("AnnotationDbi", "org.Hs.eg.db"))
    genes <- unique(unlist(lapply(out, function(r) r$audit$gene_symbol)))
    org_keys <- AnnotationDbi::keys(org.Hs.eg.db::org.Hs.eg.db, "SYMBOL")
    valid <- org_keys[toupper(org_keys) %in% genes]
    if (length(valid)) {
      mapping <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db, valid, c("ENTREZID", "ENSEMBL"), "SYMBOL"))
      mapping$annotation_symbol <- mapping$SYMBOL
      mapping$SYMBOL <- toupper(mapping$SYMBOL)
      mapping <- merge(data.frame(SYMBOL = genes), mapping, all.x = TRUE)
      write_csv(mapping, file.path(audit_dir, "symbol_entrez_ensembl_mapping.csv"))
    }
  }
  status <- bind(lapply(out, function(r) cbind(r$job, database = db, status = r$status, message = r$message)))
  for (i in seq_along(out)) if (nrow(out[[i]]$details) && exists("mapping", inherits = FALSE)) {
    d <- out[[i]]$details
    for (id in c("ENTREZID", "ENSEMBL")) {
      groups <- split(mapping[[id]], mapping$SYMBOL)
      d[[paste0("mapped_", tolower(id))]] <- vapply(d$gene_symbol, function(s) {
        x <- groups[[s]]; paste(sort(unique(x[!is.na(x)])), collapse = "; ")
      }, character(1))
    }
    for (system in c("cancer", "senescence")) {
      col <- paste0(system, "_ensembl_id")
      if (col %in% names(d)) d[[paste0(system, "_ensembl_mapping_status")]] <- vapply(seq_len(nrow(d)), function(j) {
        source <- d[[col]][j]; mapped <- d$mapped_ensembl[j]
        if (is.na(source) || !nzchar(source)) return("NO_SOURCE_ID")
        if (is.na(mapped) || !nzchar(mapped)) return("UNMAPPED")
        if (sub("\\.[0-9]+$", "", source) %in% strsplit(mapped, "; ", fixed = TRUE)[[1]]) "MATCH" else "MISMATCH_REVIEW"
      }, character(1))
    }
    out[[i]]$details <- d
  }
  result <- list(database = db, jobs = out, grid = grid, fingerprint = fingerprint, annotation = ann[setdiff(names(ann), "associations")], status = status)
  saveRDS(result, file.path(output_root, paste0(db, ".rds")))
  write_csv(status, file.path(output_root, paste0(db, "_status.csv")))
  capture.output(sessionInfo(), file = file.path(output_root, paste0(db, "_sessionInfo.txt")))
  if (any(status$status %in% c("FAILED", "MISSING_INPUT"))) stop("Results and explicit statuses saved, but inputs/analyses are incomplete. See ", db, "_status.csv")
  invisible(result)
}
excel_safe <- function(tab, csv) {
  for (nm in names(tab)) if (is.character(tab[[nm]])) {
    # Count bytes conservatively to cover Excel's UTF-16 character limit too.
    long <- !is.na(tab[[nm]]) & nchar(tab[[nm]], type = "bytes") > 30000L
    tab[[nm]][long] <- paste0("Full value: ", basename(dirname(csv)), "/", basename(csv), "; see GENE_DETAILS")
  }
  if (nrow(tab) > 1048575L) return(data.frame(status = "EXCEL_ROW_LIMIT", message = paste("Full table:", csv)))
  tab
}
compile_results <- function() {
  need("openxlsx")
  dbs <- c("GO", "MSIGDB", "REACTOME")
  runs <- setNames(lapply(dbs, function(db) {
    f <- file.path(output_root, paste0(db, ".rds"))
    if (!file.exists(f)) stop("Run missing analysis first: ", f)
    readRDS(f)
  }), dbs)
  if (!all(vapply(runs, function(z) identical(z$fingerprint, runs[[1]]$fingerprint) && identical(z$grid, runs[[1]]$grid), logical(1)))) stop("Analysis inputs/config differ across databases; rerun 01–03")
  grid <- runs[[1]]$grid
  books <- unique(grid[c("kind", "cancer", "signature")])
  excel <- file.path(output_root, "excel")
  for (i in seq_len(nrow(books))) {
    b <- books[i, ]; wb <- openxlsx::createWorkbook()
    folder <- if (b$kind == "senescence") "Senescence_individual" else b$cancer
    filename <- if (b$kind == "overlap") b$signature else paste0(if (b$kind == "cancer") sub("^TCGA-", "", b$cancer) else b$signature, "_individual")
    dest <- file.path(excel, folder, paste0(filename, ".xlsx"))
    details <- list()
    for (db in dbs) {
      indices <- which(grid$kind == b$kind & grid$cancer == b$cancer & grid$signature == b$signature)
      for (ix in indices) {
        r <- runs[[db]]$jobs[[ix]]; sheet <- paste0(db, "_", grid$direction[ix])
        tab <- r$table
        if (!nrow(tab)) tab <- data.frame(status = r$status, message = r$message)
        else { tab$status <- r$status; tab$message <- r$message }
        if (nrow(r$details)) details[[length(details) + 1L]] <- cbind(database = db, direction = grid$direction[ix], r$details)
        # Full tables always survive in CSV/RDS, even when Excel cannot hold a cell.
        csv <- file.path(excel, folder, paste0(filename, "_full"), paste0(sheet, ".csv"))
        write_csv(tab, csv)
        tab <- excel_safe(tab, csv)
        openxlsx::addWorksheet(wb, sheet)
        openxlsx::writeData(wb, sheet, tab, withFilter = TRUE)
        openxlsx::freezePane(wb, sheet, firstRow = TRUE)
        openxlsx::setColWidths(wb, sheet, cols = seq_len(ncol(tab)), widths = 20)
        longcols <- which(names(tab) %in% c("gene_details", "matched_genes"))
        if (length(longcols)) {
          openxlsx::setColWidths(wb, sheet, cols = longcols, widths = 65)
          openxlsx::addStyle(wb, sheet, openxlsx::createStyle(wrapText = TRUE, valign = "top"), rows = seq_len(nrow(tab)) + 1L, cols = longcols, gridExpand = TRUE)
        }
      }
    }
    d <- bind(details)
    if (!nrow(d)) d <- data.frame(status = "NO_MATCHED_GENES", message = "See analysis sheets for missing/skipped/failed status")
    sidecar <- file.path(excel, folder, paste0(filename, "_full"), "GENE_DETAILS.csv")
    write_csv(d, sidecar)
    d <- excel_safe(d, sidecar)
    openxlsx::addWorksheet(wb, "GENE_DETAILS")
    openxlsx::writeData(wb, "GENE_DETAILS", d, withFilter = TRUE)
    openxlsx::freezePane(wb, "GENE_DETAILS", firstRow = TRUE)
    openxlsx::setColWidths(wb, "GENE_DETAILS", seq_len(ncol(d)), widths = 23)
    dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)
    openxlsx::saveWorkbook(wb, dest, overwrite = TRUE)
    expected <- if (b$kind == "overlap") 13L else 7L
    if (length(openxlsx::getSheetNames(dest)) != expected) stop("Workbook validation failed: ", dest)
    message("Saved ", dest)
  }
  invisible(TRUE)
}
