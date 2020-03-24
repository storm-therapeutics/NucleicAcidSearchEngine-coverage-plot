## © 2020 STORM Therapeutics Limited - All rights reserved
## Author: Hendrik Weisser (hendrik.weisser@stormtherapeutics.com)

library(viridisLite) # color palettes

## Directory where modification data files can be found; change the value to
## ".../OpenMS/share/OpenMS/CHEMISTRY" of your OpenMS installation:
openms.chemistry.dir <- NULL

## Utility functions:

#' Read an mzTab file containing NucleicAcidSearchEngine results
#'
#' @param path Path to the mzTab file
#'
#' @return List of data frames for the different mzTab parts
read.mzTab <- function(path) {
  process.buffer <- function(buffer, col.names, result) {
    data <- read.table(text=buffer, sep="\t", col.names=col.names, header=FALSE,
                       row.names=NULL, na.strings="null")
    result[[as.character(data[1, 1])]] <- data[, -1]
    result
  }
  content <- readLines(path)
  prefixes <- substr(content, 1, 3)
  result <- list()
  buffer <- c()
  col.names <- c("MTH", "name", "value")
  for (i in 1:length(content)) {
    if (content[i] == "")
      next
    if ((substr(prefixes[i], 3, 3) == "H") &&
        ((i == 1) || (content[i - 1] == ""))) {
      if (length(buffer) > 0)
        result <- process.buffer(buffer, col.names, result)
      col.names <- unlist(strsplit(content[i], "\t", fixed=TRUE))
      buffer <- c()
    }
    else {
      buffer <- c(buffer, content[i])
    }
  }
  process.buffer(buffer, col.names, result)
}

#' Split a (modified) oligonucleotide sequence given as a string into a list of
#' (modified) nucleotides
split.sequence <- function(seq) {
  stopifnot(length(seq) == 1)
  split <- unlist(strsplit(seq, ""))
  bracket.open <- which(split == "[")
  bracket.close <- which(split == "]")
  stopifnot(length(bracket.open) == length(bracket.close))
  if (length(bracket.open) > 0) {
    mods <- sapply(1:length(bracket.open), function(i) {
      paste(split[(bracket.open[i] + 1):(bracket.close[i] - 1)], collapse="")
    })
    split[bracket.open] <- mods
    remove <- unlist(lapply(1:length(bracket.open), function(i) {
      (bracket.open[i] + 1):(bracket.close[i])
    }))
  }
  else {
    remove <- c()
  }
  if (length(remove) > 0)
    split <- split[-remove]
  if (split[1] == "p") {
    split[2] <- paste0("p", split[2])
    split <- split[-1]
  }
  if (split[length(split)] == "p") {
    split[length(split) - 1] <- paste0(split[length(split) - 1], "p")
    split <- split[-length(split)]
  }
  split
}

#' Resolve ambiguities in sequence assignments due to different localisations
#' of modifications
resolve.ambiguities <- function(mztab, groups=list("m5C?"=c("C", "m5C"))) {
  mztab$OSM$sequence <- as.character(mztab$OSM$sequence)
  mztab$OLI$sequence <- as.character(mztab$OLI$sequence)
  parts <- split(mztab$OSM, mztab$OSM$spectra_ref, drop=TRUE)
  nrows <- sapply(parts, nrow)
  ambig <- parts[nrows > 1]
  unexpected <- c()
  resolved <- lapply(ambig, function(part) {
    split.seqs <- lapply(part$sequence, split.sequence)
    lengths <- sapply(split.seqs, length)
    stopifnot(all(lengths[1] == lengths[-1]))
    mat <- do.call(rbind, split.seqs)
    consensus <- apply(mat, 2, function(nucs) {
      unique.nucs <- sort(unique(nucs))
      ## ensure there are only expected differences:
      if (length(unique.nucs) > 1) {
        matches <- sapply(groups, function(group) all(unique.nucs %in% group))
        if (any(matches))
          return(names(matches)[which(matches)[1]])
        else {
          signature <- paste(unique.nucs, collapse=", ")
          ## only warn once about every unexpected group:
          if (!signature %in% unexpected) {
            warning("unexpected ambiguity: ", signature)
            unexpected <<- c(unexpected, signature)
          }
          return(paste(unique.nucs, collapse="|"))
        }
      }
      return(unique.nucs)
    })
    con.seq <- sapply(consensus, function(nuc)
      ifelse((nchar(nuc) > 2) | ((nchar(nuc) == 2) & !grepl("^p|p$", nuc)),
             paste0("[", nuc, "]"), nuc))
    con.seq <- paste(con.seq, collapse="")
    if (!(con.seq %in% mztab$OLI$sequence)) {
      rows <- mztab$OLI[mztab$OLI$sequence %in% part$sequence, ]
      rows$sequence <- con.seq
      mztab$OLI <<- rbind(mztab$OLI, unique(rows))
    }
    part[1, "sequence"] <- con.seq
    part[1, , drop=FALSE]
  })
  unambig <- parts[nrows == 1]
  mztab$OSM <- do.call(rbind, c(unambig, resolved))
  mztab$OLI <- mztab$OLI[mztab$OLI$sequence %in% mztab$OSM$sequence, ]
  mztab$OSM$sequence <- as.factor(mztab$OSM$sequence)
  mztab$OLI$sequence <- as.factor(mztab$OLI$sequence)
  mztab
}

#' Read in RNA modification definitions (MODOMICS format)
#'
#' @param dir Directory containing modification data file(s)
#' @param include.custom Include custom mod. definitions?
#'
#' @return Data frame containing modification definitions
read.RNA.modifications <- function(dir=openms.chemistry.dir,
                                   include.custom=TRUE) {
  if (is.null(dir))
    stop("directory name required")
  col.classes <- rep("character", 9)
  col.classes[c(4, 7)] <- "factor"
  col.classes[8:9] <- "numeric"
  mod.info <- read.delim(file.path(dir, "Modomics.tsv"), quote="",
                         encoding="UTF-8", skip=4, colClasses=col.classes,
                         na.strings="None")
  if (include.custom) {
    col.classes <- c(col.classes, "character")
    custom.path <- file.path(dir, "Custom_RNA_modifications.tsv")
    custom.mods <- read.delim(custom.path, quote="", encoding="UTF-8", skip=4,
                              colClasses=col.classes)
    mod.info$alternatives <- ""
    mod.info <- rbind(mod.info, custom.mods)
  }
  ## use just "Q" (instead of "QtRNA") for queuosine:
  mod.info$short_name <- sub("QtRNA$", "Q", mod.info$short_name)
  mod.info
}


## Coverage plotting functions:

#' Generate a discrete color scale for a range of values
#'
#' @param nmax Maximum value
#' @param nmin Minimum value
#' @param n.colors Number of colors to use
#'
#' @return Vector of break points named with color codes
get.color.scale <- function(nmax, nmin=1, n.colors=8) {
  stopifnot(nmin <= nmax)
  nmax <- nmax - nmin + 1 # shift range to start at 1
  if (nmax <= n.colors) {
    scale <- 1:nmax
  }
  else {
    exp <- n.colors - 1
    if (nmax < ceiling(1.5^exp)) { # use a linear scale
      scale <- 1
      parts <- exp
      nums <- 2:nmax
      while (parts >= 2) {
        size <- floor(length(nums) / parts)
        scale <- c(scale, nums[size])
        nums <- nums[-(1:size)]
        parts <- parts - 1
      }
      scale <- c(scale, nmax)
    }
    else { # use an exponential scale
      base <- nmax^(1/exp)
      scale <- ceiling(round(base^(0:exp), 1))
    }
  }
  scale <- scale + nmin - 1 # shift range back if necessary
  names(scale) <- rev(magma(length(scale)))
  scale
}

#' Generate a table of "stacked bars" representing oligonucleotides
#'
#' @param oligo.data Oligonucleotide data (from mzTab)
#' @param osm.data Spectrum-match data (from mzTab)
#' @param nuc.length Length of the full RNA sequence
#' @param color.scale Color scale (from [get.color.scale()])
#' @param mod.info Table with RNA modification data (from MODOMICS)
#' @param mods User-defined mapping of modifications to symbols
#'
#' @return Matrix of HTML table cells
make.coverage.table <- function(oligo.data, osm.data, nuc.length, color.scale,
                                mod.info, mods=c()) {
  html.esc <- matrix(c("&", "&amp;", "<", "&lt;", ">", "&gt;", "\"", "&quot;"),
                     ncol=2, byrow=TRUE)
  oligo.data <- oligo.data[order(oligo.data$start, -oligo.data$end,
                                 oligo.data$sequence), ]
  oligo.seqs <- unique(as.character(oligo.data$sequence))
  split.seqs <- lapply(oligo.seqs, split.sequence)
  names(split.seqs) <- oligo.seqs
  table <- matrix("", nrow=1, ncol=nuc.length)
  for (i in seq_len(nrow(oligo.data))) {
    seq <- as.character(oligo.data[i, "sequence"])
    count <- sum(osm.data$sequence == seq)
    bucket <- which(color.scale >= count)[1]
    color <- names(color.scale)[bucket]
    split.seq <- split.seqs[[seq]]
    start <- oligo.data[i, "start"]
    end <- oligo.data[i, "end"]
    oligo.length <- end - start + 1
    if (length(split.seq) != oligo.length) {
      print(oligo.data[i, ])
      stop("length mismatch")
    }
    row.ind <- 0
    for (j in 1:nrow(table)) { # find a "free" row
      if (all(table[j, start:end] == "")) {
        row.ind <- j
        break
      }
    }
    if (row.ind == 0) { # no "free" row found -> add a new row
      row <- rep("", nuc.length)
      table <- rbind(table, row)
      row.ind <- nrow(table)
    }
    if (start == end) {
      table[row.ind, start] <- "<td class=\"left right\""
    }
    else {
      table[row.ind, start] <- "<td class=\"left\""
      table[row.ind, end] <- "<td class=\"right\""
    }
    if (oligo.length > 2) {
      table[row.ind, (start + 1):(end - 1)] <- "<td class=\"inner\""
    }
    for (j in 1:length(split.seq)) {
      if (grepl("^p?[ACGU]p?$", split.seq[j])) # no mod
        mod <- ""
      else {
        if (split.seq[j] %in% names(mods)) # special mod from user list
          char <- mods[split.seq[j]]
        else if (grepl("\\|", split.seq[j])) # ambiguous mod
          char = "&nbsp;"
        else {
          split.seq[j] <- sub("p$", "", split.seq[j]) # remove 3'-p, if present
          pos <- match(split.seq[j], mod.info$short_name)
          if (is.na(pos)) # unrecognized mod
            char <- "@"
          else {
            char <- mod.info[pos, "html_abbrev"]
          }
        }
        mod <- paste0("<div class=\"mod\" title=\"", split.seq[j], "\">", char,
                      "</div>")
      }
      table[row.ind, start + j - 1] <-
        paste0(table[row.ind, start + j - 1], " style=\"background-color:",
               color, "\" title=\"spectral count: ", count, "\">", mod,
               "</td>\n")
    }
  }
  table
}

#' Return the HTML page header for the coverage plot
get.html.header <- function() {
"<!doctype html>
<html>
<head>
<style>
body {
  font-family: sans-serif;
}
table, th, td {
  border: none;
  border-spacing: 0px 4px;
}
.scale {
  height: 10px;
  width: 40px;
}
.nums {
  writing-mode: sideways-lr;
  font-size: 50%;
}
.left {
  border-left: 2px solid white;
}
.right {
  border-right: 2px solid white;
}
.left, .inner, .right {
  height: 10px;
}
.mod {
  max-height: 10px;
  text-align: center;
  background-color: green;
  font-weight: bold;
  color: white;
  font-size: 10px;
  line-height: 10px;
}
</style>
</head>
<body>
"
}

#' Return HTML code for the color scale (from [get.color.scale()])
get.scale.html <- function(color.scale) {
  html <- "<table style=\"text-align:center\">\n<tr>\n<td>Spectral counts:</td>\n"
  lower <- c(color.scale[1], color.scale[-length(color.scale)] + 1)
  for (i in 1:length(color.scale)) {
    lo <- lower[i]
    hi <- color.scale[i]
    text <- ifelse(lo == hi, as.character(hi), paste0(lo, "-", hi))
    html <- paste0(html, "<td>", text, "</td>\n")
  }
  html <- paste0(html, "</tr>\n<tr>\n<td/>\n")
  for (color in names(color.scale)) {
    html <- paste0(html, "<td class=\"scale\" style=\"background-color:", color,
                   "\"></td>\n")
  }
  html <- paste0(html, "</tr>\n</table>\n")
  html
}

#' Helper function to generate HTML code for the coverage plot of a single RNA
#'
#' @param accession Database accession of the RNA
#' @param mztab1 First mzTab file containing RNA data
#' @param mztab2 Optional second mzTab file containing RNA data
#' @param labels Labels to show for two mzTab files
#' @param break.at Add line breaks after this many bases in the sequence
#' @param color.scale Color scale (from [get.color.scale()])
#' @param mod.info Table with RNA modification data (from MODOMICS)
#' @param description Description of the RNA
#'
#' @return HTML code
make.coverage.html.single <- function(accession, mztab1, mztab2=NULL,
                                      labels=NULL, break.at=Inf, color.scale,
                                      mod.info, description="") {
  index <- match(accession, mztab1$NUC$accession)
  if (!is.na(index))
    split.seq <- split.sequence(as.character(mztab1$NUC[index, "opt_sequence"]))
  else {
    index <- match(accession, mztab2$NUC$accession)
    split.seq <- split.sequence(as.character(mztab2$NUC[index, "opt_sequence"]))
  }

  nums <- rep("", length(split.seq))
  fives <- seq(5, length(split.seq), by=5)
  nums[fives] <- fives

  oligo.data1 <- mztab1$OLI[mztab1$OLI$accession == accession, ]
  coverage.table1 <- make.coverage.table(oligo.data1, mztab1$OSM,
                                         length(split.seq), color.scale,
                                         mod.info)
  covered <- sum(apply(coverage.table1, 2, function(col) any(col != "")))
  coverage1 <- covered / length(split.seq)

  if (!is.null(mztab2)) {
    oligo.data2 <- mztab2$OLI[mztab2$OLI$accession == accession, ]
    coverage.table2 <- make.coverage.table(oligo.data2, mztab2$OSM,
                                           length(split.seq), color.scale,
                                           mod.info)
    covered <- sum(apply(coverage.table2, 2, function(col) any(col != "")))
    coverage2 <- covered / length(split.seq)
    row.labels <- paste0(labels, ":")
  }
  else {
    coverage.table2 <- coverage.table1
    row.labels <- c("", "")
  }

  ## write header:
  html <- paste0("<h1>", accession, "</h1>\n")
  if (!is.na(description) && (description != ""))
    html <- paste0(html, "<p>", description, "</p>\n")
  html <- paste0(html, "<p>Sequence coverage: ", round(coverage1 * 100, 2),
                 "%")
  if (!is.null(mztab2))
    html <- paste0(html, " (", labels[1], ") / ", round(coverage2 * 100, 2),
                   "% (", labels[2], ")")
  html <- paste0(html, "</p>\n")
  html <- paste0(html, "<table>\n")
  ## line breaks:
  break.at <- min(break.at, length(split.seq))
  parts <- ceiling(length(split.seq) / break.at)
  row.label.style <- "max-height:10px; line-height:10px; white-space:nowrap"
  for (part in 0:(parts - 1)) {
    current.cols <- (1:break.at) + (part * break.at)
    current.cols <- current.cols[current.cols <= length(split.seq)]
    if (!is.null(mztab2)) {
      ## oligonucleotides (1):
      html <- paste0(html, "<tr>\n<td style=\"", row.label.style,
                     "\" rowspan=\"", nrow(coverage.table1), "\"><b>",
                     row.labels[1], "</b></td>\n")
      for (i in nrow(coverage.table1):1) {
        for (j in current.cols) {
          s <- coverage.table1[i, j]
          html <- paste0(html, if (s == "") "<td/>\n" else s)
        }
        html <- paste0(html, "</tr>\n")
      }
      ## position numbers (1):
      html <- paste0(html, "<tr class=\"nums\">\n<td/>\n",
                     paste0("<td>", nums[current.cols], "</td>", collapse="\n"),
                     "\n</tr>\n")
      html <- paste0(html, "<tr>\n<td>Sequence:</td>\n")
    }
    else
      html <- paste0(html, "<tr>\n<td/>\n")
    ## RNA sequence:
    html <- paste0(html, paste0("<th>", split.seq[current.cols], "</th>",
                                collapse="\n"), "\n</tr>\n")
    ## position numbers (2):
    html <- paste0(html, "<tr class=\"nums\">\n<td/>\n",
                   paste0("<td>", nums[current.cols], "</td>", collapse="\n"),
                   "\n</tr>\n")
    ## oligonucleotides (2):
    html <- paste0(html, "<tr>\n<td style=\"", row.label.style, "\" rowspan=\"",
                   nrow(coverage.table2), "\"><b>", row.labels[2],
                   "</b></td>\n")
    for (i in 1:nrow(coverage.table2)) {
      for (j in current.cols) {
        s <- coverage.table2[i, j]
        html <- paste0(html, if (s == "") "<td/>\n" else s)
      }
      html <- paste0(html, "</tr>\n")
    }
    html <- paste0(html, "<tr>\n<td colspan=\"", break.at + 1,
                   "\"><hr/></td>\n</tr>\n")
  }
  html <- paste0(html, "</table>\n")
  html
}

#' Generate a full HTML coverage plot for RNA data in one or two mzTab files
#'
#' @param mztab1 First mzTab file containing RNA data
#' @param mztab2 Optional second mzTab file containing RNA data
#' @param labels Labels to show for two mzTab files
#' @param break.at Add line breaks after this many bases in the sequence
#' @param mod.info Table with RNA modification data (from MODOMICS)
#'
#' @return HTML code
make.coverage.html <- function(mztab1, mztab2=NULL, labels=NULL, break.at=Inf,
                               mod.info=read.RNA.modifications()) {
  html <- get.html.header()
  count.range <- range(table(mztab1$OSM$sequence))
  accessions <- sort(levels(mztab1$NUC$accession))
  if (!is.null(mztab2)) {
    stopifnot(length(labels) == 2)
    count.range <- range(c(count.range, range(table(mztab2$OSM$sequence))))
    accessions <- sort(union(accessions, levels(mztab2$NUC$accession)))
  }
  color.scale <- get.color.scale(count.range[2], count.range[1])
  html <- paste0(html, get.scale.html(color.scale))
  for (accession in accessions) {
    ## cat(accession, "\n")
    desc <- mztab1$NUC[match(accession, mztab1$NUC$accession), "description"]
    acc.html <- make.coverage.html.single(accession, mztab1, mztab2, labels,
                                          break.at, color.scale, mod.info, desc)
    html <- paste0(html, acc.html)
  }
  html <- paste0(html, "</body>\n</html>\n")
  html
}
