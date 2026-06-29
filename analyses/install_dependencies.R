args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_dir <- if (length(file_arg) > 0) {
  dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), winslash = "/", mustWork = TRUE))
} else {
  getwd()
}

project <- normalizePath(file.path(script_dir, ".."), winslash = "/", mustWork = TRUE)
old_wd <- setwd(project)
on.exit(setwd(old_wd), add = TRUE)

source(file.path(project, "setup", "install_dependencies.R"))
