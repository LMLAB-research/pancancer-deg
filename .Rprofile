source("renv/activate.R")

# Claude Code sandbox proxy fix. Inside the sandbox, https_proxy points at a local
# proxy with credentials (http://user:pass@127.0.0.1:<port>). That proxy aborts the
# CONNECT when libcurl uses its default proxy auth (ANY), which breaks every request
# made through the R curl package (curl, httr -> TCGAbiolinks, httr2, readr). Base R's
# own libcurl is unaffected. Forcing basic proxy auth fixes it. Outside the sandbox
# the condition is false and nothing is changed.
local({
  proxy <- Sys.getenv("https_proxy", Sys.getenv("HTTPS_PROXY"))
  if (!grepl("^https?://[^/@:]+:[^/@]+@127\\.0\\.0\\.1:[0-9]+/?$", proxy)) {
    return(invisible())
  }
  basic_auth <- 1L # CURLAUTH_BASIC

  # curl handles: new handles default to basic proxy auth; options passed by the
  # caller (including proxyauth) still win. httr2, readr and curl::curl() use this.
  patch_curl <- function(...) {
    ns <- asNamespace("curl")
    if (isTRUE(attr(ns$new_handle, "sandbox_proxy_auth"))) return(invisible())
    original <- ns$new_handle
    patched <- function(...) {
      h <- original(proxyauth = basic_auth)
      if (...length() > 0) curl::handle_setopt(h, ...)
      h
    }
    attr(patched, "sandbox_proxy_auth") <- TRUE
    unlockBinding("new_handle", ns)
    assign("new_handle", patched, envir = ns)
    lockBinding("new_handle", ns)
  }

  # httr resets handle options after every request, so it also needs a global config
  # (re-applied per request). TCGAbiolinks queries and downloads go through httr.
  patch_httr <- function(...) httr::set_config(httr::config(proxyauth = basic_auth))

  if ("curl" %in% loadedNamespaces()) patch_curl() else
    setHook(packageEvent("curl", "onLoad"), patch_curl)
  if ("httr" %in% loadedNamespaces()) patch_httr() else
    setHook(packageEvent("httr", "onLoad"), patch_httr)
})
