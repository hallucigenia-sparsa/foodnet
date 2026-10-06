# The two ways parameters reach R, as in grownet's package: the page pushes them to a port this package
# opens, or this package fetches them from the page. Base R sockets only, so installing this package pulls
# in nothing but jsonlite.

#' @noRd
http_response <- function(body, status = "200 OK", type = "application/json") {
    paste0("HTTP/1.1 ", status, "\r\n",
           "Content-Type: ", type, "\r\n",
           "Content-Length: ", nchar(body, type = "bytes"), "\r\n",
           "Connection: close\r\n\r\n", body)
}

# Read one HTTP request from a connection: the request line, the headers, and the body named by
# Content-Length. Headers are small, so they are read a byte at a time into a fixed buffer until the blank
# line, which keeps this free of any parsing library and its cost linear in what arrives.
#' @noRd
read_request <- function(con, limit = 64e6, header_limit = 16384L) {
    header <- raw(header_limit)
    n <- 0L
    repeat {
        byte <- readBin(con, "raw", 1L)
        if (!length(byte)) break
        if (n >= header_limit) stop_foodnet("the request headers are too long to be foodnet's")
        n <- n + 1L
        header[n] <- byte
        if (n >= 4L && identical(header[(n - 3L):n], as.raw(c(13, 10, 13, 10)))) break
    }
    lines <- strsplit(rawToChar(header[seq_len(n)]), "\r\n", fixed = TRUE)[[1]]
    first <- if (length(lines)) lines[1] else ""
    field <- function(name) {
        found <- grep(paste0("^", name, ":"), lines, ignore.case = TRUE, value = TRUE)
        if (length(found)) trimws(sub("^[^:]*:", "", found[1])) else NA_character_
    }
    size <- suppressWarnings(as.integer(field("content-length")))
    if (is.na(size)) size <- 0L
    if (size < 0 || size > limit) stop_foodnet("the request body is not a size we accept")
    body <- ""
    if (size > 0) {
        body <- rawToChar(readBin(con, "raw", size))
        Encoding(body) <- "UTF-8"
    }
    list(request = first, body = body, token = field("x-foodnet-token"), origin = field("origin"))
}

# Where foodnet_listen() leaves the secret the page must send, readable by this user only. The foodnet page
# reads the same file (foodnet.rbridge.token_path computes R's tools::R_user_dir the same way), so a request
# from another machine, or from a web page open in a browser, cannot carry it: base R cannot bind a server
# socket to 127.0.0.1 alone, so the port is reachable from the network and the secret is what keeps it ours.
#' @noRd
listen_token_path <- function() file.path(tools::R_user_dir("foodnet", "cache"), "listen-token")

#' @noRd
new_listen_token <- function() {
    bytes <- if (file.exists("/dev/urandom")) {
        source <- file("/dev/urandom", "rb", raw = TRUE)
        on.exit(close(source), add = TRUE)
        readBin(source, "raw", 24L)
    } else {
        # no system source of randomness (Windows): R's generator, leaving the session's own stream untouched
        seed <- if (exists(".Random.seed", envir = globalenv())) get(".Random.seed", envir = globalenv())
        on.exit(if (is.null(seed)) rm(".Random.seed", envir = globalenv()) else
            assign(".Random.seed", seed, envir = globalenv()), add = TRUE)
        set.seed(NULL)
        as.raw(sample.int(256L, 24L, replace = TRUE) - 1L)
    }
    paste(sprintf("%02x", as.integer(bytes)), collapse = "")
}

#' @noRd
write_listen_token <- function(token, path = listen_token_path()) {
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    if (file.exists(path)) file.remove(path)
    file.create(path)
    Sys.chmod(path, "0600")
    writeLines(token, path)
    path
}

#' @noRd
answer <- function(con, body, status = "200 OK", type = "application/json") {
    writeBin(charToRaw(http_response(body, status = status, type = type)), con)
    NULL
}

# One accepted connection: answer it, and return the parameters when it carried them. Anything else is
# answered and the listener goes on waiting, so a stray request can neither stop it nor plant parameters.
#' @noRd
serve_one <- function(server, token) {
    con <- socketAccept(server, blocking = TRUE, open = "a+b", timeout = 10)
    on.exit(try(close(con), silent = TRUE), add = TRUE)
    request <- tryCatch(read_request(con), error = function(e) NULL)
    if (is.null(request)) {
        return(answer(con, "{\"received\": false, \"error\": \"not a request foodnet sends\"}",
                      status = "400 Bad Request"))
    }
    if (!grepl("^POST ", request$request)) {
        # a browser or a port scan: say what this port is, and go on waiting for the real thing
        return(answer(con, "this port belongs to the foodnet R package; foodnet posts CRM parameters to it",
                      type = "text/plain"))
    }
    if (!is.na(request$origin) || is.na(request$token) || !identical(request$token, token)) {
        # a web page (browsers name their origin) or a request without this session's secret
        return(answer(con, paste0("{\"received\": false, \"error\": \"this listener takes parameters from the ",
                                  "foodnet page on this machine only; update foodnet if it is the page\"}"),
                      status = "403 Forbidden"))
    }
    payload <- tryCatch(fromJSON(request$body, simplifyVector = FALSE), error = function(e) NULL)
    if (is.null(payload)) {
        return(answer(con, "{\"received\": false, \"error\": \"not JSON\"}", status = "400 Bad Request"))
    }
    crm <- tryCatch(as_foodnet_crm(payload), error = function(e) NULL)
    if (is.null(crm)) {
        return(answer(con, "{\"received\": false, \"error\": \"not foodnet CRM parameters\"}",
                      status = "400 Bad Request"))
    }
    answer(con, sprintf(paste0("{\"received\": true, \"taxa\": %d, \"resources\": %d, ",
                               "\"growth_rates\": %d, \"without_a_rate\": %d}"),
                        length(crm$taxa), length(crm$resources), sum(!is.na(crm$growth_rates)),
                        length(crm$caveats$without_a_rate)))
    crm
}

#' Receive CRM parameters from the foodnet page
#'
#' Opens a port on this machine and waits for foodnet's "Send to R" to post the parameters to it. The
#' page sends them to 127.0.0.1 and nowhere else. Base R cannot limit a listening port to this machine, so
#' the listener also writes a one-time secret to a file only this user can read, and accepts parameters
#' only from a request that carries it (the foodnet page reads the same file); a request from elsewhere, or
#' from a web page, is refused and the listener goes on waiting. It closes the port, and removes the
#' secret, once the parameters have arrived.
#'
#' @param port Port to listen on. The foodnet page sends to 8794 unless told otherwise (grownet's
#'   package uses 8793, so both can listen at once).
#' @param timeout Seconds to wait for the parameters before giving up.
#' @param quiet Set to `TRUE` to keep the function from printing what it is waiting for and what arrived.
#' @return The CRM parameters, an object of class `foodnet_crm`. Printing it shows its caveats; see
#'   [crm_consumed()], [crm_efficiency()] and [as_miasim()].
#' @seealso [foodnet_crm()], which fetches the same parameters when opening a port is not possible.
#' @examples
#' \dontrun{
#' crm <- foodnet_listen()          # then press Send to R on the foodnet page
#' crm
#' E <- crm_efficiency(crm)
#' }
#' @export
foodnet_listen <- function(port = 8794, timeout = 300, quiet = FALSE) {
    token <- new_listen_token()
    path <- write_listen_token(token)
    on.exit(unlink(path), add = TRUE)
    server <- serverSocket(port)
    on.exit(close(server), add = TRUE)
    if (!quiet) {
        message("foodnet: listening on http://127.0.0.1:", port, " for up to ", timeout,
                " seconds. Press Send to R on the foodnet page.")
    }
    deadline <- Sys.time() + timeout
    repeat {
        left <- as.numeric(difftime(deadline, Sys.time(), units = "secs"))
        if (left <= 0) {
            stop_foodnet("no parameters arrived within ", timeout, " seconds. Start again with ",
                         "foodnet_listen(), then press Send to R on the foodnet page.")
        }
        if (!isTRUE(socketSelect(list(server), timeout = min(left, 1))[1])) next
        crm <- serve_one(server, token)
        if (is.null(crm)) next
        if (!quiet) print(crm)
        return(invisible(crm))
    }
}

#' Fetch CRM parameters from the foodnet page
#'
#' Reads the parameters from the page's own address, for when opening a port is not possible. The page
#' shows the address under its CRM control; it carries the session token, so it works only on the machine
#' the page runs on. The `crm.json` inside a downloaded parameter zip works too.
#'
#' @param from The CRM address of a running foodnet page, for example
#'   `"http://127.0.0.1:8791/crm.json?token=..."`, or the path of a saved `crm.json`.
#' @param timeout Seconds to wait for the page.
#' @param quiet Set to `TRUE` to keep the function from printing what arrived.
#' @return The CRM parameters, an object of class `foodnet_crm`.
#' @seealso [foodnet_listen()], which receives them from the page's Send to R button.
#' @examples
#' \dontrun{
#' crm <- foodnet_crm("http://127.0.0.1:8791/crm.json?token=PASTE_THE_TOKEN")
#' }
#' @export
foodnet_crm <- function(from, timeout = 60, quiet = FALSE) {
    if (!is.character(from) || length(from) != 1 || !nzchar(from)) {
        stop_foodnet("give the CRM address of a running foodnet page, or the path of a saved crm.json")
    }
    if (grepl("^https?://", from)) {
        old <- options(timeout = timeout)
        on.exit(options(old), add = TRUE)
        connection <- base::url(from, open = "rb")
        on.exit(try(close(connection), silent = TRUE), add = TRUE)
        text <- paste(readLines(connection, warn = FALSE), collapse = "\n")
    } else {
        text <- paste(readLines(from, warn = FALSE), collapse = "\n")
    }
    crm <- as_foodnet_crm(fromJSON(text, simplifyVector = FALSE))
    if (!quiet) print(crm)
    invisible(crm)
}
