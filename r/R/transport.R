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
# Content-Length. Headers are small, so they are read a byte at a time until the blank line, which keeps
# this free of any parsing library.
#' @noRd
read_request <- function(con, limit = 64e6) {
    header <- raw(0)
    repeat {
        byte <- readBin(con, "raw", 1L)
        if (!length(byte)) break
        header <- c(header, byte)
        n <- length(header)
        if (n >= 4L && identical(header[(n - 3L):n], as.raw(c(13, 10, 13, 10)))) break
        if (n > 1e6) stop_foodnet("the request headers are too long to be foodnet's")
    }
    lines <- strsplit(rawToChar(header), "\r\n", fixed = TRUE)[[1]]
    first <- if (length(lines)) lines[1] else ""
    named <- grep("^content-length:", lines, ignore.case = TRUE, value = TRUE)
    size <- if (length(named)) as.integer(trimws(sub("^[^:]*:", "", named[1]))) else 0L
    if (is.na(size) || size < 0 || size > limit) stop_foodnet("the request body is not a size we accept")
    body <- ""
    if (size > 0) {
        body <- rawToChar(readBin(con, "raw", size))
        Encoding(body) <- "UTF-8"
    }
    list(request = first, body = body)
}

# One accepted connection: answer it, and return the parameters when it carried them.
#' @noRd
serve_one <- function(server) {
    con <- socketAccept(server, blocking = TRUE, open = "a+b", timeout = 10)
    on.exit(try(close(con), silent = TRUE), add = TRUE)
    request <- read_request(con)
    if (!grepl("^POST ", request$request)) {
        # a browser or a port scan: say what this port is, and go on waiting for the real thing
        writeBin(charToRaw(http_response(
            "this port belongs to the foodnet R package; foodnet posts CRM parameters to it",
            type = "text/plain")), con)
        return(NULL)
    }
    payload <- tryCatch(fromJSON(request$body, simplifyVector = FALSE), error = function(e) NULL)
    if (is.null(payload)) {
        writeBin(charToRaw(http_response("{\"received\": false, \"error\": \"not JSON\"}",
                                         status = "400 Bad Request")), con)
        return(NULL)
    }
    crm <- as_foodnet_crm(payload)
    answer <- sprintf(paste0("{\"received\": true, \"taxa\": %d, \"resources\": %d, ",
                             "\"growth_rates\": %d, \"without_a_rate\": %d}"),
                      length(crm$taxa), length(crm$resources), sum(!is.na(crm$growth_rates)),
                      length(crm$caveats$without_a_rate))
    writeBin(charToRaw(http_response(answer)), con)
    crm
}

#' Receive CRM parameters from the foodnet page
#'
#' Opens a port on this machine and waits for foodnet's "Send to R" to post the parameters to it. The
#' page sends them to 127.0.0.1 and nowhere else, and this function answers one request and closes the
#' port again.
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
        crm <- serve_one(server)
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
