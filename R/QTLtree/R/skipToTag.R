
## read a tag and strip the first and last characters
.readnexttag <-
  function (conn) 
{
  a <- readChar(conn, n = 1)
  while ( length(a) > 0 && a != "<") {
    a <- readChar(conn, n = 1)
    }
  if (length(a) == 0) 
    return(FALSE)
  a <- readChar(conn, n = 1)
  b <- NULL
  while ( length(a) > 0 && a != ">") {
    b <- c(b, a)
    a <- readChar(conn, n = 1)
  }
  if (length(a) == 0) 
    return(FALSE)
  paste(b, collapse = "")
}
## read a tag and strip the first and last characters
.readtextandskip <- function(conn)
  {
    b <- NULL
    a <- readChar(conn,n=1)
    while (a != "<") {
      b <- c(b,a)
      a <- readChar(conn,n=1)
    }
    ## and skip to the end
    while (a != ">") {
      a <- readChar(conn,n=1)
    }
    paste(b,collapse="")
  }


`skipToTag` <-
function(conn,findtag)
  {
    tag <- .readnexttag(conn)
    while (tag!=findtag&&tag!=FALSE) {
      tag <- .readnexttag(conn)
    }
    if (tag==FALSE) return(FALSE)
    return(TRUE) 
  }

