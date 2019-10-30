`simplexml` <-
function(conn,texttags=NULL,numbertags=NULL,vectags=NULL,mattags=NULL,opening)
  {
    tag <- .readnexttag(conn)
    res <-  list()
    while (opening!=substring(tag,2)) {
 #     print(names(res))
  #    cat("tag = ",tag,"\n")
      if (tag %in% texttags) {
        res[[tag]] <- .readtextandskip(conn)
      } else if (tag %in% numbertags) {
        res[[tag]] <- as.numeric(.readtextandskip(conn))
      } else if (tag %in% vectags) {
        res[[tag]] <- readTNT(conn)
        endtag <- .readnexttag(conn)
        if (tag != substring(endtag,2))
          stop("expected a matching end tag to ",tag,", got ",endtag)
      } else if (tag %in% mattags) {
        res[[tag]] <- readTNT2D(conn)    
         endtag <- .readnexttag(conn)
        if (tag != substring(endtag,2))
          stop("expected a matching end tag to ",tag,", got ",endtag)
      } else {
        endtag <- .readnexttag(conn)
        if (tag != substring(endtag,2))
          stop("expected a matching end tag to ",tag,", got ",endtag)
      }
      tag <- .readnexttag(conn)
    }
    return(res)
  }

