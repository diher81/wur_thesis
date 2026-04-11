###Generates a yyyymmdd date

###input
#

###output
#date

###Description



date.now <- function(){
              return(paste(unlist(strsplit(unlist(strsplit(as.character(Sys.time()),split = " "))[1],
                                           split = "-")),
                           collapse = ""))
}