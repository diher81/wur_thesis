###fast summary function

###input
#output of lm()

###output
#matrix with -log10 transformed p-values and effects

###Description
# Fast summary function
# This function cannot handle NA's!


    fast.summ <- function(input){
                          p <- input$rank
                          rdf <- input$df.residual
                          Qr <- input$qr
                          n <- nrow(Qr$qr)

                          p1 <- 1L:p
                          r <- input$residuals
                          rss <- colSums(r^2)

                          resvar <- rss/rdf
                          R <- chol2inv(Qr$qr[p1, p1, drop = FALSE])
                          se <- sqrt(rep(diag(R),each=length(resvar)) * rep(resvar,times=length(diag(R))))
                          est <- input$coefficients[Qr$pivot[p1],]
                          tval <- t(t(est)/se)

                          pvals <- round(-log10(2 * pt(abs(tval),rdf, lower.tail = FALSE)),digits=2)
                          output <- cbind(t(pvals)[,-1],t(round(est,digits=5))[,-1])
                          return(output)
                         }