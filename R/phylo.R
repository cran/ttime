phylo <-
function (out){
  tmp <- names(out); d <- dim(out); ev <- d[1]; sp <- (d[2] - 4)/3;
  dat <- matrix(0,ev,sp); tmp2 <- array(0,sp);
    for (i in 1:sp) {
      tmp1 <- 3*(i-1)+2; dat[,i] <- out[,tmp1]; 
      tmp2[i] <- tmp[tmp1];
      rm(tmp1);
    }
   rm(tmp); dat <- data.frame(dat);
   for (i in 1:sp) colnames(dat)[i] <- tmp2[i];

   par(mfrow = c(2,1));
   tmp <- hclust(dist(t(dat)), method = "complete")
   plot(tmp, frame.plot=FALSE, lwd=2, xlab = "", ylab = "", axes = FALSE, ann = FALSE); rm(tmp);
   tmp  <- pvclust(dat, method.dist="euclidean", method.hclust = "complete", nboot = 1000);
   plot(tmp, frame.plot=FALSE, lwd=2, xlab = "", ylab = "", axes = FALSE, ann = FALSE); rm(tmp);
}
