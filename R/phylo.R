phylo <-
function (data){
  tmp <- names(data); d <- dim(data); ev <- d[1]; sp <- (d[2] - 4)/3;
  dat <- matrix(0,ev,sp); tmp2 <- array(0,sp);
    for (i in 1:sp) {
      tmp1 <- 3*(i-1)+2; dat[,i] <- data[,tmp1]; 
      tmp2[i] <- tmp[tmp1];
      rm(tmp1);
    }
   rm(tmp); dat <- data.frame(dat);
   for (i in 1:sp) colnames(dat)[i] <- tmp2[i];
   tmp  <- pvclust(dat, method.dist="euclidean", method.hclust = "complete", nboot = 1000);
   dev.new(); plot(tmp); rm(tmp); 
}

