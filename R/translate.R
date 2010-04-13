translate <-
function (yo, npsp) {
      tmp <- dim(yo); dimx1 <- tmp[2]-1; namesx1 <- names(yo);
      y = yo[,2:dimx1];
      tmp <- dim(y); r <- tmp[1]; c <- tmp[2]; sp   <- c - 2;  ev <- r; 
      
tmp <- 0; for (i in 1:sp)
      tmp <- tmp + length(which(y[,i] > 0));
mat0 <- matrix(0,tmp,5); off <- 0;
      options(warn = -1);

for (i in 1:sp){
          ind <- which(y[,i] > 0 ); 
          len <- length(ind); 
          tmp1 <- off+1; tmp2 <- off+len;
          mat0[tmp1:tmp2,1] <- ind; 
          mat0[tmp1:tmp2,2] <- i; 
          mat0[tmp1:tmp2,3] <- y[ind,i];  
          mat0[tmp1:tmp2,4] <- y[ind,sp+1]*floor(i/(npsp+1)); 
          mat0[tmp1:tmp2,5] <- y[ind,sp+2]*floor(i/(npsp+1));
          off <- off + len; 
          rm(ind); rm(tmp1); rm(tmp2);
}

       x <- mat0; len <- length(unique(x[,1])); 

      tmp1 <- len + sp + 2 + 1; tmp2 <- dim(x);  num <- len + sp + 2;
      mat <- matrix(0,tmp2[1],tmp1); off <- 0;
      rm(tmp1); rm(tmp2);

for (i in 1:sp){
  ind <- which(x[,2] == i); 
  tmp1 <- off+1; tmp2 <- off+length(ind);
  mat[tmp1:tmp2,i] <- 1; colpos <- x[ind,1];
        rm(tmp1); rm(tmp2);
  for (j in 1:length(ind)){
    tmp1 <- off+j; 
                tmp2 <- sp + colpos[j];
    mat[tmp1,tmp2] <- 1;
               rm(tmp1); rm(tmp2);
  }
  off <- off + length(ind);
}
    tmp1 <- len+sp+1; tmp2 <- len+sp+2; tmp3 <- len+sp+3;
    mat[,tmp1] = x[,4]; 
    mat[,tmp2] = x[,5]; 
    mat[,tmp3] = log(x[,3]);

    write.table(mat,'tmp.txt', row.names = FALSE, col.names = FALSE); 
    inp <- read.table('tmp.txt', header = FALSE); unlink('tmp.txt');
    dm <- dim(inp); nspecies_confi <- 3*sp
    tmp <- names(inp); p <- paste('inp$',tmp[dm[2]], sep =''); inpval <- eval(parse(text = p)); rm(p); rm(tmp);
    mn <- floor(exp(min(inpval))); tmp <- mn*sp-1;
    cval_arr <- array(0,tmp); 

    for (i in 1:tmp){
      tmp1 <- (i-1)/10;
      tmp2 <- lm(formula = log(exp(inpval)-tmp1)~., data = inp[,1:num]);  
      tmp3 <- summary(tmp2);
      tmp4 <- sqrt(tmp3$r.squared);
      cval_arr[i] <- tmp3$r.squared;
      rm(tmp1); rm(tmp2); rm(tmp3);
    }
  rm(tmp);
  tmp <- which(cval_arr == max(cval_arr)); cval <- (tmp-1)/10; rm(tmp);
  fit1 <- lm(formula = log(exp(inpval)-cval)~., data = inp[,1:num]); 

  pred_val <- matrix(0,nspecies_confi,ev); 
  tmp1<- sp+1; tmp2 <- sp+2; tmp <- y[,tmp1:tmp2]; rm(tmp1); rm(tmp2);

for (i in 1:sp){
 if (i > 1) cat("Percent Completed:",floor((i/sp)*100), "%\n");
  for (j in 1:ev){
       arr <- array(0,dm[2]-1);               
       cnt1 <- i; 
             cnt2 <- sp + j; 
             cnt3 <- sp+ev+1; 
             cnt4 <- sp+ev+2;
             arr[cnt1] <- 1; arr[sp+j] <- 1; 
             arr[cnt3] <- tmp[j,1]*floor(i/(npsp+1)); ; 
             arr[cnt4] <- tmp[j,2]*floor(i/(npsp+1)); ;
       write.table(t(arr),'tmp.txt'); rdarr <- read.table('tmp.txt'); 
       tmp22 <- predict.lm(fit1,data.frame(rdarr), interval = "confidence");
       lwr <- 3*(i-1)+1; upr <- 3*i;
       pred_val[lwr:upr,j] <- tmp22;
       rm(tmp22); rm(arr);
  }
}
cat("Done........\n");

  unlink('tmp.txt');
  tpred <- exp(t(pred_val))+cval; 

  tmp <- cbind(data.frame(yo[,1]),tpred,yo[,dimx1-1],yo[,dimx1], yo[,dimx1+1]); 
  namesx2 <- array('x',3*sp+2); namesx1[1] ='';
  for (i in 1:sp){
   tmp1 <- 3*(i-1)+2; tmp2 <- i+1;
   namesx2[tmp1] <- namesx1[tmp2]; namesx2[tmp1+1] = "Lwr"; namesx2[tmp1+2] = "Upr";
  }
  namesx2[1] <- "Event Name"; namesx2[3*sp+2] <- "Cortical"; namesx2[3*sp+3] <- "Limbic"; namesx2[3*sp+4] <- "Reference";
  names(tmp) <- namesx2;
  out <- tmp; 
}


