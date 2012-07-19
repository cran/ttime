translate <-
function (event_data, npsp) {

      dim_data <- dim(event_data);        #Dimension of the given matrix data
      col_names <- names(event_data);     #Names of the columns in data
      y = event_data[,2:(dim_data[2]-1)]; #Matrix after trimming the columns (event names and references)

      #Number of species in the given data 
      sp   <- dim_data[2] - 4;  
      #Number of events in the given data
      ev <- dim_data[1];  
      
     #Number of known event timing
     event_count <- 0; 
     for (i in 1:sp) event_count <- event_count + length(which(y[,i] > 0));

#################################################################################
#   Mapping the matrix 'y' into a matrix 'mat' suitable for the regression
#################################################################################
     mat0 <- matrix(0,event_count,5); off <- 0;
     options(warn = -1);

     for (i in 1:sp){
          ind <- which(y[,i] > 0 ); len <- length(ind); 
          tmp1 <- off+1; tmp2 <- off+len;               #temporary variables 
          mat0[tmp1:tmp2,1] <- ind; 
          mat0[tmp1:tmp2,2] <- i; 
          mat0[tmp1:tmp2,3] <- y[ind,i];  
          mat0[tmp1:tmp2,4] <- y[ind,sp+1]*floor(i/(npsp+1)); 
          mat0[tmp1:tmp2,5] <- y[ind,sp+2]*floor(i/(npsp+1));
          off <- off + len; 
          rm(ind); rm(tmp1); rm(tmp2);                 
     } # end of loop for i

     mat <- matrix(0,event_count,(ev + sp + 2 + 1));

     off <- 0;                                                #temporary variable
     for (i in 1:sp){
  ind <- which(mat0[,2] == i); 
  mat[(off+1):(off+length(ind)),i] <- 1; colpos <- mat0[ind,1];
  for (j in 1:length(ind))
        mat[(off+j),(sp+colpos[j])] <- 1;
  off <- off + length(ind);
    }

    mat[,(ev+sp+1)] = mat0[,4]; 
    mat[,(ev+sp+2)] = mat0[,5]; 
    mat[,(ev+sp+3)] = log(mat0[,3]);

# Converting 'mat' into a data frame 
    mat <- as.data.frame(mat);


#################################################################################
#   Determining the optimum value of 'k'
#################################################################################

    dm <- dim(mat); 
    p <- paste('mat$',names(mat)[dm[2]], sep =''); inpval <- eval(parse(text = p)); rm(p);
    mn <- floor(exp(min(inpval)));

    #Array containing correlation between predicted and empirical vals
    step_size <- 10; corr_arr <- array(0,(mn*step_size)); 

    for (i in 1:(mn*step_size))
    corr_arr[i] <- sqrt(summary(lm(formula = log(exp(inpval)-((i-1)/step_size))~., data = mat[,1:(ev+sp+2)]))$r.squared); 
    
    k <- (which(corr_arr == max(corr_arr))-1)/10;


#################################################################################
#   The model
#################################################################################

 fit <- lm(formula = log(exp(inpval) - k)~., data = mat[,1:(ev+sp+2)]); 

########################################################################################
# Plot of the given empirical values against its predicted counterpart in the log-scale
########################################################################################

 tmp2 <- inpval; tmp3 <- log(exp(predict(fit)) + k);  #temporary variables 
 mn <- min(min(tmp2),min(tmp3)); mx <- max(max(tmp2),max(tmp3));
 dev.new();
 plot(tmp2,tmp3, ylim=c(mn, mx), xlim=c(mn, mx), xlab = "log(Post-Conceptional Days), Empirically Derived", ylab = "log(Post-Conceptional Days), Predicted");
 rm(tmp2); rm(tmp3);
 txt1 <- expression(paste("Adjusted R"^{2}, "=                                 "));
 txt2 <- paste(round(summary(fit)$adj.r.squared,2));
 text(mx-2,mx-0.1,txt1); text(mx-2,mx-0.1,txt2); 
 title('Scatter Plot', font.main = 1);
 


#################################################################################
#Predicting the unknown event timings
#################################################################################

pred_val <- matrix(0,(3*sp),ev); 

for (i in 1:sp){
 if (i > 1) cat("Percent Completed:",floor((i/sp)*100), "%\n");
  for (j in 1:ev){
       arr <- array(0,dm[2]-1);               
       arr[i] <- 1; arr[sp+j] <- 1; 
       arr[sp+ev+1] <- y[j,(sp+1)]*floor(i/(npsp+1));
       arr[sp+ev+2] <- y[j,(sp+2)]*floor(i/(npsp+1)); 
       pred_val[(3*(i-1)+1):(3*i),j] <- predict.lm(fit,as.data.frame(t(arr)), interval = "confidence");;
       rm(arr);
  } # end of loop for j
} # end of loop for i
cat("Done........\n");

# Predicted event timing values retruned by the function translate through 'results'

tpred <- exp(t(pred_val)) + k; 

#Output of the predicted and empirical values along with their confidence intervals in "results"
results <- cbind(data.frame(event_data[,1]),tpred,event_data[,(sp+2)],event_data[,(sp+3)], event_data[,(sp+4)]); 
col_names_out <- array('x',3*sp+2); col_names[1] ='';

  for (i in 1:sp){
     col_names_out[3*(i-1)+2] <- col_names[i+1]; 
     col_names_out[3*(i-1)+3] = "Lwr"; 
     col_names_out[3*(i-1)+4] = "Upr";
  }
  col_names_out[1] <- "Event Name"; col_names_out[3*sp+2] <- "Cortical"; col_names_out[3*sp+3] <- "Limbic"; col_names_out[3*sp+4] <- "Reference";
  names(results) <- col_names_out;
  
  results; 

}
