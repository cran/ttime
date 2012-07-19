phylo <-
function (data){

   size <- dim(data);              # Size of 'data'
   ev <- size[1];                  # Number of events
   sp <- (size[2] - 4)/3;          # Number of species
   names_data <- names(data);      # Columns Names in 'data'

   #Formatting the event data for the function 'phylo'
   dat <- data.frame(matrix(0,ev,sp));   
   for (i in 1:sp) {
    dat[,i] <- data[,(3*(i-1)+2)]; 
    colnames(dat)[i] <- names_data[3*(i-1)+2];
   }

   #Plot the dendrogram representing the phylogenetic proximity
   dendro  <- pvclust(dat, method.dist="euclidean", method.hclust = "complete", nboot = 1000);
   dev.new(); plot(dendro);
}
