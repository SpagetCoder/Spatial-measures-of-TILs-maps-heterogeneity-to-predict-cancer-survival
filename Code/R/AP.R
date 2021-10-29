## clear memory and set up variables and packages ## 
rm(list = ls())
gc()
library(png)
library(apcluster)
library(imager)
library(clusterCrit)
library("writexl")
library("readxl")

## Set up variables

# path for TIL map folder
folder <- "" 
file_list <- list.files(path=folder, pattern="*.png",recursive=TRUE, full.names=FALSE)

# path to mmc2 - modified.xlsx file
my_data <- read_excel("")

wcd_sd <- c()
np_sd <- c()
np_mean <- c()
wcd_mean <- c()
tll_perc <- c()
clusters_num <- c()
patches_num <- c()
file_names <- c()
file_names_errors <- c()
error_cause <- c()
c_index <- c()
ball_hall <- c()
banfeld_raftery <- c()
counter = 1


for (i in file_list)
{
  
  tryCatch(
    {
      ## Read the data
      print(paste("Working on picture number:", counter, "out of", length(file_list), sep=" "))
      img <- readPNG(paste(folder,i,sep=""))
      img2 = img[,,1]
      img_coordinates <- which(img2==1, arr.ind = TRUE)
      
      # on the used machine 40000 was the limit for which pc could handle the calculations
      if(nrow(img_coordinates) >= 40000)
      {
        print(paste("RAM ammount is to small to conduct the calculations for image: ", i, sep = " "))
        file_names_errors <<- c(file_names_errors, i)
        error_cause <<- c(error_cause, "insufficient RAM")
        counter = counter + 1
        next
      }
      
      else
      {
        ## Save the pixels coordinates plot
        pdf(paste("", substr(i,1,23), ".pdf" , sep=""))
        plot(img_coordinates)
        dev.off()
        
        ## Calcualte the best input values based on the other excel file
        name = substr(i,1,12)
        num = as.double(my_data[[which(my_data[[1]]==name, arr.ind=TRUE),14]])/nrow(img_coordinates)
        
        if (is.null(num))
        {
          file_names_errors <- c(file_names_errors,i)
          error_cause <- c(error_cause, "not found in the excel")
          print("Patch is not in the excel")
          counter = counter + 1
          next
        }
        
        ## Calculate affinity propagation
        apres1a <- apclusterL(negDistMat(r=2), img_coordinates, frac = num, sweeps = 100) 
        #apres1a <- apcluster(negDistMat(r=2), img_coordinates, q = 0) 
        print(paste("Frac =", num, "cluster number:", length(apres1a), sep = " "))
        
        ## Save the apcluster plot
        pdf(paste("", substr(i,1,23), ".pdf" , sep=""))
        plot(apres1a,img_coordinates)
        dev.off()
        
        ### Calculate indexes
        x <- as.matrix(as.numeric(img_coordinates))
        xx <- matrix(x,nrow = nrow(img_coordinates), ncol = 2)
        cl <- kmeans(xx, length(apres1a))
        intIdx <- intCriteria(xx,cl$cluster,c("C_index","Ball_Hall","Banfeld_Raftery"))
        
        #### Save values ####
        file_names <- c(file_names, substr(i,1,23))
        patches_num <- c(patches_num, nrow(img_coordinates))
        clusters_num <- c(clusters_num, length(apres1a))
        tll_perc <- c(tll_perc, nrow(img_coordinates)/as.double(tabulate(img)) * 100)
        c_index <- c(c_index, as.double(intIdx[1]))
        ball_hall <- c(ball_hall, as.double(intIdx[2]))
        banfeld_raftery <- c(banfeld_raftery, as.double(intIdx[3]))
        wcd_mean <- c(wcd_mean, as.double(cl[5])/length(apres1a))
        wcd_sd <- c(wcd_sd, sd(cl[[4]]))
        np_mean <- c(np_mean, nrow(img_coordinates)/length(apres1a))
        np_sd <- c(np_sd, sd(cl[[7]]))
      }
    },
    
    error = function(e)
    {
      print(paste("Failed to calculate affinity propagation for:", i, sep = " "))
      file_names_errors <<- c(file_names_errors, i)
      error_cause <<- c(error_cause, "algorithm failed")
    })
  
  ## Save the data after each iteration in case of power shortage etc
  
  final_data_frame <- data.frame(file_names,patches_num,clusters_num,tll_perc,
                                 np_mean, np_sd,wcd_mean, wcd_sd,
                                 c_index,ball_hall,banfeld_raftery)
  final_errors <- data.frame(file_names_errors, error_cause)
  
  # path to where to store the results both for data and errors that occurred
  write_xlsx(final_data_frame,"")
  write_xlsx(final_errors,"")
  
  
  counter = counter + 1
  gc()
}

print("DONE")
