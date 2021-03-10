oldw <- getOption("warn")
options(warn = -1)

shhh <- suppressPackageStartupMessages # stops annoying warnings when loading libraries
library(tidyr)
library(plyr)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(MASS)
library(Matrix)
library(reshape2)
library(ape) # stats
library(vegan) # stats
library(RColorBrewer)
library(cocor)
library(DescTools)
library(reshape2)
library(grid)
library(ggplotify)



# read the csv data files into a dataframe
files = list.files(pattern="*.csv")
data = sapply(files, read.csv, simplify=FALSE) %>% bind_rows(.id = "id")

colnames(data)

# Select variables we need for analysis 
trial_vars<- c( "participant",
                "Colour_1", "Colour_2", "Colour1", "Colour2", 
                "similarity", "response_time", "catchnumber", "Ecc", "catchnumberprac", "catchresponse", "catchtrialorder", "screen_size_x","screen_size_y","viewerdistancecm", 'viewer_distance',"trialnumber","Ecc")

# Create data frame for trials 
dftrials <- subset(data, !is.na(Colour1))

# Label participant number from 1 - 15 
dftrials$ID <- NA
subjectlist <- unique(dftrials$participant)
k= 0
for (participant in subjectlist){
  k = k + 1
  dftrials$ID[dftrials$participant == participant] <- k
}

# Remove participant 10 
dftrials <- dftrials[which(dftrials$ID != 10),]

# changing color values from RGB to hex for graphing purpose
dftrials$Colour1 <- as.character(dftrials$Colour1)
dftrials$Colour1 <- revalue(dftrials$Colour1, 
                            c(  "1" = '#FF0000',
                                "2" = '#FFAA00',
                                "3" = '#AAFF00',
                                "4" = '#00FF00',
                                "5" = '#00FFA9',
                                "6" = '#00A9FF',
                                "7" = '#0000FF',
                                "8" = '#AA00FF',
                                "9" = '#FF00AA'))
dftrials$Colour2 <- as.character(dftrials$Colour2)
dftrials$Colour2 <- revalue(dftrials$Colour2, 
                            c(  "1" = '#FF0000',
                                "2" = '#FFAA00',
                                "3" = '#AAFF00',
                                "4" = '#00FF00',
                                "5" = '#00FFA9',
                                "6" = '#00A9FF',
                                "7" = '#0000FF',
                                "8" = '#AA00FF',
                                "9" = '#FF00AA'))

# colors for the labels
# red, orange, yellow, green, cyan, cyan-blue, blue, purple, pink
colors <- c('#FF0000','#FFAA00','#AAFF00','#00FF00','#00FFA9','#00A9FF','#0000FF','#AA00FF','#FF00AA')
# can change the way the plot line up
# red, pink, orange, purple, yellow, blue, green, cyan-blue, cyan
#colors <- c('#FF0000','#FF00AA','#FFAA00','#AA00FF','#AAFF00','#0000FF','#00FF00','#00A9FF','#00FFA9')
abcolors <- sort(colors) # this was messing up the asymmetry plot, maybe useful for some other stuff


# factor the dataframes for the plot function
dissimdata2 <- function(dftrials, colors){
  
  # refactor the levels so they can be plotted properly later if need be
  dftrials$Colour1 <- with(dftrials, factor(Colour1, levels = colors))
  dftrials$Colour2 <- with(dftrials, factor(Colour2, levels = colors))
  
  return(dftrials)
}


df2mat_asymmetry_temporal <- function(dftrials){
  
  datatemp <- dissimdata2(dftrials, colors)
  
  # aggregate over the remaining columns of interest
  nmdsdata <- aggregate(datatemp, by = list(datatemp$Colour1, datatemp$Colour2),FUN=mean)
  nmdsdata$Colour1 <- nmdsdata$Group.1
  nmdsdata$Colour2 <- nmdsdata$Group.2
  
  nmdsdata = subset(nmdsdata, select = c("Colour1","Colour2","similarity"))  # get rid of unnecessary columns
  nmdsmatrix <- spread(nmdsdata, Colour1, similarity) # convert the dataframe to a matrix
  nmdsmatrix <- data.matrix(nmdsmatrix) # change first column from colour to number(just some label stuff) 
  nmdsmatrix <- nmdsmatrix[,-1] # get rid of the labels in the first column, it messes up the code
  
  matdf<-  as.data.frame(nmdsmatrix - t(nmdsmatrix)) # calculate the asymmetry
  matdf$colorset <- c(colors) # adding additional column "colorset"
  num_colors <- length(colors)
  matdf <- matdf %>% gather(othercolor,asymmetry ,1:num_colors) # convert the matrix back to the data frame which has the 
  # column "colortset", "othercolor", "asymmetry"
  return(matdf)
}


# Convert dataframe into a list of asymmtery matrix for each subject
dissimgroup2matlist <- function(dftrials){
  subjectlist <- sort(unique(dftrials$ID)) # obtain a list of all the subjects
  mat.list <- list()
  k = 1
  for(ID in subjectlist){
    subjectdf <- dftrials[which(dftrials$ID==ID),]
    subject.mat <- df2mat_asymmetry_temporal(subjectdf)
    mat.list[[k]] <- subject.mat
    k = k + 1
  }
  return(mat.list)
}

mat.list <- dissimgroup2matlist(dftrials)

# Return a list of asymmtery values for each cell 
vals.list.fun <- function(mat.list){
  n.cells <- length(mat.list[[1]])
  n.subjects <- length(mat.list)
  n.qualia <- nrow(mat.list[[1]])
  k=1
  vals.list <- list()
  
  for(cell in 1:n.cells){
  v <- vector()
  for(subject in 1:n.subjects){
      v <- c(v,mat.list[[subject]][cell])
    }
  vals.list[[cell]] <- v
  }
    return(vals.list)   
}



# T-test 
t.test.mu.zero <- function(cell, out='t'){
  test <- t.test(cell,
                 alternative = c("two.sided"),
                 mu = 0,
                 conf.level = 0.95)
  t <- test$statistic
  if(out=='t'){
    return(t)
  } 
}

# Create a matrix of t-scores
t.mat <- function(mat.list){
  mat.vals <- vals.list.fun(mat.list)
  n.qualia <- nrow(mat.list[[3]])
  out.mat <- matrix(nrow=n.qualia,ncol=n.qualia)
  
  for(cell in 1:length(mat.vals)){
    t <- t.test.mu.zero(mat.vals[[cell]],'t')
    out.mat[cell] <- t
  }
  return(out.mat)
}

t.mat(mat.list)

# Plot the t-scores
t.plot <- function(mat.list){
  dissimgroup2matlist(dftrials)
  
  t.test.mat <- t.mat(mat.list)
  
  plot <- mat.summary.plot(t.test.mat, low.lim=-4,up.lim=4,
                           xlabel = '', ylabel = '', col1 = 'green', col2='red',overall.title=title,na.col='purple')
  return(plot)
}