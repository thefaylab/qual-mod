#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)

#library(data.table); library(Rpath)
library(tidyverse)
library(LoopAnalyst)
library(permute)
library(igraph)
#library(ggrepel)


#-------------------------------------------------------------------------------
#User created functions

#-------------------------------------------------------------------------------

# Using matrix inverse
adjoint3 <- function(A) det(A)*solve(A)

#library(googlesheets)  
url <- "https://docs.google.com/spreadsheets/d/1-ukFR9B6idnCR_o2G4BXd25AD1ucAl-OPkjBcuHW3UA/edit#gid=1853171097"

lookup <- data.frame(response = c("Positive","Negative","None"),
                     value = c(1,-1,0))

matmatch <- data.frame(qnum = 1:11,
                       i = c(3,4,6,2,1,1,2,6,3,4,5),
                       j = c(6,5,5,6,6,7,7,7,7,7,7))  

#gs_ls()
#be <- gs_title("GCRL_model")
#lab <- gs_ws_ls(be)
#mydat <- as_tibble(gs_read(ss=be, ws = lab[1], skip=0))
mydat <- googlesheets4::read_sheet(url)

nobs <- nrow(mydat)
nudat <- mydat %>% 
   mutate(id = 1:nrow(mydat)) %>% 
   select(id,everything()) %>% 
   gather(key = "question", value = "response", 3:(ncol(mydat)+1)) %>% 
   mutate(qnum = rep(1:11, each=nobs)) %>% 
   left_join(lookup) %>% 
   left_join(matmatch)

B <- matrix(c(-1,-1, 0, 0, 0, 0, 0,
              1,-1,-1, 0, 0, 0, 0,
              0, 1,-1,-1,-1, 0, 0,
              0, 0, 1,-1, 0, 0, 0,
              0, 0, 1, 0,-1, 0, 0,
              0, 0, 0, 0, 0,-1, 0,
              0, 0, 0, 0, 0, 0,-1), byrow=TRUE, nrow = 7)

Wstore <- matrix(NA,nrow=7,ncol=nobs)
Astore <- matrix(NA,nrow=7,ncol=nobs)
evaldat <- nudat
for (k in unique(nudat$id)) {
   A <- B
   temp <- filter(nudat,id==k)
   for (irow in 1:nrow(temp))
      A[temp$i[irow],temp$j[irow]] <- temp$value[irow]
   #}  
   
   adj_A <- adjoint3(A)
   adj_A
   Tmat <- LoopAnalyst::make.T(A,status=TRUE)
   #system.time(Tmat <- LoopAnalyst::make.T(A,status=TRUE))
   
   
   Wmat <- abs(adj_A)/Tmat
   Wmat
   
   adj_AA <- adj_A
   
   #image showing sign of adjoint and weights
   colfunc <- colorRampPalette(c("white", "steelblue"))
   image(1:7,1:7,t(Wmat[7:1,]),col = colfunc(7))
   text(1,0,"+",col="white",cex=2)
   pick <- which(adj_AA>0)
   x <- ceiling(pick/7)
   y <- 8-pick%%7
   y[y==8] <- 1
   text(x,y,"+",col="white",cex=1.5)
   pick <- which(adj_AA<0)
   x <- ceiling(pick/7)
   y <- 8-pick%%7
   y[y==8] <- 1
   text(x,y,"-",col="white",cex=1.5)
   
   Wstore[,k] <- rev(Wmat[,7])
   Astore[,k] <- rev(adj_AA[,7])
   
}


#group_pic
#image showing sign of adjoint and weights
colfunc <- colorRampPalette(c("white", "steelblue"))
image(1:nobs,1:7,t(Wstore),col = colfunc(5))
#text(1,0,"+",col="white",cex=2)
pick <- which(Astore>0)
x <- ceiling(pick/7)
y <- 8-pick%%7
y[y==8] <- 1
text(x,y,"+",col="white",cex=1.5)
pick <- which(Astore<0)
x <- ceiling(pick/7)
y <- 8-pick%%7
y[y==8] <- 1
text(x,y,"-",col="white",cex=1.5)



g1 <- graph_from_adjacency_matrix( abs(t(A)) , diag = FALSE)
V(g1)$size = 20
plot(g1,edge.arrow.size=.4)



# Define UI for application that draws a histogram
ui <- fluidPage(
   
   # Application title
   titlePanel("GCRL: Qualitative Network Modeling"),
   
      mainPanel(
        
        # Output: Tabset w/ plot, summary, and table ----
        tabsetPanel(type = "tabs",
                    tabPanel("Networks", plotOutput("networkplot")),
                    #tabPanel("Perturbations", plotOutput("perplot")),
                    tabPanel("Evaluation", plotOutput("evalplot"))#,
                    # tabPanel("Summary", verbatimTextOutput("summary")),
                    #tabPanel("Summary Table", tableOutput("table"))
        )
        #plotOutput("EwEplot")
      )
  # )
)

# Define server logic 
server <- function(input, output) {
   
   output$networkplot <- renderPlot({
     
      if (length(unique(nudat$id))>2) 
         nudat <- filter(nudat,id>2)
      pick <- sample(unique(nudat$id),size = min(4,length(unique(nudat$id))), replace=FALSE)
      nudat <- filter(nudat, id %in% pick)
      nobs <- length(unique(nudat$id))
      #Wstore <- matrix(NA,nrow=7,ncol=nobs)
      #Astore <- matrix(NA,nrow=7,ncol=nobs)
      k <- unique(nudat$id)[1]
         A <- B
         temp <- filter(nudat,id==k)
         for (irow in 1:nrow(temp))
            A[temp$i[irow],temp$j[irow]] <- temp$value[irow]
      g1 <- graph_from_adjacency_matrix( abs(t(A)) , diag = FALSE)
      V(g1)$size = 20
      #plot(g1,edge.arrow.size=.4)
      k <- unique(nudat$id)[2]
      A <- B
      temp <- filter(nudat,id==k)
      for (irow in 1:nrow(temp))
         A[temp$i[irow],temp$j[irow]] <- temp$value[irow]
      g2 <- graph_from_adjacency_matrix( abs(t(A)) , diag = FALSE)
      V(g2)$size = 20
      #plot(g1,edge.arrow.size=.4)
      if (nobs>2) {
      k <- unique(nudat$id)[3]
      A <- B
      temp <- filter(nudat,id==k)
      for (irow in 1:nrow(temp))
         A[temp$i[irow],temp$j[irow]] <- temp$value[irow]
      g3 <- graph_from_adjacency_matrix( abs(t(A)) , diag = FALSE)
      V(g3)$size = 20
      }
      #plot(g1,edge.arrow.size=.4)
      if (nobs>3) {
      k <- unique(nudat$id)[4]
      A <- B
      temp <- filter(nudat,id==k)
      for (irow in 1:nrow(temp))
         A[temp$i[irow],temp$j[irow]] <- temp$value[irow]
      g4 <- graph_from_adjacency_matrix( abs(t(A)) , diag = FALSE)
      V(g4)$size = 20
      }
      #plot(g1,edge.arrow.size=.4)
      
      
      par(mfrow=c(2,2), mar=c(0,0,0,0))
      plot(g1,edge.arrow.size=.4)
      plot(g2,edge.arrow.size=.4)
      if (nobs >2) plot(g3,edge.arrow.size=.4)
      if (nobs >3) plot(g4,edge.arrow.size=.4)
   })
   
   output$evalplot <- renderPlot({
      
      Wstore <- matrix(NA,nrow=7,ncol=nobs)
      Astore <- matrix(NA,nrow=7,ncol=nobs)
      if (length(unique(evaldat$id))>2) 
         nudat <- filter(evaldat,id>2)
      pick <- sample(unique(evaldat$id),size = max(5,length(unique(nudat$id))), replace=FALSE)
      nudat <- filter(evaldat, id %in% pick)
      nobs <- length(unique(nudat$id))
      Wstore <- matrix(NA,nrow=7,ncol=nobs)
      Astore <- matrix(NA,nrow=7,ncol=nobs)
      for (k in unique(nudat$id)) {
         A <- B
         temp <- filter(nudat,id==k)
         for (irow in 1:nrow(temp))
            A[temp$i[irow],temp$j[irow]] <- temp$value[irow]
         #}  
         adj_A <- adjoint3(A)
         adj_A
         Tmat <- LoopAnalyst::make.T(A,status=TRUE)
         #system.time(Tmat <- LoopAnalyst::make.T(A,status=TRUE))
         Wmat <- abs(adj_A)/Tmat
         Wmat
         
         adj_AA <- adj_A
         
         #image showing sign of adjoint and weights
         #colfunc <- colorRampPalette(c("white", "steelblue"))
         #image(1:7,1:7,t(Wmat[7:1,]),col = colfunc(7))
         #text(1,0,"+",col="white",cex=2)
         #pick <- which(adj_AA>0)
         #x <- ceiling(pick/7)
         #y <- 8-pick%%7
         #y[y==8] <- 1
         #text(x,y,"+",col="white",cex=1.5)
         #pick <- which(adj_AA<0)
         #x <- ceiling(pick/7)
         #y <- 8-pick%%7
         #y[y==8] <- 1
         #text(x,y,"-",col="white",cex=1.5)
         
         Wstore[,which(unique(nudat$id)==k)] <- rev(Wmat[,7])
         Astore[,which(unique(nudat$id)==k)] <- rev(adj_AA[,7])
         
      }
      labs <- c("Manager",
                "Habitat",
                "Fish",
                "Seabirds",
                "Fish",
                "ZP",
                "PP")
      #group_pic
      #image showing sign of adjoint and weights
      colfunc <- colorRampPalette(c("white", "steelblue"))
      par(mar=c(1,5,4,3),las=0)
      image(1:nobs,1:7,t(Wstore),col = colfunc(5),
            xlab = "",
            ylab = "", axes=F)
      box()
      #axis(1,labels=rep("",nobs),at=seq(0.5,nobs+0.5,1),tcl=-0.2)
      par(las=2)
      axis(2,labels=labs,at=seq(1,7,1),tcl=-0.2,cex=0.8)
      par(las=0)
      axis(3,labels=1:nobs,at=seq(1,nobs),tcl=-0.2,cex=0.8)
      mtext(side=3,"Model",line = 3)
      #text(1,0,"+",col="white",cex=2)
      pick <- which(Astore>0)
      x <- ceiling(pick/7)
      y <- 8-pick%%7
      y[y==8] <- 1
      text(x,y,"+",col="white",cex=1.5)
      pick <- which(Astore<0)
      x <- ceiling(pick/7)
      y <- 8-pick%%7
      y[y==8] <- 1
      text(x,y,"-",col="white",cex=1.5)
   })


}

# Run the application 
shinyApp(ui = ui, server = server)

