# TODO: Add comment
# 
# Author: mg14
###############################################################################


library(shiny)
library(CoxHD)
load("predict.RData")


# Define UI for application that plots random distributions 
shinyUI(pageWithSidebar(
				
				# Application title
				headerPanel("AML prediction"),
				
				# Sidebar with a slider input for number of observations
				do.call("sidebarPanel", lapply(c("transplantRel",names(sort(apply(coxRFXCirTD$Z,2,var)*coef(coxRFXCirTD)^2, decreasing=TRUE)[-102])), function(x) numericInput(x, paste(x, " (mean=",round(mean(coxRFXCirTD$Z[,x]), 3),"; HR=",round(exp(coef(coxRFXCirTD)[x]),3),")", sep=""), NA, ))
				),
				# Show a plot of the generated distribution
				mainPanel(
						plotOutput(outputId="KM",height="700px"),
						tableOutput(outputId = "Risk"),
						tableOutput(outputId = "Tab")
				)
		))

