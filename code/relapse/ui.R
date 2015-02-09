# TODO: Add comment
# 
# Author: mg14
###############################################################################


library(shiny)
library(CoxHD)
load("predict.RData")
pdid <- "other"

# Define UI for application that plots random distributions 
shinyUI(pageWithSidebar(
				
				# Application title
				headerPanel("AML prediction"),
				
				# Sidebar with a slider input for number of observations
				do.call("sidebarPanel", 
						c(selectInput("pdid", "Select sample", c("other",rownames(coxRFXCirTD$Z)), selected = "1", multiple=FALSE),
								lapply(c("transplantRel",names(sort(apply(coxRFXCirTD$Z,2,var)*coef(coxRFXCirTD)^2, decreasing=TRUE)[-102])), function(x) numericInput(x, paste(x, " (mean=",round(mean(coxRFXCirTD$Z[,x]), 2),"; HR_CIR=",round(exp(coef(coxRFXCirTD)[x]),2),"; HR_NRM=",round(exp(coef(coxRFXNrmTD)[x]),2),"; HR_PRM=",round(exp(coef(coxRFXPrsTD)[x]),2),")", sep=""), ifelse(pdid=="other",NA,coxRFXCirTD$Z[pdid,x]), )
								))
				),
				# Show a plot of the generated distribution
				mainPanel(
						plotOutput(outputId="KM",height="600px"),
						tableOutput(outputId = "Risk"),
						tableOutput(outputId = "Tab")
				)
		))

