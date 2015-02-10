# TODO: Add comment
# 
# Author: mg14
###############################################################################


library(shiny)
library(CoxHD)
load("predict.RData")
data <- coxRFXCirTD$Z

# Define UI for application that plots random distributions 
shinyUI(fluidPage(
				
				# Application title
			    titlePanel("AML prediction tool"),
				
				fluidRow(
							# Sidebar with a slider input for number of observations
							column(3, 
									wellPanel(
												selectInput("pdid", "Select sample", c("reset",rownames(data)), selected = "1", multiple=FALSE) ## select sample
											),
									wellPanel(
												uiOutput("ui")
											)
							),
							
							# Show a plot of the generated distribution
							column(8,
									plotOutput(outputId="KM",height="600px"),
									tabsetPanel(
											tabPanel('Risk',
													dataTableOutput("Risk")),
											tabPanel("Coefficients",
													dataTableOutput("Tab"))
									))

						)
		))

