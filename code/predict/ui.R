# TODO: Add comment
# 
# Author: mg14
###############################################################################


library(shiny)
library(CoxHD)
#load("predictGG.RData", envir=globalenv())
#data <- coxRFXCirTD$Z

# Define UI for application that plots random distributions 
shinyUI(fluidPage(
				
				# Application title
				titlePanel("AML prediction tool"),
				
				fluidRow(
						# Sidebar with a slider input for number of observations
						column(3, 
								wellPanel(
										selectInput("pdid", tags$b("Select sample"), c("reset",rownames(data)), selected = "reset", multiple=FALSE), ## select sample
										submitButton("Load presets")
								),
								wellPanel(
										tags$b("Prognostic variables"),
										submitButton("Compute survival"),
										uiOutput("ui")
								),
								wellPanel(
										checkboxGroupInput("ciType", tags$b("Confidence intervals"), c("analytical","simulated"), selected = "analytical") ## CI type
								)
						),
						
						# Show a plot of the generated distribution
						column(8,
								plotOutput(outputId="KM",height="800px"),
								tabsetPanel(
										tabPanel('Risk',
												dataTableOutput("Risk")),
										tabPanel("Coefficients",
												dataTableOutput("Tab"))
								))
				
				)
		))

