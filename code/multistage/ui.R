library(shiny)
library(shinyBS)
library(CoxHD)
load("multistage.RData", envir=globalenv())


# Define UI for application that plots random distributions 
#shinybootstrap2::withBootstrap2({
fluidPage(
		tags$head(
				includeHTML("www/popup.html")
		),
		includeHTML("www/disclaimer.html"),

		# Application title
		titlePanel("AML multistage predictions (research only)"),
		#div(HTML('<h4 style="color:red;""> For research use only</h4>')),
		
		
		fluidRow(
				# Sidebar with a slider input for number of observations
				column(3, 
						wellPanel(
								actionButton("compute", "Compute survival")
								),
						wellPanel(
								tags$b("1. Select sample"),								
								selectizeInput(inputId="pdid", label="", choices=c("reset",rownames(data)[order(as.numeric(gsub("[A-z]","", rownames(data))))]),  multiple=FALSE, 
										options = list(maxOptions = nrow(data)+1,
												placeholder = 'Please select',
												onInitialize = I('function() { this.setValue(""); }'))),
								tags$em(tags$small("Data may be rounded for privacy reasons."))
								
						),
						#bsCollapse(#"2. Enter/change variables",
								uiOutput("ui"),
						#),
						wellPanel(
								radioButtons("ciType", tags$b("Confidence intervals"), choices=c("analytical (fast, CR only)"="analytical","simulated (slow)"="simulated"), selected = "analytical"), ## CI type
								#tags$hr(),
								#div(HTML('<b><a href="help.html">Help</a></b>')),
							    div(HTML('<b><a id="disclaimer">Disclaimer</a></b>')))
				),
				
				# Show a plot of the generated distribution
				column(9,
						tabsetPanel(
								tabPanel('Results',
										tags$h4("Patient summary"),
										htmlOutput(outputId="patientSummary"),
										tags$h4("Multistage probabilities"),
										plotOutput(outputId="KM",height="300px"),
										tags$h4("Outcome 3 years after diagnosis"),
										htmlOutput(outputId="absoluteRiskDiag"),
										tags$h4("Outcome 3 years after remission"),
										htmlOutput(outputId="absoluteRiskCr")),
								
								tabPanel('Log hazard',
										dataTableOutput("Risk")),
								tabPanel("Coefficients",
										dataTableOutput("Tab")),
								tabPanel("Help",
										includeHTML("www/help.html"))
						))
		
		)
)
#})
