library(shiny)
library(CoxHD)
load("multistage.RData", envir=globalenv())


# Define UI for application that plots random distributions 
fluidPage(
		tags$head(
				includeHTML("www/popup.html")
		),
		includeHTML("www/disclaimer.html"),
		# Style
		tags$header(tags$style(
						type = 'text/css',
						'.well .special { max-height: 400px; overflow-y: auto; }'
				)
		),

		# Application title
		titlePanel("AML multistage predictions (beta)"),
		div(HTML('<h4 style="color:red;""> For research use only</h4>')),
		
		
		fluidRow(
				# Sidebar with a slider input for number of observations
				column(3, 
						wellPanel(
								tags$b("Select sample"),
								tags$em(tags$small("Data may be rounded for privacy reasons.")),								
								selectizeInput(inputId="pdid", label="", choices=c("reset",rownames(data)[order(as.numeric(gsub("[A-z]","", rownames(data))))]),  multiple=FALSE, 
										options = list(maxOptions = nrow(data)+1,
												placeholder = 'Please select',
												onInitialize = I('function() { this.setValue(""); }'))),
								#tags$hr(),
								tags$br(),
								actionButton("compute", "Compute survival")
						),
								uiOutput("ui"),
						wellPanel(
								radioButtons("ciType", tags$b("Confidence intervals"), choices=c("analytical (fast, CR only)"="analytical","simulated (slow)"="simulated"), selected = "analytical"), ## CI type
								#tags$hr(),
								#div(HTML('<b><a href="help.html">Help</a></b>')),
							    div(HTML('<b><a id="disclaimer">Disclaimer</a></b>')))
				),
				
				# Show a plot of the generated distribution
				column(8,
						tabsetPanel(
								tabPanel('Results',
										tags$h4("Patient summary"),
										textOutput(outputId="patientSummary", container=pre),
										tags$h4("Multistage probabilities"),
										plotOutput(outputId="KM",height="300px"),
										tags$h4("3-year post diagnosis risk estimates"),
										textOutput(outputId="absoluteRiskDiag", container=pre),
										tags$h4("3-year post CR risk estimates"),
										textOutput(outputId="absoluteRiskCr", container=pre)),
								
								tabPanel('Log hazard',
										dataTableOutput("Risk")),
								tabPanel("Coefficients",
										dataTableOutput("Tab")),
								tabPanel("Help",
										includeHTML("www/help.html"))
						))
		
		)
)

