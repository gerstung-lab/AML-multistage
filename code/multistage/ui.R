library(shiny)
library(shinyBS)
library(CoxHD)
load("multistage.RData", envir=globalenv())

wellStyle <- "background-color:rgb(255, 255, 255); border-color:rgb(204, 205, 205); padding-bottom:9px; padding-top:9px;"


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
								tags$b("1. Select sample"),								
								selectizeInput(inputId="pdid", label="", choices=c("reset all variables",rownames(data)[order(as.numeric(gsub("[A-z]","", rownames(data))))]),  multiple=FALSE, 
										options = list(maxOptions = nrow(data)+1,
												placeholder = 'Please select',
												onInitialize = I('function() { this.setValue(""); }'))),
								tags$em(tags$small("Data may be rounded for privacy reasons."))
								
						),
						#bsCollapse(#"2. Enter/change variables",
						
								uiOutput("ui"),
						#),
						wellPanel(
								actionButton("compute", tags$b("3. Compute outcome"), class="btn btn-primary", style = "margin-bottom:20px"),
								tags$br(),
								wellPanel(
										actionLink("showOptions", 					
												tags$div("Options", HTML("&#9662;")),
												style = "color:rgb(0,0,0);"
										),
										conditionalPanel(condition = 'input.showOptions % 2',
												#tags$hr(),
												radioButtons("ciType", tags$b("Confidence intervals"), choices=c("analytical (fast, CR only)"="analytical","simulated (slow)"="simulated"), selected = "analytical"), ## CI type
												style = "overflow-y:scroll; max-height: 400px; position:relative"
										),
										style = paste(wellStyle, "margin-bottom:0px")
								)
						)
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
