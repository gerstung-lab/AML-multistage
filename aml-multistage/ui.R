library(shiny)
library(CoxHD)
load("multistage.RData", envir=globalenv())

gitLog <- system("git log --pretty=format:'Revision %h, commited by %an on %ai' -n 1", intern=TRUE)
wellStyle <- "background-color:rgb(255, 255, 255); border-color:rgb(204, 205, 205); padding-bottom:9px; padding-top:9px;"


# Define UI for application that plots random distributions 
#shinybootstrap2::withBootstrap2({
fluidPage(
		tags$head(
				includeHTML("www/popup.html")
		),
		includeHTML("www/disclaimer.html"),

		# Application title
		titlePanel(HTML("AML multistage predictions <small>research only</small>"), windowTitle = "AML multistage predictions"),
		#div(HTML('<h4 style="color:red;""> For research use only</h4>')),
		
		
		fluidRow(
				# Sidebar with a slider input for number of observations
				column(3, 
						wellPanel(
								tags$b("1. Select sample (or enter new patient's variables below)"),								
								selectizeInput(inputId="pdid", label="", choices=c("reset all variables",rownames(data)[order(as.numeric(gsub("[A-z]","", rownames(data))))]),  multiple=FALSE, 
										options = list(maxOptions = nrow(data)+1,
												placeholder = 'Please select',
												onInitialize = I('function() { this.setValue(""); }'))),
								tags$em(tags$small("Data may be rounded for privacy reasons."))
								
						),
						#bsCollapse(#"2. Enter/change variables",
						wellPanel(tags$b("2. Enter/amend variables"),	
												wellPanel(
														actionLink("showClinical", 					
																tags$div("Clinical variables", HTML("&#9662;")),
																style = "color:rgb(0,0,0);"
														),
														style = wellStyle),
												uiOutput("expandClinical"),
												wellPanel(
														actionLink("showDrivers", 					
																tags$div("Driver mutations", HTML("&#9662;")),
																style = "color:rgb(0,0,0);"
														),
														style = wellStyle),
												
												uiOutput("expandDrivers"),
												wellPanel(
														actionLink("showTreatment", 					
																tags$div("Treatment", HTML("&#9662;")),
																style = "color:rgb(0,0,0);"
														),
														style = paste(wellStyle, "margin-bottom: 0px")),
												uiOutput("expandTreatment")
												),

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
								),
								uiOutput("inputCheck"),
								uiOutput("multistageCheck")

						)
				),
				
				# Show a plot of the generated distribution
				column(9,
						tabsetPanel(
								tabPanel('Results',
										tags$h3("Patient summary"),
										htmlOutput(outputId="patientSummary"),
										tags$h3("Outcome after diagnosis"),
										#tags$h4("Multistage probabilities"),
										plotOutput(outputId="multistageDiag",height="300px", width="500px"),
										tags$h4("Outcome 3 years after diagnosis"),
										htmlOutput(outputId="absoluteRiskDiag"),
										tags$h3("Outcome after first complete remission"),
										#tags$h4("Multistage probabilities"),
										plotOutput(outputId="multistageCR",height="300px", width="500px"),
										tags$h4("Outcome 3 years after remission"),
										htmlOutput(outputId="absoluteRiskCr")),
								
								tabPanel('Log hazard',
										dataTableOutput("Risk")),
								tabPanel("Coefficients",
										dataTableOutput("Tab")),
								tabPanel("Help",
										includeHTML("www/help.html"),
										tags$h4("Version info"),
										tags$code(gitLog))
						))
		
		)
)
#})
