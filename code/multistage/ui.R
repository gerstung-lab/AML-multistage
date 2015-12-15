library(shiny)
library(CoxHD)
load("multistage.RData", envir=globalenv())


# Define UI for application that plots random distributions 

fluidPage(
		# Style
		tags$header(tags$style(
						type = 'text/css',
						'.well .special { max-height: 400px; overflow-y: auto; }'
				)
		),

		# Application title
		titlePanel("AML multistage predictions (beta)"),
		
		fluidRow(
				# Sidebar with a slider input for number of observations
				column(3, 
						wellPanel(
								selectInput("pdid", tags$b("Select sample"), c("reset",rownames(data)), selected = "reset", multiple=FALSE)
						),
						wellPanel(								
								tags$b("Prognostic variables"),
								tags$hr(),
								HTML('<div class="special">'),
								uiOutput("ui"),
								HTML("</div>")
						),
						wellPanel(actionButton("compute", "Compute survival"),
								tags$hr(),
								radioButtons("ciType", tags$b("Confidence intervals"), choices=c("analytical (fast, CR only)"="analytical","simulated (slow)"="simulated"), selected = "analytical"), ## CI type
								tags$hr(),
								div(HTML('<b><a href="help.html">Help</a></b>')))
				),
				
				# Show a plot of the generated distribution
				column(8,
						plotOutput(outputId="KM",height="300px"),
						tabsetPanel(
								tabPanel('Risk',
										dataTableOutput("Risk")),
								tabPanel("Coefficients",
										dataTableOutput("Tab"))
						))
		
		)
)

