#Clear away old session chem.physical_and_invitro.data if ran on server
if("httk" %in% (.packages())){
  detach("package:httk", unload=TRUE)
}

#Need libraries installed
#install.packages(c("ChemmineF", "remotes", "httk", "openxlsx", "rcdk", "shiny")
#remotes::install_github("miraisolutions/godmode")

#Need to have java with JAVA_HOME defined
#Can be in Renvironment.site file
#Can be done manually as well. Example: >Sys.setenv(JAVA_HOME= "C:/Users/MBEAL/Software/jdk-13.0.1")
#Required for rcdk to work (required for Pubchem and Morgan Fingerprints)

#Ensure that the current directory is the same as the app.R file
#Example setwd("c:/users/mbeal/desktop/Bioactivity_app")

#Copy and paste the all of the code below to run the application.

library(ChemmineR)
library(godmode) #required to overwrite chem.physical_and_invitro.data on server version
library(httk)
library(openxlsx)
library(rcdk)
library(shiny)

ui <- fluidPage(

	tags$title("ToxCast Bioactivity and Read-Across Tool"),
			
	hr(),
	
	sidebarPanel(

		tags$h3("HTTK Data Input"),
		
		fileInput("file1", "Upload file containing CASRN, SMILES, and HTTK data",
			multiple = FALSE,
			accept = c(".xlsx")),
		
		helpText("Click on the link below to download templates with example data. "),
		
		downloadLink("templateFile", "Download input template"),
		
		radioButtons("overwriteOption", "Overwrite existing HTTK data if already present?", c("Yes", "No")),
		
		tags$br(),
		
		tags$h3("Read-Across"),
		
		radioButtons("fpType", "Choose fingerprint type for read-across", c("Morgan (ECFP6)", "Pubchem", "Other")),
		
		uiOutput('ui.CustomFingerprintTarget'),
		
		uiOutput('ui.CustomFingerprintTargetExample'),
		
		uiOutput('ui.CustomFingerprintTargetHelp'),
	
		uiOutput('ui.CustomFingerprintAnalog'),
		
		uiOutput('ui.CustomFingerprintAnalogExample'),
		
		uiOutput('ui.CustomFingerprintAnalogHelp'),
		
		sliderInput("s-value", "Choose structural similarity cut-off (i.e., Tanimoto coefficient or Jaccard Index):", 
			value=0.3, min=0.05, max=1, step=0.05),
		
		numericInput("k-value", "Choose maximum number of analogs for each target", value=10),
		
		numericInput("numActiveAssays", "Choose minimum number of active assays required for analogs", value=6),
		
		numericInput("numActiveBits", "Choose minimum number of active bits required for analogs and targets", value=6),
		
		uiOutput('ui.action')
	
	),
	
	mainPanel(
		
			
	)
)

server <- function(input, output) {

	options(shiny.maxRequestSize=30*1024^2, shiny.trace=TRUE)

	output$templateFile <- downloadHandler(
		filename = function() {
		"Input_Example_Template.xlsx"
		},
		content = function(file) {
			write.xlsx(read.xlsx("data/Input_Example_Template.xlsx"), file)
		}
	)

	output$ui.CustomFingerprintTarget <- renderUI({
		if (input$fpType != "Other") return()
		
		fileInput("fpCustomTargets", "Upload fingerprints file for targets",
		multiple = FALSE,
		accept = c(".xlsx"))
	})

	output$ui.CustomFingerprintTargetExample <- renderUI({
		if (input$fpType != "Other") return()
			
		downloadLink("templateTarget", "Download Target Data Template Example")
	})	

	output$templateTarget <- downloadHandler(
		filename = function() {
		"Input_FP_Example_ToxPrint.xlsx"
		},
		content = function(file) {
			write.xlsx(read.xlsx("data/Input_FP_Example_ToxPrint.xlsx"), file)
		}
	)

	output$ui.CustomFingerprintTargetHelp <- renderUI({
		if (input$fpType != "Other") return()
			
		helpText("The CASRN in the fingerprint file should match CASRN provided in HTTK data input. Columns starting at 3 and up should contain the bitwise comparisons for the fingerprint. Example provided contains the 729 ToxPrint fingerprints in column 3 to column ABC.")
	})		
	
	output$ui.CustomFingerprintAnalog <- renderUI({
		if (input$fpType != "Other") return()
			
		fileInput("fpCustomAnalogs", "Upload fingerprints file for analogs",
		multiple = FALSE,
		accept = c(".xlsx"))
	})
	
	output$ui.CustomFingerprintAnalogExample <- renderUI({
		if (input$fpType != "Other") return()
			downloadLink("templateAnalog", "Download Analog Data Template Example")
	})
	
	output$templateAnalog <- downloadHandler(
		filename = function() {
		"ToxCast_FP_Example_ToxPrint.xlsx"
		},
		content = function(file) {
			write.xlsx(read.xlsx("data/ToxCast_FP_Example_ToxPrint.xlsx"), file)
		}
	)

	output$ui.CustomFingerprintAnalogHelp <- renderUI({
		if (input$fpType != "Other") return()
			
		helpText("Example file download is slow, one can find ToxCast_FP_Example_ToxPrint.xlsx file in data folder instead. The CHID and SMILES columns should not be altered. SMILES provided can be used to generate other fingerprints.  Columns starting at 3 and up should contain the bitwise comparisons for the fingerprint. Example provided contains the 729 ToxPrint fingerprints in column 3 to column ABC.")
	})

	output$ui.action <- renderUI({
		if (is.null(input$file1)) return()
		
		if (input$fpType == "Other") {
			if (is.null(input$fpCustomTargets) | is.null(input$fpCustomAnalogs)) return()
		}
		
		downloadButton("analyzeData", "Analyze Data and Download Results")
	})

	##########################################
	#####Analyze data and generate report#####
	##########################################
	output$analyzeData <- downloadHandler(
    	filename = function() {
      	"AEDs.xlsx"
    	},
    	content = function(file) {
			#Read in data
			inFile = input$file1
			inData = read.xlsx(inFile$datapath, colNames=TRUE)
			
			#Import data into HTTK, overwrite existing data if chosen
			if (input$overwriteOption == "Yes") {
				overwriteBoolean <- TRUE
			} else {
				overwriteBoolean <- FALSE
			}
			
			chem.physical_and_invitro.data2 <- chem.physical_and_invitro.data
			chem.physical_and_invitro.data2 <- add_chemtable(inData,current.table=chem.physical_and_invitro.data2, data.list=list(Compound="Compound", CAS="CASRN", Funbound.plasma="Human.Funbound.plasma"), species="Human", reference="User Provided", overwrite=overwriteBoolean)
			chem.physical_and_invitro.data2 <- add_chemtable(inData,current.table=chem.physical_and_invitro.data2, data.list=list(Compound="Compound", CAS="CASRN", Clint="Human.Clint"), species="Human", reference="User Provided",overwrite=overwriteBoolean)
			chem.physical_and_invitro.data2 <- add_chemtable(inData,current.table=chem.physical_and_invitro.data2, data.list=list(Compound="Compound", CAS="CASRN", logp="logP"), reference="User Provided", overwrite=overwriteBoolean)
			chem.physical_and_invitro.data2 <- add_chemtable(inData,current.table=chem.physical_and_invitro.data2, data.list=list(Compound="Compound", CAS="CASRN", MW="MW"), reference="User Provided", overwrite=overwriteBoolean)
			
			#In the HTTK environment overwrite chem.physical_and_invitro.data with the new data
			originalData <- chem.physical_and_invitro.data
			godmode:::assignAnywhere("chem.physical_and_invitro.data", chem.physical_and_invitro.data2)
			
			##Calculate Css
			Css <- vector()

			for (i in inData$CASRN) {

			  tempCss = "NA" #reset tempCSS each time
			  
			  if (i %in% get_cheminfo("CAS")) { #checks if Css can be modelled
				tempCss = calc_mc_css(chem.cas=i,
				  which.quantile=c(0.95), model="3compartmentss", output.units="uM", species="Human", 
				  suppress.messages=FALSE,
				  return.samples=FALSE)
			  
				if(tempCss < 0.1) { tempCss <- 0.1 }
			  }
			  
			  Css <- c(Css, tempCss)
			}

			inData <- cbind(inData, Css)
			
			##Import bioactivity data
			toxCastData <- read.xlsx("data/Bioactivity_data.xlsx")
			
			##Get bioactivity data for compounds in ToxCast
			bioactivityConcentration <- rep(NA, nrow(inData))
			bioactivitySource <- rep(NA, nrow(inData))
			inData <- cbind(inData, bioactivityConcentration, bioactivitySource)
			
			rownames(inData) <- inData$CASRN
			rownames(toxCastData) <- toxCastData$CASRN
			overlappingCASRN <- inData$CASRN[inData$CASRN %in% toxCastData$CASRN]
			
			inData[overlappingCASRN,"bioactivityConcentration"] <- toxCastData[overlappingCASRN,"modl_ga_0.05_uM"]
			inData[overlappingCASRN,"bioactivitySource"] <- "ToxCast"
			
			##Use read-across to get bioactivity concentrations for compounds outside ToxCast
			
			###Prepare fingerprints
			if (input$fpType == "Morgan (ECFP6)") {
				toxCastFP <- read.xlsx("data/toxcastECFP.xlsx", colNames=TRUE)
				rownames(toxCastFP) <- toxCastFP[,1]
				toxCastFP <- toxCastFP[,-c(1)]

				inFP <- matrix(rep( 0, len=1024*nrow(inData)), nrow=nrow(inData))
				rownames(inFP) <- inData$CASRN

				for (i in 1:length(inData$SMILES)) {
					inFP[i,attributes(get.fingerprint(parse.smiles(inData$SMILES[i])[[1]], type="circular"))$bits] <- 1
				}
			}
			if (input$fpType == "Pubchem") {
				toxCastFP <- read.xlsx("data/toxcastPubChemFP.xlsx", colNames=TRUE)
				rownames(toxCastFP) <- toxCastFP[,1]
				toxCastFP <- toxCastFP[,-c(1)]

				inFP <- matrix(rep( 0, len=881*nrow(inData)), nrow=nrow(inData))
				rownames(inFP) <- inData$CASRN

				for (i in 1:length(inData$SMILES)) {
					inFP[i,attributes(get.fingerprint(parse.smiles(inData$SMILES[i])[[1]], type="pubchem"))$bits] <- 1
				}

			}
			if (input$fpType == "Other") {
				inFileTargets <- input$fpCustomTargets
				inFP <- read.xlsx(inFileTargets$datapath, colNames=TRUE)
				rownames(inFP) <- inFP$CASRN
				inFP <- inFP[,-which(colnames(inFP) %in% c("CASRN", "SMILES"))]
				
				inFileAnalogs = input$fpCustomAnalogs
				toxCastFP <- read.xlsx(inFileAnalogs$datapath, colNames=TRUE)
				rownames(toxCastFP) <- toxCastFP$CHID
				toxCastFP <- toxCastFP[,-which(colnames(toxCastFP) %in% c("CHID", "SMILES"))]
			}
			
			###Retain only chemicals with >= user defined active bit threshold
			toxCastFP <- toxCastFP[rowSums(toxCastFP) >= input$numActiveBits,]
			inFP <- inFP[rowSums(inFP) >= input$numActiveBits,,drop=F]
			
			###Retain only ToxCast chemicals with >= user defined active assays threshold
			toxCastFP <- toxCastFP[rownames(toxCastFP) %in% toxCastData$chid[as.numeric(toxCastData$hitc) >= input$numActiveAssays],]
			
			###Eliminate targets with bioactivity data from ToxCast
			inFP <- inFP[!rownames(inFP) %in% toxCastData$CASRN,,drop=F]
			
			###Check if there are any compounds where read-across can be applied, if not skip
			if(nrow(inFP != 0)) {
			
				###Convert fingerprints to FPset for comparison
				toxCastFPset <- as(as.matrix(toxCastFP), "FPset")
				inFPset <- as(as.matrix(inFP), "FPset")
			
				###Perform generalized read-across
				referenceChemical <- vector()
				genRA <- vector()
				similarChemicals <- vector()
				similarScores <- vector()
				bioactivityValues <- vector()
				
				rownames(toxCastData) <- toxCastData$chid
				
				for (i in 1:nrow(inFP)) {
					temp = NA
					temp = fpSim(inFPset[i], toxCastFPset, sorted=TRUE, method="Tanimoto", addone=0)
					temp = temp[temp>=input$"s-value"]

					lengthTemp <- length(temp)
					if (lengthTemp > input$"k-value") { lengthTemp <- input$"k-value" }
					temp = temp[1:lengthTemp]
					
					referenceChemical <- c(referenceChemical, rownames(inFP)[i])
	  
					genRA <- c(genRA, 10^as.numeric(sum(as.numeric(temp)*log10(toxCastData[names(temp),"modl_ga_0.05_uM"]))/sum(as.numeric(temp))))
				}
				
				##Add results to summary table
				inData[referenceChemical,"bioactivityConcentration"] <- genRA
				inData[referenceChemical,"bioactivitySource"] <- "Read-Across"
			
			}
			
			##Calculate AEDs
			AED <- as.numeric(inData$bioactivityConcentration)/as.numeric(inData$Css)
			inData <- cbind(inData, AED)
			
			##Write results
			write.xlsx(inData, file)			
			
		}
	)
}

shinyApp(ui, server)