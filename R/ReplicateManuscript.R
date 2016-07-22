ReplicateManuscript <-
function() {
	if (!requireNamespace("CellMapperData")) {
		stop("The 'CellMapperData' package must be installed before using this function")
	} else { 
		# Load the datasets and query genes
		message("* Loading the pre-processed datasets...")
		data(ExampleQueryGenes, envir = environment())
		QGs <- get("ExampleQueryGenes", envir = environment())
		data(list = unique(unlist(strsplit(QGs$Dataset, ", "))), package = "CellMapperData", envir = environment())

		# Run CellMapper for each cell type
		environ <- environment()
		Results <- lapply(seq_len(dim(QGs)[1]), function(i) {
			message("* Running CellMapper for ", QGs$CellType[i], "...")

			dataset <- mget(strsplit(QGs$Dataset[i], ", ")[[1]], envir = environ)
			if (length(dataset) == 1) { dataset <- dataset[[1]] }

			return(suppressWarnings(CM.search(dataset, query.genes = QGs$EntrezID[i], verbose = FALSE)))
		})
		names(Results) <- QGs$CellType

		return(Results)
	}
}
