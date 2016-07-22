CM.prep <-
function(Data, verbose = TRUE) {

	if (is.matrix(Data) & is.numeric(Data)) {
		dataMat <- Data
		rm(Data)
	} else if (class(Data) == "ExpressionSet") {
		if (requireNamespace("Biobase")) {
			dataMat <- Biobase::exprs(Data)
			rm(Data)
		} else { stop("'Biobase' must be installed to access an 'ExpressionSet' object")
		}
	} else {
		stop("'Data' must be either an Affy 'ExpressionSet' object or a numeric matrix\n")
	}

	if (any(is.na(dataMat))) {
		stop("There are NAs in the input data. At this point, CellMapper cannot accept datasets with missing values. Consider using a packages such as 'impute' to estimate the missing expression values.")
	}


	if (verbose) {
		message("* Scaling the data matrix...\n")
	}
	dataMat <- sweep(dataMat, 1, apply(dataMat, 1, mean), "-")
	dataMat <- sweep(dataMat, 1, sqrt(apply(dataMat, 1, var)), "/")
		

	if (verbose) {
		message("* Computing the singular-value decomposition of the data matrix...\n")
	}
	sv <- svd(dataMat)


	if (verbose) {
		message("* Preparing data for CellMapper...\n")
	}
	keep <- sum(length(sv$d)*sv$d/sum(sv$d) > 1)
	sv$u <- sv$u[,seq_len(keep)]
	sv$d = sv$d[seq_len(keep)]
	rownames(sv$u) <- rownames(dataMat)

	out <- list(B = t(sv$u), d = sv$d)
	return(out)
}
