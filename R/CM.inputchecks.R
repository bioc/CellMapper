## Check which query genes match gene names in the dataset
chk.Genes <-
function(query.genes, Data, type = 'query.genes', multi = FALSE) {
	if (multi) {
		if (is.character(query.genes)) {
			QGs = matrix(NA, nrow = length(Data), ncol =
				length(query.genes))
			colnames(QGs) = query.genes
			rownames(QGs) = names(Data)
			for (i in names(Data)) {
				QGs[i,] <- match(query.genes,
					colnames(Data[[i]]$B))
			}

			if (all(is.na(QGs))) {
				if (type == 'query.genes') {
					stop("None of the 'query.genes' match ",
					"any gene identifiers in the provided ",
					"datasets\n")
				} else {
					warning("None of the 'control.genes' ",
					"match any gene identifiers in the ",
					"provided datasets\n")
					QGs = NULL
				}
			} else {
				if (any(rowMeans(is.na(QGs)) == 1)) {
					warning("None of the '", type, "' ",
					"match any gene identifiers in the ", 
					"dataset(s): '",
paste(rownames(QGs)[rowMeans(is.na(QGs)) == 1], collapse = ", "),
					if (type == 'query.genes') {
"'. These datasets will be excluded from the analysis.\n" } else {
"No 'control.genes' were used for this dataset\n" })
				}

				if (any(is.na(QGs))) {
		for (D in rownames(QGs)[!(rowMeans(is.na(QGs)) %in% c(0,1))]) {
			warning("The following '", type,
			"' did not match any gene identifiers in '", D, "': '",
		paste(colnames(QGs)[is.na(QGs[D,])], collapse = "', '"), "'\n")
		}
				}
			}
		} else {
			stop("The '", type, "' must be provided as a character",
			" vector. See ?as.character\n")
		}
	} else {
		if (is.character(query.genes)) {
			QGs <- match(query.genes, colnames(Data$B))
			names(QGs) = query.genes
			if (all(is.na(QGs))) {
				if (type == 'query.genes') {
					stop("None of the 'query.genes' match ",
			"any gene identifiers in the provided dataset\n")
				} else {
					warning("None of the 'control.genes' ",
"match any gene identifiers in the provided dataset. No 'control.genes' were ",
"used for this analysis\n")
					QGs = NULL
				}
			} else if (any(is.na(QGs))) {
				warning("The following '", type, "' did not ",
"match any gene identifiers in the provided dataset: '",
paste(query.genes[is.na(QGs)], collapse = "', '"), "'\n")
				QGs <- na.omit(QGs)
			}
		} else {
			stop("The '", type, "' must be provided as a character",
" vector of gene identifiers (e.g. Entrez IDs). See ?as.character\n")
		}
	}

	return(QGs)
}

