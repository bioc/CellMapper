CMsearch <-
function(Data, query.genes = NULL, control.genes = NULL, QDW = TRUE,
alpha = 0.5, verbose = TRUE, raw.pvals = FALSE) {
	if (class(Data) == "CellMapperList") {
		query.genes <- chk.Genes(query.genes, Data)
		if (!is.null(control.genes)) { control.genes <-
			chk.Genes(control.genes, Data, type = 'control.genes') }

		if (verbose) {
		message("* Running CellMapper with query gene(s) '",
		paste(names(query.genes), collapse = "', '"),
		if (!is.null(control.genes)) {
			paste("' and control genes '",
			paste(names(control.genes), collapse = "', '"),
			sep = "") } else { "" }, "'...\n")
		}

		out <- CM.search.sub(Data, QGs = query.genes,
CGs = control.genes, QDW = QDW, alpha = alpha, verbose = verbose)

	} else if (is.list(Data)) {
		if (all(sapply(Data, function(x) {
		class(x) == "CellMapperList" }))) {
			if (is.null(names(Data))) {
				names(Data) <- paste("Dataset",
				seq_len(length(Data)))
			} else if (any(names(Data) == "")) {
				names(Data)[names(Data) == ""] <-
				paste("Dataset", which(names(Data) == ""))
			}

			query.genes <-
				chk.Genes(query.genes, Data, multi = TRUE)
			if (!is.null(control.genes)) {
				control.genes <-
				chk.Genes(control.genes, Data,
					type = 'control.genes', multi = TRUE)
			}

			out <- lapply(
rownames(query.genes)[rowMeans(is.na(query.genes)) < 1], function(D) {
				if (!is.null(control.genes)) {
					if (all(is.na(control.genes[D,]))) {
					controls <- NULL
					} else {
					controls <- na.omit(control.genes[D,])
					}
				} else { controls <- NULL }

				if (verbose) {
					message("* Running CellMapper on '", D,
					"' with query gene(s) '",
paste(colnames(query.genes)[!is.na(query.genes[D,])], collapse = "', '"),
					if (!is.null(controls)) {
paste("' and control genes '",
paste(colnames(control.genes)[!is.na(control.genes[D,])], collapse = "', '"),
sep = "")
					} else { "" }, "'...\n")
				}

				return(CM.search.sub(Data[[D]], QGs =
				na.omit(query.genes[D,]), CGs = controls,
				QDW = QDW, alpha = alpha, verbose = verbose))
			})

			if (verbose) { message("* Pooling results...\n") }
			out <- CM.poolResults(out)
		} else {
			stop("Please provide a dataset that has been ",
			"preprocessed with the 'CMprep' function, or a list of",
			" pre-processed datasets\n")
		}
	} else {
		stop("Please provide a dataset that has been preprocessed with",
		" the 'CMprep' function, or a list of pre-processed datasets\n")
	}

	out <- data.frame(Gene = names(out),
		FDR = p.adjust(exp(out), method = 'BH'),
		p_unadjusted = exp(out))
	rownames(out) = seq_len(dim(out)[1])
	if (!raw.pvals) { out <- out[,-3] }

	return(out)
}

