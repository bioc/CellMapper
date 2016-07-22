CM.search.sub <-
function(Data, QGs = NULL, CGs = NULL, QDW = TRUE, alpha = 0.5, verbose = TRUE) {
	if (length(QGs) > 1) {
		QGcors <- t(Data$B[,QGs]*Data$d^2) %*% Data$B[,QGs]
		diags <- sqrt(colSums((Data$B[,QGs]*Data$d)^2))
		QGcors <- sweep(QGcors/diags, 2, diags, '/')
		QG0 <- names(which(apply(QGcors < 0, 1, any)))
		QG5 <- names(which(apply(QGcors < 0.8, 1, any)))
		QG5 <- QG5[!(QG5 %in% QG0)]
		if (length(QG0) > 0) {
			warning("The following 'query.genes' are negatively correlated with other 'query.genes': '", paste(QG0, collapse = "', '"), "'\n")
		}
		if (length(QG5) > 0) {
			warning("The following 'query.genes' are weakly correlated (< 0.5) with other 'query.genes': '", paste(QG5, collapse = "', '"), "'\n")
		}
	}  # Warn the user if multiple query genes are used that are poorly correlated

	w <- Data$d^alpha
	if (QDW) {
		if (length(QGs) > 1) {
			w2 <- sweep(sweep(Data$B[,QGs], 1, rowMeans(Data$B), '-'), 1, apply(Data$B, 1, sd), '/')
			w2 <- rowMeans(tanh(w2))
		} else if (length(QGs) == 1) {
			w2 <- (Data$B[,QGs] - rowMeans(Data$B))/apply(Data$B, 1, sd)
			w2 <- tanh(w2)
		}
		w <- abs(w2)*w
	}


	if ((length(QGs) + length(CGs)) > 1) {
		cov.mat <- t(sweep(Data$B[, c(QGs, CGs)], 1, w^2, '*')) %*% Data$B
	} else {
		cov.mat <- t(Data$B[, QGs]*w^2) %*% Data$B
	}

	diags <- sqrt(colSums((Data$B*w)^2))
	cor.mat <- cov.mat/diags[c(QGs, CGs)]
	cor.mat <- sweep(cor.mat, 2, diags, '/')

	colnames(cor.mat) <- colnames(Data$B)
	rownames(cor.mat) <- colnames(Data$B)[c(QGs, CGs)]

	if (!is.null(CGs)) { # Adjust for positive partial correlation with 'control.genes'
		cor.mat <- adjust.pcor(cor.mat, QGs, CGs, verbose)
	}

	if (length(QGs) > 1) {
		cor.vec <- colMeans(atanh(cor.mat[seq_len(length(QGs)),-c(QGs, CGs)]))
	} else if (!is.null(dim(cor.mat))) {
		cor.vec <- atanh(cor.mat[1,-c(QGs, CGs)])
	} else {
		cor.vec <- atanh(cor.mat[-c(QGs, CGs)])
	}

	S <- (cor.vec - median(cor.vec, na.rm = TRUE)) / mad(cor.vec, na.rm = TRUE)

	probs <- pnorm(S, lower.tail = FALSE, log.p = TRUE)
	probs <- probs[order(probs)]

	return(probs)
}

adjust.pcor <-
function(cor.mat, QGs, CGs, verbose = TRUE) {
	if (length(QGs) > 1) {
		CG2 <- which(rownames(cor.mat) %in% names(which(colMeans(cor.mat[seq_len(length(QGs)),CGs] > 0)==1)))
		Positive.cors <- colMeans(sign(cor.mat[seq_len(length(QGs)),])) == 1
	} else {
		CG2 <- which(rownames(cor.mat) %in% names(which(cor.mat[1,CGs] > 0)))
		Positive.cors <- cor.mat[1,] > 0
	}  # Only consider 'control genes' that are positively correlated with all 'query genes', and only calculate partial correlation for genes that are positively correlated with all 'query genes'

	cor.mat <- cor.mat[c(seq_len(length(QGs)), CG2),]

	if (length(CG2) > 0) {
		if (verbose) { message("     Adusting for partial correlation with gene(s): '", paste(rownames(cor.mat)[-c(seq_len(length(QGs)))], collapse = "', '"), "'\n") }

		cor.mat.rd <- cor.mat
		for (z in dim(cor.mat.rd)[1]:(length(QGs)+1)) {
			p.zy <- cor.mat.rd[z,Positive.cors]
			for (x in seq_len(z-1)) {
				p.xz <- cor.mat.rd[x, rownames(cor.mat)[z]]
				denom <- (1 - p.xz^2)*(1 - p.zy^2)
				denom[denom < 10^-12] <- NA
				cor.mat.rd[x,Positive.cors] <- (cor.mat.rd[x,Positive.cors] - p.xz*p.zy)/sqrt(denom)
			}
		}
		cor.mat <- sweep(cor.mat.rd, 1, diag(cor.mat.rd[,rownames(cor.mat.rd)]), '/')
	}
	return(cor.mat)
}

CM.poolResults <-
function(CM.output) {
	CM.output <- CM.output[!sapply(CM.output, is.null)]

	if (length(CM.output) == 1) {
		out <- CM.output[[1]]
	} else {

		genes <- unique(unlist(lapply(seq_len(length(CM.output)), function(i) { names(CM.output[[i]]) })))

		p.vals <- sapply(CM.output, "[", genes)
		out <- apply(p.vals, 1, combine.stouffer)
		out <- out[order(out)]
	}

	return(out)
}

combine.stouffer <-
function(p.vals) {
	p.vals <- p.vals[!is.na(p.vals)]
	if (length(p.vals) > 0) {
		p.stouffer <- pnorm(sum(qnorm(p.vals, log.p = TRUE))/sqrt(length(p.vals)), log.p = TRUE)
		return(p.stouffer)
	} else { return(NA) }
}
