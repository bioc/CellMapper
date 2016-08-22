.CellMapperList <- setClass("CellMapperList", contains="SimpleList")

setValidity("CellMapperList", function(object) {
    valid <- TRUE
    if (identical(sort(names(object)), sort(c("B", "d")))) {
        if (!(is.matrix(object$B) & is.numeric(object$B) & is.vector(object$d) &
		 is.numeric(object$d))) {
            valid <- "'B' must be numeric matrix and 'd' must be numeric vector"
        } else if (length(object$d) != dim(object$B)[1]) {
            valid <- "The length of 'd' must equal the number of rows in 'B'"
        }
    } else {
        valid <- "'B' and 'd' must be present"
    }
    if (!is.list(metadata(object))) { valid <- paste0("metadata must be provid",
	"ed as a list") }
    valid
})

setMethod("show", "CellMapperList",
	function(object)
	{
		meta = metadata(object)
		metaValid = if (is.null(names(meta))) { FALSE } else { 
			identical(sort(names(meta)), sort(c("nrow", "ncol",
				"DataSource", "GeneIDType"))) }
		cat("An object of class \"", class(object), "\"\n",
		"#  Provide as input to the 'CMsearch' function of the",
		" 'CellMapper' package\n",
		if (metaValid) {
			paste0("#  Derived from an expression dataset with ",
			     meta$nrow, " genes and ", meta$ncol, " samples\n",
			"#     Dataset source:  '", if (meta$DataSource == "") {
			     "not provided" } else { meta$DataSource }, "'\n",
			"#  The type of gene ID used is: '",
			if (meta$GeneIDType == "") { "unknown" } else {
				meta$GeneIDType }, "'\n",
			"#     Example gene IDs:  '",
			paste(head(colnames(object$B)), collapse = "', '"),
			"', ...\n")
		} else { "#  No metadata available\n" }, sep = '')
		if (!metaValid) {
			warning("This object was not prepared with the",
			" 'CMprep' function")
		}
	}
)

CellMapperList <- function(B, d, meta = list()) {
    .CellMapperList(listData = list(B = B, d = d), metadata = meta)
}

