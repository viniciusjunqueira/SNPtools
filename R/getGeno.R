#' Leitura flexível e eficiente de genótipos com autodetecção usando fread
#'
#' @param path Caminho para o FinalReport.txt
#' @param fields Lista com as colunas: sample, snp, allele1, allele2, confidence
#' @param codes Vetor com os códigos de alelos (ex: c("A", "B"))
#' @param threshold Corte para qualidade (confidence)
#' @param sep Separador usado no arquivo
#' @param skip Linhas a pular no topo
#' @param verbose Exibir progresso?
#' @param every Frequência de progresso
#'
#' @return Objeto da classe SNPDataLong ou NULL em caso de erro
#' @export
setGeneric("getGeno", function(...) standardGeneric("getGeno"))

setMethod("getGeno", signature(),
  function(path,
           fields = list(sample = 2, snp = 1, allele1 = 7, allele2 = 8, confidence = 9),
           codes = c("A", "B"),
           threshold = 0.15,
           sep = "\t",
           skip = 0,
           verbose = TRUE,
           every = NULL) {

    if (!requireNamespace("data.table", quietly = TRUE)) {
      stop("O pacote 'data.table' é necessário. Instale com install.packages('data.table').")
    }

    if (!file.exists(file.path(path, "FinalReport.txt"))) {
      warning("Arquivo FinalReport.txt não encontrado em: ", path)
      return(NULL)
    }

    # Leitura rápida com fread
    cols_to_read <- unique(unlist(fields[c("sample", "snp")]))
    col_select <- as.integer(cols_to_read)

    fread_result <- tryCatch({
      data.table::fread(
        file = file.path(path, "FinalReport.txt"),
        sep = sep,
        skip = skip,
        header = TRUE,
        select = col_select,
        data.table = FALSE,
        colClasses = "character"
      )
    }, error = function(e) {
      warning("Falha ao ler FinalReport.txt com fread: ", e$message)
      return(NULL)
    })

    if (is.null(fread_result)) return(NULL)

    sample.col <- match(fields$sample, col_select)
    snp.col    <- match(fields$snp, col_select)

    sample.id <- unique(fread_result[[sample.col]])
    snp.id    <- unique(fread_result[[snp.col]])

    if (is.null(every)) every <- length(snp.id)

    # Leitura com snpStats
    data <- tryCatch({
      snpStats::read.snps.long(
        file = file.path(path, "FinalReport.txt"),
        sample.id = sample.id,
        snp.id = snp.id,
        fields = fields,
        codes = codes,
        threshold = threshold,
        sep = sep,
        skip = skip,
        verbose = verbose,
        every = every
      )
    }, error = function(e) {
      warning("Erro ao executar read.snps.long: ", e$message)
      return(NULL)
    })

    if (is.null(data)) return(NULL)

    # Leitura do mapa
    map_file <- file.path(path, "SNP_Map.txt")
    if (!file.exists(map_file)) {
      warning("Arquivo SNP_Map.txt não encontrado em: ", path)
      return(NULL)
    }

    map <- tryCatch({
      read.table(map_file, colClasses = "character", sep = sep, header = TRUE)
    }, error = function(e) {
      warning("Erro ao ler SNP_Map.txt: ", e$message)
      return(NULL)
    })

    if (is.null(map)) return(NULL)

    map <- map[, c("Name", "Chromosome", "Position", "GenTrain.Score")]

    new("SNPDataLong",
        geno = data,
        map  = map,
        path = path)
  }
)

#' Importa e combina múltiplas configurações de genótipos
#'
#' @param object Objeto do tipo SNPImportList
#'
#' @return Objeto SNPDataLong combinado
#' @export
setGeneric("importAllGenos", function(object) standardGeneric("importAllGenos"))

setMethod("importAllGenos", "SNPImportList", function(object) {
  all_genos <- lapply(object@configs, function(cfg) {
    tryCatch(getGeno(cfg), error = function(e) {
      warning("Erro em getGeno(): ", e$message)
      NULL
    })
  })

  all_genos <- Filter(Negate(is.null), all_genos)
  if (length(all_genos) == 0) stop("Nenhum genótipo foi importado com sucesso.")
  combinarSNPData(all_genos)
})
