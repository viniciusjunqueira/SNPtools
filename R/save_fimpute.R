#' Método S4 para salvar arquivos no formato FImpute
#'
#' @param object Objeto da classe FImputeExport
#' @export
setGeneric("saveFImpute", function(object, ...) standardGeneric("saveFImpute"))

#' @export
setMethod("saveFImpute", "FImputeExport", function(object, ...) {
  save_fimpute_raw(object@geno, object@map, object@path, object@name)
})

#' Exporta genótipos e mapa no formato FImpute a partir de argumentos simples
#'
#' @param geno Objeto do tipo SnpMatrix (do pacote snpStats)
#' @param map Data frame com colunas 'Name', 'Chromosome' e 'Position'
#' @param path Caminho onde os arquivos serão salvos
#' @param name Nome base para os arquivos de saída
#'
#' @export
saveFImputeRaw <- function(geno, map, path, name) {
  export <- new("FImputeExport", geno = geno, map = map, path = path, name = name)
  saveFImpute(export)
}

#' Função interna: salva os arquivos .gen e .map no formato FImpute
#'
#' @noRd
save_fimpute_raw <- function(genotype, map, caminho, name) {
  if (!requireNamespace("snpStats", quietly = TRUE)) {
    stop("O pacote 'snpStats' é necessário. Instale com install.packages('snpStats').")
  }

  if (!inherits(genotype, "SnpMatrix")) {
    stop("O objeto 'genotype' deve ser da classe 'SnpMatrix'.")
  }

  if (!is.data.frame(map)) {
    stop("'map' deve ser um data.frame.")
  }

  required_cols <- c("Name", "Chromosome", "Position")
  missing_cols <- setdiff(required_cols, colnames(map))
  if (length(missing_cols) > 0) {
    stop("O mapa está faltando as colunas obrigatórias: ", paste(missing_cols, collapse = ", "))
  }

  if (!is.character(caminho) || length(caminho) != 1) {
    stop("'caminho' deve ser uma string indicando o diretório de saída.")
  }

  if (!is.character(name) || length(name) != 1) {
    stop("'name' deve ser uma string para o nome base dos arquivos.")
  }

  if (!dir.exists(caminho)) {
    dir.create(caminho, recursive = TRUE)
  }

  # --- Arquivo .gen ---
  gen_file <- file.path(caminho, paste0(name, ".gen"))
  con <- file(gen_file, "wt")

  smp <- rownames(genotype)
  if (is.null(smp)) stop("O objeto 'genotype' deve ter nomes de linha (indivíduos).")

  nC <- max(nchar(smp), na.rm = TRUE)
  writeLines("ID    Chip                   Call...", con)

  for (i in seq_len(nrow(genotype))) {
    temp <- as(genotype[i, ], "character")
    temp[temp == "A/A"] <- "0"
    temp[temp == "A/B"] <- "1"
    temp[temp == "B/B"] <- "2"
    temp[is.na(temp) | temp == "NA"] <- "5"
    linha <- paste(format(smp[i], width = nC, justify = "left"), "1", paste(temp, collapse = ""))
    writeLines(linha, con)
  }

  close(con)

  # --- Arquivo .map ---
  map_file <- file.path(caminho, paste0(name, ".map"))
  write("SNP_ID           Chr      Pos   Chip1", file = map_file)

  map <- map[map$Name %in% colnames(genotype), ]
  if (nrow(map) == 0) {
    warning("Nenhum SNP do mapa corresponde às colunas do genótipo.")
  }

  map_out <- data.frame(map[, required_cols, drop = FALSE], 1:nrow(map))
  write.table(
    map_out,
    file = map_file,
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE,
    sep = " ",
    append = TRUE
  )

  invisible(NULL)
}

#' Exporta genótipos e mapa de um objeto SNPDataLong no formato FImpute
#'
#' Esta função cria internamente um objeto da classe FImputeExport e exporta os arquivos `.gen` e `.map`.
#'
#' @param object Objeto da classe `SNPDataLong`
#' @param path Diretório onde os arquivos serão salvos (default = "fimpute_run")
#' @param name Nome base dos arquivos gerados (default = "gen_data")
#'
#' @return NULL (arquivos são salvos no disco)
#' @export
setMethod("saveFImpute", "SNPDataLong", function(object, path = "fimpute_run", name = "gen_data") {
  # --- Validações ---
  if (!requireNamespace("snpStats", quietly = TRUE)) {
    stop("O pacote 'snpStats' é necessário. Instale com install.packages('snpStats').")
  }

  if (!inherits(object@geno, "SnpMatrix")) {
    stop("O slot 'geno' do objeto precisa ser da classe 'SnpMatrix'.")
  }

  if (!is.data.frame(object@map)) {
    stop("O slot 'map' do objeto precisa ser um data.frame.")
  }

  fimpute_export <- new("FImputeExport",
                        geno = object@geno,
                        map = object@map,
                        path = path,
                        name = name)

  saveFImpute(fimpute_export)
})
