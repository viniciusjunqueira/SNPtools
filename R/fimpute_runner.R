#' Executa o FImpute a partir de um objeto FImputeRunner
#'
#' Esta função executa o software externo FImpute a partir de um objeto da classe FImputeRunner,
#' garantindo que os arquivos necessários estejam disponíveis e que os resultados possam ser lidos.
#'
#' @param object Objeto da classe FImputeRunner
#'
#' @return Objeto FImputeRunner com o slot `results` preenchido
#' @export
setGeneric("runFImpute", function(object) standardGeneric("runFImpute"))

#' @export
setMethod("runFImpute", "FImputeRunner", function(object) {
  dir <- object@export@path
  par_file <- object@par_file
  exec <- object@exec_path

  # --- Validações iniciais ---
  if (!dir.exists(dir)) stop("Diretório de trabalho não existe: ", dir)
  if (!file.exists(par_file)) stop("Arquivo fimpute.par não encontrado: ", par_file)
  if (!file.exists(exec)) stop("Executável FImpute não encontrado: ", exec)

  # --- Copiar arquivo fimpute.par ---
  destino_par <- file.path(dir, "fimpute.par")
  file.copy(par_file, destino_par, overwrite = TRUE)

  # --- Remover saída anterior, se houver ---
  output_dir <- file.path(dir, "output_fimpute")
  if (dir.exists(output_dir)) {
    unlink(output_dir, recursive = TRUE)
  }

  # --- Executar o FImpute ---
  command <- paste("cd", shQuote(dir), "&&", shQuote(exec), "fimpute.par")
  message("Executando FImpute com comando: ", command)

  status <- system(command, intern = FALSE, ignore.stdout = TRUE, ignore.stderr = TRUE)
  if (status != 0) stop("Erro ao executar o FImpute. Verifique se o executável está funcionando corretamente.")

  # --- Importar resultados ---
  if (!dir.exists(output_dir)) {
    stop("Pasta de saída 'output_fimpute' não foi criada. A execução do FImpute falhou.")
  }

  if (!exists("read.fimpute", mode = "function")) {
    stop("Função 'read.fimpute()' não está disponível no escopo. Verifique sua definição.")
  }

  res <- read.fimpute(file = output_dir)

  if (!is.data.frame(res)) {
    stop("A função read.fimpute() não retornou um data.frame. Verifique o formato de saída.")
  }

  object@results <- res
  return(object)
})
