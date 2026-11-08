library(irace)

# ----------------------------------------------------------------------
# Definition of the functions that return the measure of the box
# according to its size (box_type).
# ----------------------------------------------------------------------
get_area <- function(box_type) {
  if (box_type == 1) return(4)    # T=2
  if (box_type == 2) return(36)   # T=6
  if (box_type == 3) return(100)  # T=10
}

get_volume <- function(box_type) {
  if (box_type == 1) return(8)      # T=2
  if (box_type == 2) return(216)    # T=6
  if (box_type == 3) return(1000)   # T=10
}

get_hypervolume <- function(box_type) {
  if (box_type == 1) return(16)     # T=2
  if (box_type == 2) return(1296)   # T=6
  if (box_type == 3) return(10000)  # T=10
}

# ----------------------------------------------------------------------
# Evaluation function. It runs the executable, selects the measure
# (area, volume or hypervolume) according to the problem dimension and
# returns the cost (fbest).
# ----------------------------------------------------------------------
evaluate_montecarlo <- function(problem_id, configuration) {
  exe <- "/home/user/Descargas/mc_box_esphere"

  densidad <- as.numeric(configuration$densidad)
  box_type <- as.numeric(configuration$box_type)

  dim2_instances <- c(4, 10, 16)
  dim3_instances <- c(1, 3, 5, 12)
  dim4_instances <- c(17, 11, 8, 9)

  pid <- as.numeric(problem_id)
  if (pid %in% dim2_instances) {
    medida <- get_area(box_type)
  } else if (pid %in% dim3_instances) {
    medida <- get_volume(box_type)
  } else if (pid %in% dim4_instances) {
    medida <- get_hypervolume(box_type)
  } else {
    medida <- 1
  }

  # Number of evaluations = point density * measure of the box
  numeval <- as.integer(densidad * medida)

  output <- suppressWarnings(
    system2(exe,
            args = c(problem_id, box_type, numeval),
            stdout = TRUE, stderr = FALSE))

  idx <- grep("^\\s*fbest\\s*=", output)

  cost <- as.numeric(sub("^.*=\\s*", "", output[idx[1]]))
  

  return(list(cost = cost))
}

# ----------------------------------------------------------------------
# Parameter space X
# ----------------------------------------------------------------------
parameters <- readParameters(text = '
densidad "" c (1, 10, 100, 1000)
box_type "" c (1, 2, 3)
')

# ----------------------------------------------------------------------
# Runner
# ----------------------------------------------------------------------
runner <- function(experiment, scenario) {
  C <- experiment$configuration
  evaluate_montecarlo(experiment$instance, C)
}

# ----------------------------------------------------------------------
# Instances (dimension <= 4)
# ----------------------------------------------------------------------
trainInstances <- as.character(c(1, 3, 4, 5, 8, 9, 10, 11, 12, 17, 16))
testInstances  <- trainInstances

# ----------------------------------------------------------------------
# irace scenario
# ----------------------------------------------------------------------
scenario <- list(
  targetRunner        = runner,
  instances           = trainInstances,
  testInstances       = testInstances,
  maxExperiments      = 200,
  logFile             = "irace-results.Rdata",
  debugLevel          = 0,
  deterministic       = TRUE,
  testNbElites        = 1,
  testIterationElites = TRUE,
  elitist             = 1,
  minNbSurvival       = 4
)

set.seed(123456)
irace.output <- irace::irace(scenario = scenario, parameters = parameters)

# ----------------------------------------------------------------------
# Results
# ----------------------------------------------------------------------
elites <- getFinalElites(logFile = "irace-results.Rdata", n = 3)
print(elites)
