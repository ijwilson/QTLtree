test_examples <- function(path = "../../man") {
  man <- dir(path, "\\.Rd$", full.names = TRUE)
  lapply(man, test_example)
}
