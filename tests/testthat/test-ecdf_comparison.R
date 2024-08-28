
test_that("Check ecdf_comparison when using all genes and both tails", {
  # Simulate the input and reference data
  input_simulated_data <- simulate_data(n = 50, pfh = 20, pnh = 20, mnh = 3,
                                        pfi = 20, pni = 2000, sdNoise = 2)
  reference_simulated_data <- simulate_data(n = 200, pfh = 20, pnh = 20, mnh = 3,
                                            pfi = 20, pni = 2000, sdNoise = 1)

  # Compare the CDFs using all genes
  ecdf_comparison_results <- ecdf_comparison(input_simulated_data$x, reference_simulated_data$x)

  # Check that the probabilities are always between 0 and 1
  expect_true(min(ecdf_comparison_results$probabilities) >= 0 &&
                max(ecdf_comparison_results$probabilities) <= 1)

  # Check that the proper tail is returned according the user inputs
  expect_equal(ecdf_comparison_results$tail, "both")

  # Check that the proper tail is returned according the user inputs
  expect_equal(ecdf_comparison_results$gs, 1:nrow(input_simulated_data$x))

})




test_that("Check ecdf_comparison when using only the first 500 genes and the left tail", {
  # Simulate the input and reference data
  input_simulated_data <- simulate_data(n = 75, pfh = 5, pnh = 20, mnh = 5,
                                        pfi = 20, pni = 2000, sdNoise = 2)
  reference_simulated_data <- simulate_data(n = 50, pfh = 5, pnh = 20, mnh = 5,
                                            pfi = 20, pni = 2000, sdNoise = 1)

  # Compare the CDFs using all genes
  use_genes_for_comparison <- rownames(input_simulated_data$x)[1:500]
  ecdf_comparison_results <- ecdf_comparison(input_simulated_data$x, reference_simulated_data$x,
                                             gs = use_genes_for_comparison, tail = "left")

  # Check that the probabilities are always between 0 and 1
  expect_true(min(ecdf_comparison_results$probabilities) >= 0 &&
                max(ecdf_comparison_results$probabilities) <= 1)

  # Check that the proper tail is returned according the user inputs
  expect_equal(ecdf_comparison_results$tail, "left")

  # Check that the proper tail is returned according the user inputs
  expect_equal(ecdf_comparison_results$gs, use_genes_for_comparison)
})
