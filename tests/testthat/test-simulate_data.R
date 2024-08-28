
test_that("Check that simulate_data dimensions match expectations", {
  # Simulate the data
  sim_df <- simulate_data(n = 300, pfh = 15, pnh = 20, mnh = 3,
                          pfi = 5, pni = 1000, sdNoise = 2)

  # What are the number of features in the simulated data?
  number_of_features <- 15 + 20*3 + 5 + 1000

  # Check the items returned
  expect_equal(names(sim_df), c("x", "y", "featureID"))

  # Check the dimensions of the expression dataset
  expect_equal(dim(sim_df$x),
               c(number_of_features, 300))

  # Check the dimensions of the outcome
  expect_equal(length(sim_df$y),
               300)

  # Check the dimensions of the featureID
  expect_equal(length(sim_df$featureID),
               number_of_features)
})
