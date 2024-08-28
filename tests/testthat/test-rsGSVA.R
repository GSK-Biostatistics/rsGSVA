

test_that("Check that rsGSVA and classic GSVA equal each other when the input dataset is the reference", {
  # Simulate the data
  sim_df <- simulate_data()

  # What genes should be used for the gene sets?
  gene_sets <- list(gs1 = rownames(sim_df$x)[sim_df$featureID=="fh"],
                    gs2 = rownames(sim_df$x)[sim_df$featureID=="nh"])

  # Check the classic GSVA
  classic_gsva <- GSVA::gsva(sim_df$x,
                             gene_sets)

  # Check rsGSVA with the same dataset as the reference
  rs_gsva <- rsGSVA(sim_df$x, gene_sets, train_expr = sim_df$x)

  # Check that these equal
  expect_equal(classic_gsva, rs_gsva)
})




test_that("Check that the GSVA scores are different when an input dataset is used and the dimensions are appropriate", {
  # Simulate the data
  sim_df <- simulate_data()
  train_df <- simulate_data(n = 300)

  # What genes should be used for the gene sets?
  gene_sets <- list(gs1 = rownames(sim_df$x)[sim_df$featureID=="fh"],
                    gs2 = rownames(sim_df$x)[sim_df$featureID=="nh"])

  # Check the classic GSVA
  classic_gsva <- GSVA::gsva(sim_df$x,
                             gene_sets)

  # Check rsGSVA with the same dataset as the reference
  rs_gsva <- rsGSVA(sim_df$x, gene_sets, train_expr = train_df$x)

  # Check that these equal
  expect_false(all.equal(classic_gsva, rs_gsva)==TRUE)

  # Check the dimensions
  expect_equal(dim(rs_gsva),
               c(length(gene_sets), ncol(sim_df$x)))
})
