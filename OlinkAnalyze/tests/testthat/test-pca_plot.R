set.seed(10)
#Load reference results
refRes_file <- '../data/refResults.RData'
load(refRes_file)

pca_plot <- npx_data1 %>%
  mutate(SampleID = paste(SampleID, "_", Index, sep = "")) %>%
  olink_pca_plot()

pca_plot_treatCol <- npx_data1 %>%
  mutate(SampleID = paste(SampleID, "_", Index, sep = "")) %>%
  filter(!is.na(Treatment)) %>% #Or else, a warning shows up in the test results
  olink_pca_plot(color_g = 'Treatment')

pca_plot_treatCol_topLoadings <- npx_data1 %>%
  mutate(SampleID = paste(SampleID, "_", Index, sep = "")) %>%
  filter(!is.na(Treatment)) %>% #Or else, a warning shows up in the test results
  olink_pca_plot(color_g = 'Treatment',
                 loadings_list = {ref_results$t.test_results %>%
                     head(5) %>%
                     pull(OlinkID)})

test_that("olink_pca_plot works", {

  # Two Warnings thrown: for dropped assays and droppes samples
  expect_warning(
    expect_warning(
      pca_plot_drop <- npx_data1 %>%
      mutate(SampleID = paste(SampleID, "_", Index, sep = "")) %>%
      olink_pca_plot(drop_assays = TRUE, drop_samples = TRUE)
    )
  )

  skip_on_ci()
  vdiffr::expect_doppelganger('PCA plot', pca_plot)
  vdiffr::expect_doppelganger('PCA plot color by treatment', pca_plot_treatCol)
  vdiffr::expect_doppelganger('PCA plot with loadings', pca_plot_treatCol_topLoadings)
  vdiffr::expect_doppelganger('PCA plot drop_assays and drop_samples', pca_plot_drop)
})
