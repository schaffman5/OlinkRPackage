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


# To use expect_snapshot_file() you'll typically need to start by writing
# a helper function that creates a file from your code, returning a path
save_png <- function(plot, width = 400, height = 400) {
  path <- tempfile(fileext = ".png")
  png(path, width = width, height = height)
  on.exit(dev.off())
  ggplot2::update_geom_defaults("point", list(shape = 17))
  print(
    plot +
      set_plot_theme(font = "")
    )

  path
}

# You'd then also provide a helper that skips tests where you can't
# be sure of producing exactly the same output
expect_snapshot_plot <- function(name, plot) {
  # Other packages might affect results
  skip_if_not_installed("ggplot2", "2.0.0")
  # You'll need to carefully think about and experiment with these skips

  name <- paste0(name, ".png")

  # Announce the file before touching `plot`. This way, if `plot`
  # unexpectedly fails or skips, testthat will not auto-delete the
  # corresponding snapshot file.
  announce_snapshot_file(name = name)

  path <- save_png(plot)
  expect_snapshot_file(path, name)
}




test_that("olink_pca_plot works", {

  # Two Warnings thrown: for dropped assays and droppes samples
  expect_warning(
    expect_warning(
      pca_plot_drop <- npx_data1 %>%
      mutate(SampleID = paste(SampleID, "_", Index, sep = "")) %>%
      olink_pca_plot(drop_assays = TRUE, drop_samples = TRUE)
    )
  )

  expect_snapshot_plot('PCA_plot', pca_plot)
  expect_snapshot_plot('PCA_plot_color_by_treatment', pca_plot_treatCol)
  expect_snapshot_plot('PCA_plot_with_loadings', pca_plot_treatCol_topLoadings)
  expect_snapshot_plot('PCA_plot_drop_assays_and_drop_samples', pca_plot_drop)
})
