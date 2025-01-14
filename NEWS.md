# OlinkAnalyze 1.3.0

## Major Changes

* DESCRIPTION file updated to include all authors and maintainers (@kathy-nevola, #66)
* Unit testing was added (@marisand, @simfor, #65, #55, #47)
* Continuous Integration (CI) was added (@AskPascal, #57, #53, #40, #31)
* `read_NPX` now converts Panel Names to Sentence case except Roman Numerals (@kathy-nevola, #44)
* Output order of `olink_bridgeselector()` is now randomized (@marisand, #60)
* PlateSize argument was added to `displayPlateLayout()` and `olink_plate_randomizer` for T48 and T96. Hardcoded sizes have been removed. (@marisand, #47)
* Changed `olink_pal()` to have gray instread of light blue (@marisand, #22)

## Bug Fixes

* `set_plot_theme()` will now load Swedish Gothic Thin if available on all OS (@marisand, @AskPascal, #70, #39)
* Help documentation was updated to correct typos and clarify notation for ANOVA and LME models (@kathy-nevola, @marisand,@AskPascal #67, #33, #37)
* Fix `olink_dist_plot()` from showing multiple bars when sample has QC warning in some assays (@marisand, @AskPascal, #64)
* Fix `olink_qc_plot()` from showing multiple points/labels when sample has QC warning in some assays (@marisand, @AskPascal, #63)
* Panel names in npx_data1 and npx_data2 were fixed to match new `read_NPX()` formatting (@kathy-nevola, #62)
* `olink_bridgeselector()` QC warning filtering changed to match `olink_qc_plot()` and `olink_dist_plot()` (@marisand, @AskPascal, #61)
* global variables are now specified (@kathy-nevola, #45)
* Function documentation updated to remove notes on R CMD check (@kathy-nevola, #39)
* Fix summary casting "tukey" to "sidak" adjustment warning in `olink_anova_posthoc()` and `olink_lmer_posthoc()` functions (@marisand, #38)
* Update functions to import selectively (@kathy-nevola, @Orbmac, @AskPascal, @marisand, #21, #20, #19, #18, #15, #29)
* Fix guides size argument in `olink_pca_plot()` (#17) (@AskPascal, #24)






