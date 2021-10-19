#'Function which plots boxplots of a selected variable
#'
#'Generates faceted boxplots of NPX vs. grouping variable for a given list of proteins (OlinkIDs) using ggplot and ggplot2::geom_boxplot.
#'
#' @param df NPX data frame in long format with at least protein name (Assay), OlinkID (unique), UniProt and a grouping variable.
#' @param TtestResult The output dataframe from Olink_ttest function.
#' @param variable  A character value indiciating which column to use as the x-axis grouping variable
#' @param olinkid_list Character vector indicating which proteins (OlinkIDs) to plot.
#' @param pair_id Character value indicating which column indicates the paired sample identifier.
#' @param verbose Boolean. If the plots are shown as well as returned in the list (default is false).
#' @param ... coloroption passed to specify color order
#'
#' @return A list of objects of class “ggplot” (the actual ggplot object is entry 1 in the list).
#' @export
#' @examples
#' \donttest{
#' npx_df <- npx_data1 %>% filter(!grepl('control',SampleID, ignore.case = TRUE))
#' ttest_results <- olink_ttest(df=npx_df,
#'                              variable = 'Treatment',
#'                              alternative = 'two.sided')
#' significant_assays <- ttest_results %>%
#'     filter(Threshold == 'Significant') %>%
#'     pull(OlinkID)
#' olink_boxplot_ttest(npx_df, 
#' TtestResult = ttest_results, 
#' variable = 'Treatment',
#' olinkid_list = significant_assays)}

olink_boxplot_ttest <- function(df,
                                TtestResult, 
                                variable,
                                olinkid_list,
                                verbose = F,
                                pair_id,
                                ...){
  #checking ellipsis
  if(length(list(...)) > 0){
    
    ellipsis_variables <- names(list(...))
    
    if(length(ellipsis_variables) == 1){
      
      if(!(ellipsis_variables == 'coloroption')){
        
        stop(paste0('The ... option only takes the coloroption argument. ... currently contains the variable ',
                    ellipsis_variables,
                    '.'))
        
      }
      
    }else{
      
      stop(paste0('The ... option only takes one argument. ... currently contains the variables ',
                  paste(ellipsis_variables, collapse = ', '),
                  '.'))
    }
  }
  
  
  #Filtering on valid OlinkID
  df <- df %>%
    filter(stringr::str_detect(OlinkID,
                               "OID[0-9]{5}"))
  
  #Column setup
  columns_for_npx_data <- c("OlinkID","UniProt","Assay", "NPX", eval(variable))
  
  #Testing that needed columns are correct
  if(!(all(columns_for_npx_data %in% colnames(df)))){
    
    
    stop(paste0("Column(s) ",
                paste(setdiff(columns_for_npx_data,
                              colnames(df)),
                      collapse=", "),
                " not found in NPX data frame!"))
    
  }
  
  #Ttest Results
  p.val <- TtestResult %>% mutate(Star = case_when(Adjusted_pval <0.05 & Adjusted_pval > 0.01 ~ "*",
                                                   Adjusted_pval <= 0.01 & Adjusted_pval > 0.005 ~ "**",
                                                   Adjusted_pval <= 0.005  ~ "***",
                                                   Adjusted_pval >= 0.05 ~ NA_character_))
  #Setup
  topX <- length(olinkid_list)
  
  list_of_plots <- list()
  COUNTER <- 1
  
  for (i in c(1:topX)){
    
    assays_for_plotting <- olinkid_list[i]
    
    
    npx_for_plotting <- df %>%
      filter(OlinkID %in% assays_for_plotting) %>%
      mutate(OlinkID = factor(OlinkID, levels = assays_for_plotting)) %>%
      select(OlinkID, UniProt, Assay, NPX, eval(variable))
    
    #Star for groups
    star_group <-  npx_for_plotting %>% pull(eval(variable)) %>% unique()
    
    #Star to add for the assays
    anno <- as.character(p.val[p.val$OlinkID == olinkid_list[i],"Star"])
    
    boxplot <- npx_for_plotting %>%
      ggplot(aes(y = NPX,
                 x = !!rlang::ensym(variable))) +
      geom_boxplot(aes(fill = !!rlang::ensym(variable))) +
      
      ggsignif::geom_signif(comparisons=list(c(star_group[1], star_group[2])), annotation=anno,
                            y_position = max(npx_for_plotting$NPX)+diff(range(npx_for_plotting$NPX))*c(.05), tip_length = 0, vjust=0.4) +
      
      
      set_plot_theme() +
      olink_fill_discrete(...)+
      theme(axis.text.x = element_blank(),
            legend.title = element_blank(),
            axis.ticks.x = element_blank(),
            legend.text=element_text(size=13))
    
    list_of_plots[[COUNTER]] <- boxplot
    COUNTER <- COUNTER + 1
    
  }
  
  if(verbose){
    show(boxplot)
  }
  
  return(invisible(list_of_plots))
  
}
