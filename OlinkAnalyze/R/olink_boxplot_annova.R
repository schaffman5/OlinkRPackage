#'Function which plots boxplots of a selected variable
#'
#'Generates faceted boxplots of NPX vs. grouping variable for a given list of proteins (OlinkIDs) using ggplot and ggplot2::geom_boxplot.
#'
#' @param df NPX data frame in long format with at least protein name (Assay), OlinkID (unique), UniProt and a grouping variable.
#' @param variable  A character value indiciating which column to use as the x-axis grouping variable
#' @param olinkid_list Character vector indicating which proteins (OlinkIDs) to plot.
#' @param posthoc.results  Results from olink_anova_posthoc function.
#' @param verbose Boolean. If the plots are shown as well as returned in the list (default is false).
#' @param ... coloroption passed to specify color order
#'
#' @return A list of objects of class “ggplot” (the actual ggplot object is entry 1 in the list).
#' @export
#' @examples
#' \donttest{
#' # calculate the p-value for the ANOVA
#' anova_results_oneway <- olink_anova(df = npx_data1, 
#'                                    variable = 'Site')
#' # extracting the significant proteins
#' anova_results_oneway_significant <- anova_results_oneway %>%
#'  filter(Threshold == 'Significant') %>%
#'  pull(OlinkID)
#' anova_posthoc_oneway_results <- olink_anova_posthoc(df = npx_data1,
#'                                                    olinkid_list = anova_results_oneway_significant,
#'                                                    variable = 'Site',
#'                                                    effect = 'Site')
#' olink_boxplot_annova(npx_data1,
#'                     variable="Site",
#'                     olinkid_list=anova_results_oneway_significant,
#'                     posthoc.results= anova_posthoc_oneway_results,
#'                     verbose = F)}




olink_boxplot_annova <- function(df,
                                 variable,
                                 olinkid_list,
                                 posthoc.results,
                                 verbose = F,
                                 ...){

  
  myRound <- function(x){
    if(x >= .00009){
      return(as.character(round(x,4)))
    }else{
      out <- as.character(x)
      if(nchar(out)>8){
        out <- paste0(substring(out,1,4),substring(out,nchar(out)-3,nchar(out)))
      }
      return(return(out))
    }
  }
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
  
  list_of_plots <- list()

  if(is.null(posthoc.results)){
    
    for (i in olinkid_list){
    
      assays_for_plotting <- olinkid_list[i]
      
      npx_for_plotting <- df %>%
        filter(OlinkID %in% assays_for_plotting) %>%
        mutate(OlinkID = factor(OlinkID, levels = assays_for_plotting)) %>%
        select(OlinkID, UniProt, Assay, NPX, eval(variable))
      
      
      boxplot <- npx_for_plotting %>%
        ggplot(aes(y = NPX,
                   x = !!rlang::ensym(variable))) +
        geom_boxplot(aes(fill = !!rlang::ensym(variable))) +
        set_plot_theme() +
        olink_fill_discrete(...)+
        theme(axis.text.x = element_blank(),
              legend.title = element_blank(),
              axis.ticks.x = element_blank(),
              legend.text=element_text(size=13)) 
      
      list_of_plots[[i]] <- boxplot
    }
  }else{
    for (i in olinkid_list){
      assays_for_plotting <- i 
      posthoc.results_temp <- posthoc.results %>% filter(OlinkID==i)
      npx_for_plotting <- df %>%
        filter(OlinkID==i) %>%
        select(OlinkID, UniProt, Assay, NPX, eval(variable)) %>% filter(!is.na(!!rlang::ensym(variable)))
      
      if(!any(posthoc.results_temp$Threshold=="Significant")){
        
        boxplot <- npx_for_plotting %>%
          ggplot(aes(y = NPX,
                     x = !!rlang::ensym(variable))) +
          geom_boxplot(aes(fill = !!rlang::ensym(variable))) +
          set_plot_theme() +
          olink_fill_discrete(...)+
          theme(axis.text.x = element_blank(),
                legend.title = element_blank(),
                axis.ticks.x = element_blank(),
                legend.text=element_text(size=13))
      }else{
        
        star.info <- data.frame(x.vals = levels(npx_for_plotting %>% pull(eval(variable)) %>% as.factor()), 
                                id = 1:length(levels(npx_for_plotting %>% pull(eval(variable)) %>% as.factor())))
        #significant assays
        line.data <-  posthoc.results_temp %>% 
          mutate(C1=sapply(strsplit(contrast," - "),function(x) x[1]),
                 C2=sapply(strsplit(contrast," - "),function(x) x[2])) %>%
          group_by(contrast) %>% 
          mutate(c.sort=min(C1,C2)) %>% 
          mutate(p.value=paste0(myRound(Adjusted_pval)," Contrast: ", contrast)) %>% 
          ungroup() %>% 
          arrange(c.sort) %>% 
          mutate(rowNum=n():1,
                 y.anchor = max(npx_for_plotting$NPX) + rowNum * diff(range(npx_for_plotting$NPX)) *(.5)/max(rowNum)) %>%
          select(contrast,Adjusted_pval,C1,C2,p.value,Threshold,c.sort,y.anchor) %>% 
          pivot_longer(-c(Threshold,contrast,Adjusted_pval,p.value,c.sort,y.anchor),names_to="tmp",values_to = "x.vals") %>% 
          mutate(Star = case_when(Adjusted_pval <0.05 & Adjusted_pval > 0.01 ~ "*",
                                  Adjusted_pval <= 0.01 & Adjusted_pval > 0.005 ~ "**",
                                  Adjusted_pval <= 0.005  ~ "***",
                                  Adjusted_pval >= 0.05 ~ NA_character_)) %>%
          
          #group_by(x.vals) %>% mutate(id=cur_group_id()) %>% ungroup() %>%
          left_join(star.info, by = "x.vals") %>% 
          group_by(contrast) %>% mutate(x.m = sum(id)/2) %>% ungroup() %>% 
          filter(Threshold=="Significant") 
        

            
            
            
        
        boxplot <- npx_for_plotting %>%
          ggplot(aes(y = NPX,
                     x = !!rlang::ensym(variable))) +
          geom_boxplot(aes(fill = !!rlang::ensym(variable))) +
          
          geom_line(data=line.data,aes(x=x.vals,y=y.anchor,group=p.value))+
          geom_text(data=line.data %>% filter(tmp == "C1"),aes(group=p.value,x=x.m,y=y.anchor+0.1,label = Star)) +
          
          set_plot_theme() +
          olink_fill_discrete(...)+
          theme(text=element_text(size=20),legend.position = "none")+
          labs(x=eval(variable),title=paste0(unique( npx_for_plotting$Assay)," - ",assays_for_plotting)) +
          theme(plot.title = element_text(hjust = 0.5),axis.title=element_text(size=14,face="bold"),legend.title=element_text(size=12), 
                legend.text=element_text(size=10))
        
      }
      list_of_plots[[i]] <- boxplot
    }     
  }
  
  
  if(verbose){
    show(boxplot)
  }
  
  return(invisible(list_of_plots))
  
}
