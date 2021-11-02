#'Function which performs a Kruskal-Wallis Test or Friedman Test per protein
#'Performs an Kruskal-Wallis Test for each assay (by OlinkID) in every panel using kruskal.test.
#'Performs an Friedman Test for each assay (by OlinkID) in every panel using friedman_test.
#'The function handles factor variable. \cr\cr
#'Samples that have no variable information or missing factor levels are automatically removed from the analysis (specified in a messsage if verbose = T).
#'Character columns in the input dataframe are automatically converted to factors (specified in a message if verbose = T).
#'Numerical variables are not converted to factors.
#'If a numerical variable is to be used as a factor, this conversion needs to be done on the dataframe before the function call. \cr\cr
#'Inference is specified in a message if verbose = T. \cr
#'The formula notation of the final model is specified in a message if verbose = T. \cr\cr
#'Adjusted p-values are calculated by stats::p.adjust according to the Benjamini & Hochberg (1995) method (“fdr”).
#'The threshold is determined by logic evaluation of Adjusted_pval < 0.05.
#'
#'
#' @param df NPX data frame in long format with at least protein name (Assay), OlinkID, UniProt, Panel and a factor with at least 3 levels.
#' @param variable Single character value.
#' @param outcome Character. The dependent variable. Default: NPX.
#' @param p_adjust_method Adjust P-value for Multiple Comparisons.
#' @param dependence Boolean. Default: FALSE. When the groups are independent, the kruskal-Wallis will run, when the groups are dependent, the Friedman test will run.
#' @param verbose Boolean. Default: True. If information about removed samples, factor conversion and final model formula is to be printed to the console.
#'
#' @return A tibble containing the Kruskal-Wallis Test results for every protein.
#' The tibble is arranged by ascending p-values.
#' @export
#' @examples
#' \donttest{
#' npx_df <- npx_data1 %>% filter(!grepl('control','SampleID', ignore.case = T))
#' #One-way Kruskal-Wallis Test.
#' #Results in a model NPX~Time
#' Kruskal_results <- olink_non_parametric_one(df = npx_df, variable = "Time")}
#' #One-way Friedman Test.
#' #Results in a model NPX~Time
#' Friedman_results <- olink_non_parametric_one(df = npx_df, variable = "Time", dependence = TRUE)}
#' @import dplyr stringr tidyr rstatix broom

olink_non_parametric_one <- function(df,
                        variable,
                        outcome="NPX",
                        p_adjust_method="BH",
                        dependence = FALSE,
                        verbose=T
){

  if(missing(df) | missing(variable)){
    stop('The df and variable arguments need to be specified.')
  }

  withCallingHandlers({

    #Filtering on valid OlinkID
    df <- df %>%
      filter(stringr::str_detect(OlinkID,
                                 "OID[0-9]{5}"))

    #Variables to check
    variable_testers <- intersect(c(variable), names(df))

    ##Remove rows where variables is NA (cant include in analysis anyway)
    removed.sampleids <- NULL

    removed.sampleids <- unique(c(removed.sampleids,df$SampleID[is.na(df[[variable_testers]])]))

    df <- df[!is.na(df[[variable_testers]]),]

    #Not testing assays that have all NA:s
    all_nas <- df  %>%
      group_by(OlinkID) %>%
      summarise(n = n(), n_na = sum(is.na(NPX)),.groups = "drop") %>%
      filter(n == n_na) %>%
      pull(OlinkID)


    if(length(all_nas) > 0) {

      warning(paste0('The assays ',
                     paste(all_nas, collapse = ', '),
                     ' have only NA:s. They will not be tested.'),
              call. = F)

    }

    ##Convert character vars to factor
    converted.vars <- NULL
    num.vars <- NULL
    for(i in variable_testers){
      if(is.character(df[[i]])){
        df[[i]] <- factor(df[[i]])
        converted.vars <- c(converted.vars,i)
      } else if(is.numeric(df[[i]])){
        num.vars <- c(num.vars,i)
      }
    }


    #Not testing assays that have all NA:s in one level
    #Every sample needs to have a unique level of the factor

    nas_in_var <- character(0)


    single_fixed_effects <- variable



    for(effect in single_fixed_effects){

      current_nas <- df %>%
        filter(!(OlinkID %in% all_nas)) %>%
        group_by(OlinkID, !!rlang::ensym(effect)) %>%
        summarise(n = n(), n_na = sum(is.na(NPX)),.groups="drop") %>%
        filter(n == n_na) %>%
        distinct(OlinkID) %>%
        pull(OlinkID)

      if(length(current_nas) > 0) {

        nas_in_var <- c(nas_in_var, current_nas)

        warning(paste0('The assay(s) ',
                       current_nas,
                       ' has only NA:s in atleast one level of ',
                       effect,
                       '. It will not be tested.'),
                call. = F)
      }

      number_of_samples_w_more_than_one_level <- df %>%
        group_by(SampleID, Index) %>%
        summarise(n_levels = n_distinct(!!rlang::ensym(effect), na.rm = T),.groups = "drop") %>%
        filter(n_levels > 1) %>%
        nrow(.)

      if (number_of_samples_w_more_than_one_level > 0) {
        stop(paste0("There are ",
                    number_of_samples_w_more_than_one_level,
                    " samples that do not have a unique level for the effect ",
                    effect,
                    ". Only one level per sample is allowed."))
      }


    }

    formula_string <- paste0(outcome, "~", paste(variable,collapse="*"))


    #Get factors
    fact.vars <- sapply(variable_testers, function(x) is.factor(df[[x]]))
    fact.vars <- names(fact.vars)[fact.vars]


    #Print verbose message
    if(verbose){
      if(!is.null(removed.sampleids) & length(removed.sampleids) >0){
        message("Samples removed due to missing variable ",
                paste(removed.sampleids,collapse=", "))
      }
      if(!is.null(converted.vars)){
        message(paste0("Variables converted from character to factors: ",
                       paste(converted.vars,collapse = ", ")))
      }
      if(!is.null(num.vars)){
        message(paste0("Variables treated as numeric: ",
                       paste(num.vars,collapse = ", ")))
      }
    }


    if (dependence){
      if(verbose){message(paste("Friedman model fit to each assay: "),formula_string)}
      formula_string <- paste0(formula_string,"|","group_id")
      p.val <- df %>%
              filter(!(OlinkID %in% all_nas)) %>%
              filter(!(OlinkID %in% nas_in_var)) %>%
              group_by(Assay,OlinkID, UniProt, Panel,!!!rlang::syms(variable)) %>%
              mutate(group_id = as.factor (1:n()))%>%
              convert_as_factor(!!!rlang::syms(variable)) %>%
              ungroup(!!!rlang::syms(variable))%>%
              group_by(Assay, OlinkID, UniProt, Panel) %>%

              do(friedman_test(as.formula(formula_string),
                                   data=.,na.action=na.omit)) %>%
              ungroup() %>%
              mutate(Adjusted_pval = p.adjust(p, method = p_adjust_method)) %>%
              mutate(Threshold  = ifelse(Adjusted_pval<0.05, "Significant","Non-significant")) %>%
              ungroup() %>%
              arrange(Adjusted_pval)
 }else{
   if(verbose){message(paste("Kruskal model fit to each assay: "),formula_string)}
      p.val <- df %>%
        filter(!(OlinkID %in% all_nas)) %>%
        filter(!(OlinkID %in% nas_in_var)) %>%
        group_by(Assay, OlinkID, UniProt, Panel) %>%
        do(tidy(kruskal.test(as.formula(formula_string),
                              data=.))) %>%
        ungroup() %>%
        mutate(Adjusted_pval=p.adjust(p.value,method=p_adjust_method)) %>%
        mutate(Threshold  = ifelse(Adjusted_pval<0.05,"Significant","Non-significant")) %>%
        ungroup() %>%
        arrange(Adjusted_pval)}

  }, warning = function(w) {
    if (grepl(x = w, pattern = glob2rx("*contains implicit NA, consider using*")))
      invokeRestart("muffleWarning")
  })
}




#'Function which performs an posthoc test per protein.
#'
#'Performs a post hoc test using emmeans::emmeans with Tukey p-value adjustment per assay (by OlinkID) for each panel at confidence level 0.95.
#'See \code{olink_kruskal} for details of input notation. \cr\cr
#'The function handles both factor and numerical variables.
#'The posthoc test for a numerical variable compares the difference in means of the outcome variable (default: NPX) for 1 standard deviation difference in the numerical variable, e.g.
#'mean NPX at mean(numerical variable) versus mean NPX at mean(numerical variable) + 1*SD(numerical variable).
#'
#' @param df NPX data frame in long format with at least protein name (Assay), OlinkID, UniProt, Panel and a factor with at least 3 levels.
#' @param olinkid_list Character vector of OlinkID's on which to perform post hoc analysis. If not specified, all assays in df are used.
#' @param variable Single character value or character array.
#' @param outcome Character. The dependent variable. Default: NPX.
#' @param comparisons A list of length-2 vectors specifying the groups of interest to be compared. For example to compare groups "A" vs "B" and "B" vs "C", the argument is as follow: comparisons = list(c("A", "B"), c("B", "C")).
#' @param p_adjust_method Adjust P-value for Multiple Comparisons.
#' @param verbose Boolean. Deafult: True. If information about removed samples, factor conversion and final model formula is to be printed to the console.
#' @return Tibble of posthoc tests for specicified effect, arranged by ascending adjusted p-values.
#' @export
#' @examples \donttest{
#' test_results <- olink_olink_non_parametric_one(df, "Group")
#' significant_assays <- test_results %>%
#' filter(Threshold == 'Significant') %>%
#' pull(OlinkID)
#'
#' test_posthoc_results <- olink_non_parametric_one_posthoc(df,
#' variable = c("A"),
#' olinkid_list = significant_assays)}

olink_non_parametric_one_posthoc <- function(df,
                                             olinkid_list = NULL,
                                             variable,
                                             outcome="NPX",
                                             comparisons = NULL,
                                             p_adjust_method="BH",
                                             verbose=T
                                             ){



  if(missing(df) | missing(variable)){
    stop('The df and variable and effect arguments need to be specified.')
  }

  withCallingHandlers({

    #Filtering on valid OlinkID
    df <- df %>%
      filter(stringr::str_detect(OlinkID,
                                 "OID[0-9]{5}"))

    if(is.null(olinkid_list)){
      olinkid_list <- df %>%
        select(OlinkID) %>%
        distinct() %>%
        pull()
    }



    variable_testers <- intersect(c(variable), names(df))

    removed.sampleids <- NULL

    removed.sampleids <- unique(c(removed.sampleids,df$SampleID[is.na(df[[variable_testers]])]))

    df <- df[!is.na(df[[variable_testers]]),]


    ##Convert character vars to factor
    converted.vars <- NULL
    num.vars <- NULL

    if(is.character(df[[variable_testers]])){
        df[[variable_testers]] <- factor(df[[variable_testers]])
        converted.vars <- c(converted.vars,variable_testers)
        } else if(is.numeric(df[[variable_testers]])){
        num.vars <- c(num.vars,variable_testers)
        }

    formula_string <- paste0(outcome, "~", paste(variable,collapse="*"))

    #Print verbose message
    if(verbose){
      if(!is.null(removed.sampleids) & length(removed.sampleids) >0){
        message("Samples removed due to missing variable: ",
                paste(removed.sampleids,collapse=", "))
      }
      if(!is.null(converted.vars)){
        message(paste0("Variables converted from character to factors: ",
                       paste(converted.vars,collapse = ", ")))
      }
      if(!is.null(num.vars)){
        message(paste0("Variables treated as numeric: ",
                       paste(num.vars,collapse = ", ")))
      }
      if(any(variable %in% num.vars)){
        message(paste0("Numeric variables post-hoc performed using Mean and Mean + 1SD: ",
                       paste(num.vars[num.vars%in%variable],collapse = ", ")))
      }
      message(paste("Means estimated for each assay from non-parametric model: ",formula_string))
    }



    p.hoc_val <- df %>%
      filter(OlinkID %in% olinkid_list) %>%
      mutate(OlinkID = factor(OlinkID, levels = olinkid_list)) %>%
      group_by(Assay, OlinkID, UniProt, Panel) %>%
      do(wilcox_test(data =., as.formula(formula_string), p.adjust.method = p_adjust_method)) %>%
      ungroup() %>%
      mutate("variable" = variable)

    anova_posthoc_results <- p.hoc_val %>%
      select(all_of(c("Assay", "OlinkID", "UniProt", "Panel","variable","group1","group2","p.adj","p.adj.signif"))) %>%
      mutate(Threshold = if_else(p.adj < 0.05,
                                'Significant',
                                'Non-significant'))



  }, warning = function(w) {
    if (grepl(x = w, pattern = glob2rx("*contains implicit NA, consider using*")))
      invokeRestart("muffleWarning")
  })
}




