# SWG DATA CLEANING & ANALYSIS PRELIMINARIES
# LR 2025

###############################################################################
#-------------------PACKAGES, FUNCTIONS, & GLOBAL VARIABLES--------------- ----
###############################################################################
# a. Packages                                                              ----

# 1. data cleaning / manipulation
library(tidyverse) # includes dplyr, ggplot2, lubridate, purrr, readxl
library(corrplot) # run correlation matrices & visualize collinearity in predictor variables
# library(janitor)
library(timeDate)
# library(pollen)
library(data.table) # PRISM solution 2
library(slider) # PRISM solution 3

# 2. GLMM & model selection
library(glmmTMB)
library(car)
#library(multcomp)
library(performance) # running VIFs (variance inflation factor tests)
library(MASS) # model selection - check model fit & generate best matched model
library(MuMIn) # model selection

# 3. community composition 
library(vegan)

# 4. network modeling
library(bipartite)

# b. Global Variables                                                      ----
# create strings with different bee category values
Method = as.factor(c("netting", "Pan Trap"))
Order = as.factor(c("Hymenoptera"))
Family = as.factor(c("Mellitidae", "Megachilidae", "Halictidae", "Andrenidae", "Apidae", "Colletidae"))
year = as.factor(c("2017", "2018","2019","2021","2022"))

# d. Custom functions - LR                                                 ----
site.plot <- function(x) {
  
  x <- x |>
    mutate(site = recode_factor(site, "Cross Timbers" = "CrossTimbers")) |>
    filter(site != c("Oka")) |> # Filter out Oka and Union Grove plots
    filter(plot != c("Union Grove"))  |>
    mutate_if(is.factor, funs(factor(trimws(.)))) # trim white space from factors (b/c of "Wagley" and "Wagley ") 
  x$site <- as.factor(x$site)
  x$plot <- as.factor(x$plot)
  x$site.plot <- as.factor(paste5(x$site, x$plot, sep=" ", na.rm=T))
  droplevels(x)
  return(x)
}
reclass_year <- function(x) {
  x$year <- year(as.Date(as.character(x$year), format = '%Y'))
  return(x)
}

# e. Custom functions - EL (do not edit)                                   ----
# Credit: Elinor Lichtenberg
                                 
# [paste5] combine multiple character columns
paste5 <- function(..., sep = " ", collapse = NULL, na.rm = F) 
{
  if (na.rm == F)
    paste(..., sep = sep, collapse = collapse)
  else
    if (na.rm == T) {
      paste.na <- function(x, sep) {
        x <- gsub("^\\s+|\\s+$", "", x)
        ret <- paste(na.omit(x), collapse = sep)
        is.na(ret) <- ret == ""
        return(ret)
      }
      df <- data.frame(..., stringsAsFactors = F)
      ret <- apply(df, 1, FUN = function(x) paste.na(x, sep))
      
      if (is.null(collapse))
        ret
      else {
        paste.na(ret, sep = collapse)
      }
    }
}
    # calculate diversity ----
Evar <- function(x) # subfunction of {diversity calculate} 
{
  x1 <- data.frame(abundance=x[x>0]) #Remove taxa with zero abundance
  
  S <- nrow(x1)
  x1$v1 <- log(x1$abundance)
  x1$v2 <- x1$v1/S
  s1 <- sum(x1$v2)
  x1$v3 <- ((x1$v1 - s1)^2)/S
  s2 <- sum(x1$v3)
  evar <- 1 - 2/pi * atan(s2)
  
  return(evar)
}

# [diversity_calculate] Calc abundance, richness, evenness 
diversity_calculate <- function(x, grps, evenness=T, Chao=F, rarefaction=F, 
                                permutations=1000)
{
  require(tidyverse)
  
  #Remove taxon from list of grouping variables
  #grps2 <- syms(grps)[1:(length(grps)-1)]
  
  #Calculate abundance for each group
  x1 <- x %>%  
    #group_by(!!!syms(grps)) %>% 
    group_by("grps") |> # doesn't work yet...
    tally() %>% 
    summarise(abun=sum(n))
  
  #Calculate richness, and evenness when specified, for each group
  if (evenness)
  {
    x2 <- x %>% 
      #filter(!diversity_exclude) %>% 
      droplevels() %>% 
      #group_by(!!!syms(grps)) %>% 
      group_by("grps") %>% 
      tally() %>% 
      summarise(rich=n_distinct("Genus_species"),
                Evar3=Evar(n)) %>% 
      mutate(Evar3=ifelse(rich>=3, Evar3, NA))
  } else
  {
    x2 <- x %>% 
      #filter(!diversity_exclude) %>% 
      droplevels() %>% 
      #group_by(!!!syms(grps)) %>% 
      group_by("grps") %>% 
      tally() %>% 
      summarise(rich=n_distinct("Genus_species"))
  }
  
  #Merge abundance & diversity
  x3 <- x1 %>% 
    #left_join(x2) %>% 
    cbind(x2) %>% 
    ungroup()
  
  #Calculate estimated richness (Chao1) and merge with raw values, when specified
  if (Chao)
  {
    #Select data, calculate Chao
    xchao <- x %>% 
      #filter(!diversity_exclude) %>% 
      droplevels() %>% 
      #arrange(!!!syms(grps)) %>% 
      arrange("grps") %>% 
      #group_by(!!!syms(grps)) %>% 
      group_by("grps") %>% 
      tally() %>%
      do(getChao(.))
    
    #Merge with abundance & diversity
    x3 <- x3 %>% 
      left_join(xchao)
  }
  
  if (rarefaction)
  {
    
    #Re-calculate abundances without the diversity exclude individuals
    x1B <- x %>%  
      #filter(!diversity_exclude) %>% 
      group_by(!!!syms(grps)) %>% 
      tally() %>% 
      summarise(abun=sum(n))
    
    #Determine minimum number of individuals (>=20) at any one plot, add to the list of all individuals, conduct rarefaction on each specified group; execute this in different ways if there is no grouping variable (other than site) versus if diversity is being calculated for subsets of the data
    if (length(grps)==2){
      raremin <- x1B %>% 
        filter(abun>=20) %>% 
        summarise(minindivs=min(abun, na.rm=T)) %>% 
        pull(minindivs)
      xrare <- x %>% 
        filter(!diversity_exclude) %>% 
        mutate(minindivs=raremin) %>%
        group_by(!!!syms(grps)) %>%
        summarise(n=n(), minindivs=first(minindivs)) %>%
        do(indiv.rarefy(.$taxon, .$n, .$minindivs, permutations))
    } else if (length(grps)>=2)
    {
      xrare <- x1B %>% 
        filter(abun>=20) %>%
        group_by(!!!syms(grps)[2:(length(grps)-1)]) %>%
        summarise(minindivs=min(abun, na.rm=T)) %>%
        right_join(x) %>%
        filter(!diversity_exclude) %>% 
        group_by(!!!syms(grps)) %>%
        summarise(n=n(), minindivs=first(minindivs)) %>%
        do(indiv.rarefy(.$taxon, .$n, .$minindivs, permutations))
    }
    
    #Remove rarefied evenness when the user does not want evenness
    if (!evenness) xrare <- xrare %>% select(-rareEvar3)
    
    #Merge with raw abundance & diversity
    x3 <- x3 %>% 
      left_join(xrare)
  }
  
  return(x3)
}

# [indiv.rarefy] - subfunction of {diversity calculate} 
indiv.rarefy <- function(taxa, ntaxon, samplemin, permutations)
{
  require(tidyverse)
  
  #Check that have more than one species and at least 20 individuals
  if((length(taxa)<=1) | sum(ntaxon, na.rm=T)<20)
  {
    return(tibble(rarerich=NA, rareEvar3=NA))
  }
  
  #Expand to one row per individual
  eachflower <- rep(taxa, times=ntaxon)
  
  #Set up tibble as long as the number of permutations desired; generate random subsamples for each iteration (as a list tibble variable); calculate richness and evenness of each subsample; set evenness to NA if fewer than 3 species in the subset; calculate average rarefied richness and evenness across subsamples
  perm <- tibble(whichperm=seq(1:permutations)) %>% 
    group_by(whichperm) %>% 
    do(., inflorescences=sample(eachflower, samplemin, replace=F)) %>% 
    unnest() %>% 
    group_by(whichperm, inflorescences) %>% 
    tally() %>% 
    summarise(rarerich.perm=n(),
              rareEvar3.perm=ifelse(rarerich.perm>=3, Evar(n), NA)) %>% 
    summarise(rarerich=mean(rarerich.perm, na.rm=T),
              rareEvar3=mean(rareEvar3.perm, na.rm=T))
  
  return(perm)
}

getChao <- function(xC) # subfunction of {diversity calculate} 
{
  require(tidyverse)
  require(vegan)
  
  #Store grouping variables
  #grps2 <- grps[1:(length(grps)-1)]
  grps2 <- names(xC %>% dplyr::select(group_cols(.)))
  
  #Make site by species tibble, remove grouping variables
  xC2 <- xC %>%
    select(taxon, n) %>%
    spread(taxon, n, fill=0) %>%
    ungroup() %>%
    select(-grps2) %>%
    do(as.data.frame(t(estimateR(.)))) %>%
    select(2:3) %>%
    rename(rich_Chao1=S.chao1,
           rich_Chao1.se=se.chao1) %>%
    replace(is.nan(.), NA)
  
  #Return tibble with Chao1 values
  return(xC2)
}

is.nan.data.frame <- function(x) # detect NaNs in df
{
  do.call(cbind, lapply(x, is.nan))
}

# [cor.test.nofail] store correlation test output
cor.test.nofail <- function(x, y, alternative="two.sided", method, exact=NULL, conf.level=0.95, continuity=F)
{
  curtest <- try(cor.test(x, y, alternative=alternative, method=method, exact=exact, conf.level=conf.level, continuity=continuity), silent=T)
  if (inherits(curtest, "try-error"))
  {
    curtest <- list(statistic=NA, parameter=NA, p.value=NA, estimate=NA, null.value=0, alternative=alternative, method=method, data.name=NA, conf.int=NA)
    class(curtest) <- "htest"
  }
  
  return(curtest)
}

# [summary_to_file] visual/numerical Summary
summary_to_file <- function(df, wd, filename, cont_vars=NULL, cat_vars=NULL, 
                            no_graph=NULL, no_table=NULL, DEBUG=F)
{
  #Load necessary libraries
  require(tidyverse)
  
  #Store original working directory
  wd0 <- getwd()
  #Set the working directory to where the output should go
  setwd(wd)
  #Store the current system time to use in the file name
  curtime <- format(Sys.time(), "%b%d%H%M%S")
  
  #Store current par() values
  parcurrent <- par(no.readonly=T)
  
  #Store continuous variables
  if (is.null(cont_vars))
  {
    df_cont <- df %>% 
      select_if(is.numeric) %>% 
      data.frame()
  } else
  {
    df_cont <- df %>% 
      dplyr::select(all_of(cont_vars)) %>%
      data.frame()
    #Check that all columns are indeed numeric. If not, re-set the working directory and stop the function.
    cont_check <- all(sapply(df_cont, is.numeric))
    if (!cont_check)
    {
      setwd(wd0)
      stop("Error: Some of the specified continuous variables are not numeric.")
    }
  }
  
  #Store continuous variables for which histograms should be plotted & convert tibble to data frame
  if (ncol(df_cont)>0 & !is.null(no_graph))
  {
    df_cont_hist <- df_cont %>% 
      dplyr::select(-all_of(no_graph)) %>%
      data.frame()
  } else
  {
    df_cont_hist <- df_cont %>% data.frame()
  }
  
  #Store categorical variables
  if (is.null(cat_vars))
  {
    df_cat <- df %>% 
      select_if(negate(is.numeric)) %>% 
      mutate_if(is.character, as.factor) %>% 
      droplevels() %>% 
      data.frame()
  } else
  {
    df_cat <- df %>% 
      dplyr::select(all_of(cat_vars)) %>% 
      data.frame()
    #Check that all columns are indeed non-numeric. If not, re-set the working directory and stop the function.
    cat_check <- all(sapply(df_cat, negate(is.numeric)))
    if (!cat_check)
    {
      setwd(wd0)
      stop("Error: Some of the specified categorical variables are not categorical.")
    }
  }
  
  #Generate flags for whether have continuous and categorical variables, variables the user wants histograms of
  if (ncol(df_cont)>0) cont_flag <- T else cont_flag <- F
  if (ncol(df_cont_hist)>0) hist_flag <- T else hist_flag <- F
  if (ncol(df_cat)>0) cat_flag <- T else cat_flag <- F
  if (DEBUG)
  {
    print(paste0("cont_flag: ", cont_flag, "; cat_flag: ", cat_flag))
  }
  
  if (DEBUG)
  {
    print("Continuous variables:", quote=F)
    print(names(df_cont))
    print("Categorical variables: ", quote=F)
    print(names(df_cat))
  }
  
  #Start sending text output to a file
  sink(file=paste(filename, curtime, "summaries.txt", sep="_"), type="output")
  
  #Summarize the data
  cat("Overall summary:", "\n", "\n")
  print(summary(df))
  cat("\n", rep("*", 30), "\n", "\n")
  if (cat_flag)
  {
    cat("Categorical variable level frequencies:", "\n", "\n")
    #Tables of categorical variables (removing any variables the user specifies not to include)
    cat_temp <- df_cat %>% 
      select_if(is.factor) %>% 
      data.frame()
    if (!is.null(no_table))
    {
      cat_temp <- df_cat %>% 
        dplyr::select(-all_of(no_table)) %>% 
        data.frame()
    }
    options(tibble.print_max = Inf, max.print=9999)
    cat_temp %>% 
      droplevels() %>% 
      gather("var", "value", factor_key=T) %>% 
      count(var, value) %>% 
      print()
    #print(n = Inf)
    options(tibble.print_max = 20, max.print=1000)
  } else cat("No categorical variables to summarize")
  
  #Stop sending text output to a file
  sink()
  
  #Graph continuous variables
  if (hist_flag)
  {
    #Histogram function that is called by lapply
    myhistogram <- function(x)
    {
      hist(df_cont_hist[, x], main="", xlab=x)
    }
    
    #Start sending graphic output to a pdf file
    pdf(file=paste(filename, curtime, "graphs.pdf", sep="_"), width=6, height=9, onefile=T, title="Descriptive plots", paper="letter")
    
    #Two plots per page
    par(mfrow=c(2, 1))
    
    #Histograms of continuous variables
    if (DEBUG) lapply(names(df_cont_hist), function(x) {print(class(df_cont_hist[, x]))})
    lapply(names(df_cont_hist), myhistogram)
    
    #Stop sending plots to file
    dev.off()
  }
  
  #Re-set the working directory to the default
  setwd(wd0)
  
  #Re-set par to what it was when the function was called
  par(parcurrent)
}  

# [vars_cors2] visual/statistical relationship
var_cors2 <- function(df, cortype, tabletype, sampletype, wd, filename, 
                      cont_vars=NULL, cat_vars=NULL, no_cont=NULL, no_cat=NULL, DEBUG=F)
{
  #Load necessary libraries
  require(tidyverse)
  require(PerformanceAnalytics)
  require(Hmisc)
  require(corrr)
  require(plotrix)
  #require(PMCMRplus)
  
  #Store original working directory
  wd0 <- getwd()
  #Set the working directory to where the output should go
  setwd(wd)
  #Store the current system time to use in the file name
  curtime <- format(Sys.time(), "%b%d%H%M%S")
  
  #Store current par() values
  #parcurrent <- par(no.readonly=T)
  
  #Store continuous variables
  if (is.null(cont_vars))
  {
    df_cont <- df %>% 
      select_if(is.numeric) %>% 
      data.frame()
  } else
  {
    df_cont <- df %>% 
      dplyr::select(all_of(cont_vars)) %>% 
      data.frame()
    #Check that all columns are indeed numeric. If not, re-set the working directory and stop the function.
    cont_check <- all(sapply(df_cont, is.numeric))
    if (!cont_check)
    {
      setwd(wd0)
      stop("Error: Some of the specified continuous variables are not numeric.")
    }
  }
  #Remove any variables the user specifies to exclude, first checking that they are in the current tibble
  cont_overlap <- intersect(names(df_cont), no_cont)
  if (!is.null(cont_overlap))
  {
    df_cont <- df_cont %>% 
      dplyr::select(-all_of(cont_overlap)) %>% 
      data.frame()
  }
  #Store number of continuous variables
  n_cont <- ncol(df_cont)
  
  #Function that returns TRUE if the variable is a factor or is logical
  is.factor_logical <- function(x) {is.factor(x) | is.logical(x)}
  
  #Store categorical variables & number of categorical variables
  if (is.null(cat_vars))
  {
    df_cat <- df %>% 
      select_if(is.factor_logical) %>% 
      droplevels() %>% 
      data.frame()
  } else
  {
    df_cat <- df %>% 
      dplyr::select(all_of(cat_vars)) %>% 
      data.frame()
    #Check that all columns are indeed non-numeric. If not, re-set the working directory and stop the function.
    cat_check <- all(sapply(df_cat, is.factor_logical))
    if (!cat_check)
    {
      setwd(wd0)
      stop("Error: Some of the specified categorical variables are not categorical.")
    }
  }
  #Remove any variables the user specifies to exclude, first checking that they are in the current tibble
  cat_overlap <- intersect(names(df_cat), no_cat)
  if (!is.null(cat_overlap))
  {
    df_cat <- df_cat %>% 
      dplyr::select(-all_of(cat_overlap)) %>% 
      data.frame()
  }
  #Store number of categorical variables
  n_cat <- ncol(df_cat)
  
  if (DEBUG)
  {
    print("Continuous variables:", quote=F)
    print(names(df_cont))
    print("Categorical variables: ", quote=F)
    print(names(df_cat))
  }
  
  #Generate flags for whether have continuous (any, and at least 2) and categorical variables
  if (n_cont>0) cont_flag <- T else cont_flag <- F
  if (n_cont>1) cont_pair_flag <- T else cont_pair_flag <- F
  if (n_cat>0) cat_flag <- T else cat_flag <- F
  if (n_cat>1) cat_pair_flag <- T else cat_pair_flag <- F
  if (DEBUG)
  {
    print(paste0("cont_flag: ", cont_flag, "; cont_pair_flag: ", cont_pair_flag, "; cat_flag: ", cat_flag, "; cat_pair_flag: ", cat_pair_flag))
  }
  
  #Graphing functions
  #Scatter plot
  myscatter <- function(x, y)
  {
    plot(df_cont[, x], df_cont[, y], main="", xlab=x, ylab=y)
  }
  #Boxplot
  myboxplot <- function(x, y)
  {
    boxplot(df_cont[, y]~addNA(df_cat[, x]), main="", xlab=x, ylab=y)
  }
  
  #Start sending graphic outputs, if any, to a PDF file
  if(cont_pair_flag | (cont_flag & cat_flag))
  {
    pdf(file=paste(filename, curtime, "correlations.pdf", sep="_"), onefile=T, title="Correlation and box plots", paper="letter")
    
    if (cont_pair_flag)
    {
      #Convert continuous variables tibble to matrix
      m_cont <- as.matrix(df_cont)
      #Correlations among potential response variables
      chart.Correlation(m_cont, histogram=T, method=cortype)
      
      #Make list of continuous-continuous variable pairs
      point_pairs <- df_cont %>% 
        rownames_to_column(var="IDs") %>%
        gather(contvar, contvalue, names(df_cont), -IDs) %>%
        inner_join(., ., by="IDs") %>%
        mutate(temp=ifelse(as.character(contvar.x)<as.character(contvar.y), paste(contvar.x, contvar.y), paste(contvar.y, contvar.x))) %>%
        distinct(temp, .keep_all=T) %>% 
        filter(contvar.x!=contvar.y) %>% 
        mutate(Var1=ifelse(as.character(contvar.x)<as.character(contvar.y), contvar.x, contvar.y),
               Var2=ifelse(as.character(contvar.x)<as.character(contvar.y), contvar.y, contvar.x)) %>%
        arrange(Var1, Var2) %>%
        dplyr::select(all_of(Var1, Var2))
      #Scatter plots of each pair of continuous variables
      mapply(myscatter, point_pairs$Var1, point_pairs$Var2)
    }
    
    if (cont_flag & cat_flag)
    {
      #Make list of categorical-continuous variable pairs
      box_pairs <- crossing(names(df_cat), names(df_cont)) %>% 
        rename(Var1="names(df_cat)", Var2="names(df_cont)") %>%
        arrange(Var2, Var1)
      #Box plots of each continuous variable against each categorical variable
      mapply(myboxplot, box_pairs$Var1, box_pairs$Var2)
    }
    
    #Stop sending graphic output to a PDF file
    dev.off()
  }
  
  #This function stores the output of correlation tests, storing NA values when the test fails (rather than breaking)
  cor.test.nofail <- function(x, y, alternative="two.sided", method, exact=NULL, conf.level=0.95, continuity=F)
  {
    curtest <- try(cor.test(x, y, alternative=alternative, method=method, exact=exact, conf.level=conf.level, continuity=continuity), silent=T)
    if (inherits(curtest, "try-error"))
    {
      curtest <- list(statistic=NA, parameter=NA, p.value=NA, estimate=NA, null.value=0, alternative=alternative, method=method, data.name=NA, conf.int=NA)
      class(curtest) <- "htest"
    }
    
    return(curtest)
  }
  
  #Send correlations to CSV file
  #Calculate pairwise correlations, removing duplicated variable pairs and sorting to make the file easier to read
  if (cont_pair_flag)
  {
    corrtable <- df_cont %>%
      rownames_to_column(var="IDs") %>%
      gather(contvar, contvalue, names(df_cont), -IDs) %>%
      inner_join(., ., by="IDs") %>%
      mutate(temp=ifelse(as.character(contvar.x)<as.character(contvar.y), paste(contvar.x, contvar.y), paste(contvar.y, contvar.x))) %>% 
      distinct(IDs, temp, .keep_all=T) %>% 
      filter(contvar.x!=contvar.y) %>% 
      group_by(contvar.x, contvar.y, .drop=T) %>% 
      summarise(teststat=cor.test.nofail(contvalue.x, contvalue.y, method=cortype)$statistic,
                r.estimate=cor.test.nofail(contvalue.x, contvalue.y, method=cortype)$estimate,
                pvalue=cor.test.nofail(contvalue.x, contvalue.y, method=cortype)$p.value) %>%
      mutate(sig05=ifelse(pvalue<0.05, T, F),
             variable1=ifelse(as.character(contvar.x)<as.character(contvar.y), contvar.x, contvar.y),
             variable2=ifelse(as.character(contvar.x)<as.character(contvar.y), contvar.y, contvar.x)) %>% 
      ungroup() %>% 
      arrange(variable1, variable2) %>% 
      dplyr::select(all_of(variable1, variable2, teststat:sig05))
    #Output correlations to CSV file
    write.csv(corrtable, file=paste(filename, curtime, "correlations.csv", sep="_"), row.names=F)
  }
  
  #Start sending text output to file
  sink(file=paste(filename, curtime, "correlations.txt", sep="_"), type="output")
  
  #Correlations among potential response variables
  if (cont_pair_flag)
  {
    cat("Correlations among continuous variables: see CSV file", "\n")
  } else cat("No pairs of continuous variables")
  cat("\n", rep("*", 30), "\n", "\n")
  
  if (cat_pair_flag)  #Change from cat_flag 10/24/18
  {
    cat("Relationships among categorical variables", "\n", "\n")
    for (i in 1:(n_cat-1))
      for (j in (i+1):n_cat)
      {
        #Print variable names
        cat("** ", names(df_cat)[i], " and ", names(df_cat)[j], "\n", "\n")
        #Contingency table
        print(table(df_cat[, i], df_cat[, j], useNA="ifany"))
        
        #Test for independence
        #Avoid function stopping when a variable has only one level
        if (tabletype=="chisq")
        {
          curtest <- try(chisq.test(df_cat[, i], df_cat[, j]), silent=T)
          curresid <- curtest
        } else if (tabletype=="fisher")
        {
          if (nrow(df_cat)==2 & ncol(df_cat==2))
          {
            curtest <- try(fisher.test(df_cat[, i], df_cat[, j]), silent=T)
          } else #Simulate p-values when table larger than 2x2, to avoid errors
          {
            curtest <- try(fisher.test(df_cat[, i], df_cat[, j], simulate.p.value=T), silent=T)
          }
          curresid <- try(chisq.test(df_cat[, i], df_cat[, j]), silent=T)
        } else
        {
          cat("Unsupported test of independence", "\n")
        }
        if (inherits(curtest, "try-error"))
        {
          cat("No test of independence run", "\n")
        } else
        {
          print(curtest)
          cat("\n", "Standardized residuals (from chi-squared test)", "\n")
          print(curresid$stdres)
        }
        
        cat("\n")
      }
  } else cat("No pairs of categorical variables")
  cat("\n", rep("*", 30), "\n", "\n")
  
  if (cont_flag & cat_flag)
  {
    cat("Relationships between categorical and continuous variables", "\n", "\n")
    for (i in 1:n_cat)
      for (j in 1:n_cont)
      {
        #Print variable names
        cat(names(df_cat)[i], " and ", names(df_cont)[j], "\n", "\n")
        
        #Print summary info about each group
        df_cont %>%
          rownames_to_column(var="IDs") %>% 
          left_join(df_cat %>% select_at(i) %>% rownames_to_column(var="IDs"), by="IDs") %>% 
          dplyr::select(-all_of(IDs)) %>% 
          filter_at(j, all_vars(!is.na(.))) %>%
          group_by_at(ncol(.)) %>% 
          summarise_at(j, list(n=length, min=min, median=median, mean=mean, max=max, se=std.error)) %>%
          print(n=Inf)
        
        #Test for independence, and conduct post-hoc tests
        #Avoid function stopping when a variable has only one level
        # if (sampletype=="t")
        # {
        # 	curtest <- try(t.test(df_cat[, i], df_cat[, j], alternative="two.sided", paired=F, na.action=na.omit), silent=T)
        # } else if (sampletype=="wilcox")
        # {
        # 	curtest <- try(wilcox.test(df[, i], df[, j], alternative="two.sided", paired=F, na.action=na.omit), silent=T)
        # } else
        # {
        # 	cat("Unsupported test of group differences", "\n")
        # }
        # if (inherits(curtest, "try-error"))
        # {
        # 	cat("No test of group differences run", "\n")
        # } else print(curtest)
        if (sampletype=="ANOVA")
        {
          curtest <- try(aov(df_cont[, j]~df_cat[, i], na.action=na.omit), silent=T)
          if (inherits(curtest, "try-error"))
          {
            cat("No test of group differences run", "\n")
          } else
          {
            print(anova(curtest))
            cat("\n")
            curposthoc <- try(TukeyHSD(curtest), silent=T)
            if (inherits(curposthoc, "try-error"))
            {
              cat("No post-hoc test run", "\n")
            } else
            {
              print(curposthoc)
            }
          }
        } else if (sampletype=="Kruskal")
        {
          curtest <- try(kruskal.test(df_cont[, j]~df_cat[, i], na.action=na.omit), silent=T)
          if (inherits(curtest, "try-error"))
          {
            cat("No test of group differences run", "\n")
          } else
          {
            print(curtest)
            # cat("\n")
            # curposthoc <- try(kwAllPairsDunnTest(df_cont[, j]~df_cat[, i], na.action=na.omit, p.adjust.method="holm"), silent=T)
            # if (inherits(curposthoc, "try-error"))
            # {
            # 	cat("No post-hoc test run", "\n")
            # } else
            # {
            # 	print(curposthoc)
            # }
          }
        } else
        {
          cat("Unsupported test of group differences", "\n")
        }
        # if (inherits(curtest, "try-error"))
        # {
        # 	cat("No test of group differences run", "\n")
        # } else print(curtest)
        
        cat("\n")			
      }
  } else cat("No pairs of continuous and categorical variables")
  
  #Stop sending text output to file
  sink()
  
  #Re-set the working directory to the default
  setwd(wd0)
  
  #Re-set par to what it was when the function was called
  #par(parcurrent)
}

lm_output <- function(out, wd, filename, posthocs=NULL) # linear reg. output 
{
  #Load necessary libraries
  require(car)
  require(multcomp)
  require(lmerTest)
  require(MuMIn)
  
  #Function that combines a Q-Q plot with a Q-Q line
  qq.EML <- function(resids)
  {
    qqnorm(resids)
    qqline(resids)
  }
  
  #Copy and pasted from https://github.com/aufrank/R-hacks/blob/master/mer-utils.R
  #Written by Austin F. Franks
  vif.mer <- function (fit) {
    ## adapted from rms::vif
    
    v <- vcov(fit)
    nam <- names(fixef(fit))
    
    ## exclude intercepts
    ns <- sum(1 * (nam == "Intercept" | nam == "(Intercept)"))
    if (ns > 0) {
      v <- v[-(1:ns), -(1:ns), drop = FALSE]
      nam <- nam[-(1:ns)]
    }
    
    d <- diag(v)^0.5
    v <- diag(solve(v/(d %o% d)))
    names(v) <- nam
    v
  }
  
  #Store original working directory
  wd0 <- getwd()
  #Set the working directory to where the output should go
  setwd(wd)
  #Store the current system time to use in the file name
  curtime <- format(Sys.time(), "%b%d%H%M%S")
  
  #Start sending graphic output to a PDF file
  pdf(file=paste(filename, curtime, "lm_graphs.pdf", sep="_"), width=6.5, height=8, onefile=T, title="Regression diagnostic plots", paper="letter")
  
  #Put all 2/4 plots on the same page
  #if (class(out)=="lm")
  if ("lm" %in% class(out))
  {
    par(mfrow=c(2, 2))
    plot(out, ask=F)
  } else if (class(out) %in% c("lmerMod", "lmerModLmerTest", "lmerTest", "glmerMod", "lme"))
  {
    par(mfrow=c(2, 2))
    plot(fitted(out), resid(out, type="pearson"))
    abline(h=0, lty=2)
    qq.EML(resid(out, type="pearson"))
  }
  
  #Stop sending graphic output to PDF file
  dev.off()
  
  #Start sending text output to file
  sink(file=paste(filename, curtime, "lm_output.txt", sep="_"), type="output")
  
  #Diagnostics
  #Collinearity when model has at least 2 terms
  if (length(attr(terms(out), "term.labels"))>=2)
  {
    cat("Collinearity?", "\n", "\n")
    #if (class(out)=="lm")
    #if (class(out) %in% c("lm", "lme"))
    if (length(intersect(class(out), c("lm", "lme"))))
    {
      print(vif(out))
    } else if (class(out) %in% c("lmerMod", "lmerModLmerTest", "lmerTest", "glmerMod"))
    {
      print(vif.mer(out))
    }
    cat("\n", rep("*", 30), "\n", "\n")
  }
  
  #Results (summary & anova tables)
  #if (class(out)=="lm")
  if ("lm" %in% class(out))
  {
    print(summary(out))
    cat("AICc = ", AICc(out), "\n", "\n")
    print(Anova(out, type="II"))
    # for (i in posthocs)
    # {
    # 	curcol <- noquote(i)
    # 	if (T %in% grepl(curcol, colnames(getME(out, "X"))))
    # 	{
    # 		cat("\n")
    # 		print(summary(glht(out, linfct=mcp(curcol="Tukey"))))
    # 	} else cat(noquote(i), " is not in the model.", "\n")
    # }
  } else if (class(out) %in% c("lmerMod", "lmerModLmerTest", "lmerTest", "glmerMod"))
  {
    print(summary(out))
    cat("\n", "AICc = ", AICc(out), "\n", "\n")
    print(anova(out, type="II"))
    cat("\n", rep("*", 30), "\n", "\n")
    print(drop1(out, ddf="lme4", test="Chi"))
    # for (i in posthocs)
    # {
    # 	curcol <- noquote(i)
    # 	if (T %in% grepl(curcol, colnames(getME(out, "X"))))
    # 	{
    # 		cat("\n")
    # 		print(summary(glht(out, linfct=mcp(curcol="Tukey"))))
    # 	} else cat(noquote(i), " is not in the model.", "\n")
    # }
  } else if (class(out) %in% c("lme"))
  {
    #Re-run regression w/ ML
    outML <- update(out, method="ML")
    
    print(summary(out))
    cat("\n", "AICc = ", AICc(outML), "\n", "\n")
    print(anova(out, type="marginal"))  #Pinheiro & Bates says conditional F-tests should be done on REML models.
  }
  
  #Stop sending text output to file
  sink()
  
  #Re-set the working directory to the default
  setwd(wd0)
}

# [logit_output] 
logit_output <- function(out, wd, filename, posthocs=NULL) #logistic reg. output 
{
  #Load necessary libraries
  require(car)
  require(multcomp)
  require(lmerTest)
  
  #This function was written by Elinor Lichtenberg. Following a suggestion on p. 253 of Zuur et al. 2009 (mixed models in R book): "In cases where you have a large data set, like we have in this example (1254 observations), it may be an option to extract the residuals, put them in groups of, say 10, calculate an average of the residuals per group, and use these in graphical validation plots. The groups can be based on the order of the fitted values, or on the order of a covariate."
  #This function is useful for visually assessing the fit of a logistic regression. It produces two graphs: fitted values vs. residuals, QQ plot of residuals to assess normality.  Currently the two graphs are produced in the same window, requiring that you go back to see the first one.
  aggresids.general <- function(object, predicted, resids)
  {
    temp <- data.frame(predicted, resids)
    temp <- temp[order(temp$predicted),]
    numreps <- ceiling(nrow(temp)/10)-1
    temp$index <- as.factor(c(rep(1:numreps, each=10), rep(numreps+1, (nrow(temp)-numreps*10))))
    temp2 <- aggregate(temp[, 1:2], list(temp$index), mean, na.action=na.omit)
    plot(temp2$predicted, temp2$resids)
    qqnorm(temp2$resids)
    qqline(temp2$resids)
  }
  
  #Function that combines a Q-Q plot with a Q-Q line
  qq.EML <- function(resids)
  {
    qqnorm(resids)
    qqline(resids)
  }
  
  #Copy and pasted from https://github.com/aufrank/R-hacks/blob/master/mer-utils.R
  #Written by Austin F. Franks
  vif.mer <- function (fit) {
    ## adapted from rms::vif
    
    v <- vcov(fit)
    nam <- names(fixef(fit))
    
    ## exclude intercepts
    ns <- sum(1 * (nam == "Intercept" | nam == "(Intercept)"))
    if (ns > 0) {
      v <- v[-(1:ns), -(1:ns), drop = FALSE]
      nam <- nam[-(1:ns)]
    }
    
    d <- diag(v)^0.5
    v <- diag(solve(v/(d %o% d)))
    names(v) <- nam
    v
  }
  
  # #This function tests for overdispersion. The main portion was written by Alex Forde. Elinor Lichtenberg made it a function, and modified what is outputted and returned.
  # overdisp_AF <- function(out)
  # {
  #   X2 <- sum(residuals(out, type = "pearson")^2)
  #   phi <- X2/out$df.residual
  #   p.val <- pchisq(X2, out$df.residual, df=1, lower.tail=F)
  #   return(c(phi, X2, p.val))
  # }
  
  #This function tests for overdispersion. It is from the GLMM Wiki
  overdisp_fun <- function(model) {
    rdf <- df.residual(model)
    rp <- residuals(model,type="pearson")
    Pearson.chisq <- sum(rp^2)
    prat <- Pearson.chisq/rdf
    pval <- pchisq(Pearson.chisq, df=rdf, lower.tail=FALSE)
    c(chisq=Pearson.chisq,ratio=prat,rdf=rdf,p=pval)
  }
  
  #This function calculates c-hat for computation of QAIC values
  chat <- function(out)
  {
    c <- out$deviance/out$df.residual
    return(c)
  }
  
  #Store original working directory
  wd0 <- getwd()
  #Set the working directory to where the output should go
  setwd(wd)
  #Store the current system time to use in the file name
  curtime <- format(Sys.time(), "%b%d%H%M%S")
  
  #Start sending graphic output to a PDF file
  pdf(file=paste(filename, curtime, "logit_graphs.pdf", sep="_"), width=6.5, height=8, onefile=T, title="Regression diagnostic plots", paper="letter")
  
  #Put all 2/4 plots on the same page
  if ("glm" %in% class(out))
  {
    par(mfrow=c(2, 2))
    plot(out, ask=F)
  } else if (class(out) %in% c("glmerMod", "glmerModLmerTest", "glmerTest")) #Are there cases where would be lmer rather than glmer?
  {
    par(mfrow=c(2, 2))
    plot(fitted(out), resid(out, type="pearson"))
    abline(h=0, lty=2)
    if (length(resid(out, type="pearson"))>50)
    {
      aggresids.general(out, fitted(out), resid(out, type="pearson"))
    } else
    {
      qq.EML(resid(out, type="pearson"))
    }
  }
  
  #Stop sending graphic output to PDF file
  dev.off()
  
  #Start sending text output to file
  sink(file=paste(filename, curtime, "logit_output.txt", sep="_"), type="output")
  
  #Diagnostics
  
  #Collinearity when model has at least 2 terms
  if (length(attr(terms(out), "term.labels"))>=2)
  {
    cat("Collinearity?", "\n", "\n")
    if ("glm" %in% class(out))
    {
      print(vif(out))
    } else if (class(out) %in% c("glmerMod", "glmerModLmerTest", "glmerTest"))
    {
      print(vif.mer(out))
    }
    cat("\n", rep("*", 30), "\n", "\n")
  }
  
  #Goodness of fit test (H0: Logistic regression provides adequate fit)
  if ("glm" %in% class(out))
  {
    gof.p <- with(out, 1-pchisq(deviance, df=df.residual))
    #Move below two lines outside the else once get else working
    cat("Goodness of fit p-value (rejecting H0 means poor fit): ", round(gof.p, 4))
    cat("\n", rep("*", 30), "\n", "\n")
  } else if (class(out) %in% c("lmerMod", "lmerModLmerTest", "lmerTest"))
  {
    #See https://github.com/lme4/lme4/issues/375 for info
  }
  
  
  #  if (!("quasibinomial" %in% family(out)))
  if (length(intersect(c("quasibinomial", "x.quasibinomial"), family(out)))==0)
  {
    cat("Overdispersion?", "\n", "\n")
    #overdisp <- overdisp_AF(out)  #Not sure if this function works with glmer objects - need to test
    #cat("Overdispersion parameter is ", round(overdisp[1], 4), ", which yields a chi-square value of ", round(overdisp[2], 4), " and a p-value of ", round(overdisp[3], 4))
    overdisp <- overdisp_fun(out)
    cat("Overdispersion parameter is ", round(overdisp[2], 4), ", which yields a chi-squared value of ", round(overdisp[1], 4), " and a p-value of ", round(overdisp[4], 4), " at ", round(overdisp[3], 0), " degrees of freedom")
    cat("\n", rep("*", 30), "\n", "\n")
  }
  
  #Results (summary & anova tables)
  if ("glm" %in% class(out))
  {
    print(summary(out))
    if (("quasibinomial" %in% family(out)))
    {
      require(MuMIn)
      
      cat("\nQAICc = ", QAICc(out, chat=chat(out)), "\n")
    }
    cat("\n", rep("*", 30), "\n", "\n")
    print(Anova(out, type="II"))
    # for (i in posthocs)
    # {
    # 	curcol <- noquote(i)
    # 	if (T %in% grepl(curcol, colnames(getME(out, "X"))))
    # 	{
    # 		cat("\n")
    # 		print(summary(glht(out, linfct=mcp(curcol="Tukey"))))
    # 	} else cat(noquote(i), " is not in the model.", "\n")
    # }
  } else if (class(out) %in% c("glmerMod", "glmerModLmerTest", "glmerTest"))
  {
    print(summary(out))
    if (length(intersect(c("quasibinomial", "x.quasibinomial"), family(out)))==0)
    {
      require(MuMIn)
      
      cat("\nAICc = ", AICc(out), "\n")
    }
    cat("\n", rep("*", 30), "\n", "\n")
    print(anova(out, type="II"))
    cat("\n", rep("*", 30), "\n", "\n")
    print(drop1(out, test="Chisq"))
    # for (i in posthocs)
    # {
    # 	curcol <- noquote(i)
    # 	if (T %in% grepl(curcol, colnames(getME(out, "X"))))
    # 	{
    # 		cat("\n")
    # 		print(summary(glht(out, linfct=mcp(curcol="Tukey"))))
    # 	} else cat(noquote(i), " is not in the model.", "\n")
    # }
  }
  
  #Stop sending text output to file
  sink()
  
  #Re-set the working directory to the default
  setwd(wd0)
}

# [poisson_output] poisson output
poisson_output <- function(out, wd, filename, posthocs=NULL, DEBUG=F)
{
  #Load necessary libraries
  require(car)
  require(multcomp)
  require(lmerTest)
  
  #Function that combines a Q-Q plot with a Q-Q line
  qq.EML <- function(resids)
  {
    qqnorm(resids)
    qqline(resids)
  }
  
  #Copy and pasted from https://github.com/aufrank/R-hacks/blob/master/mer-utils.R
  #Written by Austin F. Franks
  vif.mer <- function (fit) {
    ## adapted from rms::vif
    
    v <- vcov(fit)
    nam <- names(fixef(fit))
    
    ## exclude intercepts
    ns <- sum(1 * (nam == "Intercept" | nam == "(Intercept)"))
    if (ns > 0) {
      v <- v[-(1:ns), -(1:ns), drop = FALSE]
      nam <- nam[-(1:ns)]
    }
    
    d <- diag(v)^0.5
    v <- diag(solve(v/(d %o% d)))
    names(v) <- nam
    v
  }
  
  # #This function tests for overdispersion. The main portion was written by Alex Forde. Elinor Lichtenberg made it a function, and modified what is outputted and returned.
  # overdisp_AF <- function(out)
  # {
  #   X2 <- sum(residuals(out, type = "pearson")^2)
  #   phi <- X2/out$df.residual
  #   p.val <- pchisq(X2, out$df.residual, df=1, lower.tail=F)
  #   return(c(phi, X2, p.val))
  # }
  
  #This function tests for overdispersion. It is from the GLMM Wiki
  overdisp_fun <- function(model) {
    rdf <- df.residual(model)
    rp <- residuals(model,type="pearson")
    Pearson.chisq <- sum(rp^2)
    prat <- Pearson.chisq/rdf
    pval <- pchisq(Pearson.chisq, df=rdf, lower.tail=FALSE)
    c(chisq=Pearson.chisq,ratio=prat,rdf=rdf,p=pval)
  }
  
  #This function calculates c-hat for computation of QAIC values. It is from Ben Bolker's "Dealing with quasi- models in R" vignette: https://cran.r-project.org/web/packages/bbmle/vignettes/quasi.pdf.
  chat <- function(out)
  {
    c <- with(out, sum((weights * residuals^2)[weights > 0])/df.residual)
    return(c)
  }
  
  #Store original working directory
  wd0 <- getwd()
  #Set the working directory to where the output should go
  setwd(wd)
  #Store the current system time to use in the file name
  curtime <- format(Sys.time(), "%b%d%H%M%S")
  
  if (DEBUG) print("Filename set up")
  
  #Start sending graphic output to a PDF file
  pdf(file=paste(filename, curtime, "poisson_graphs.pdf", sep="_"), width=6.5, height=8, onefile=T, title="Regression diagnostic plots", paper="letter")
  
  #Put all 2/4 plots on the same page
  if ("glm" %in% class(out))
  {
    par(mfrow=c(2, 2))
    plot(out, ask=F)
  } else if (class(out) %in% c("glmerMod", "glmerModLmerTest", "glmerTest"))
  {
    par(mfrow=c(2, 2))
    plot(fitted(out), resid(out, type="pearson"))
    abline(h=0, lty=2)
    qq.EML(resid(out, type="pearson"))
  }
  
  #Stop sending graphic output to PDF file
  dev.off()
  
  if (DEBUG) print("Graphs sent to PDF")
  
  #Start sending text output to file
  sink(file=paste(filename, curtime, "poisson_output.txt", sep="_"), type="output")
  
  #Diagnostics
  
  #Collinearity when model has at least 2 terms
  if (length(attr(terms(out), "term.labels"))>=2)
  {
    cat("Collinearity?", "\n", "\n")
    if ("glm" %in% class(out))
    {
      print(vif(out))
    } else if (class(out) %in% c("glmerMod", "glmerModLmerTest", "glmerTest"))
    {
      print(vif.mer(out))
    }
    cat("\n", rep("*", 30), "\n", "\n")
  }
  
  if (DEBUG) print("Collinearity assessed")
  
  #Goodness of fit test (H0: Poisson regression provides adequate fit)
  if ("glm" %in% class(out))
  {
    gof.p <- with(out, pchisq(deviance, df.residual, lower.tail=FALSE))
    #From UCLA IDRE
    #Move below two lines out of else once else works
    cat("Goodness of fit p-value (rejecting H0 means poor fit): ", round(gof.p, 4))
    cat("\n", rep("*", 30), "\n", "\n")
  } else if (class(out) %in% c("glmerMod", "glmerModLmerTest", "glmerTest"))
  {
    #See https://github.com/lme4/lme4/issues/375 for info
    cat("Goodness of fit not yet implemented for mixed models")
    cat("\n", rep("*", 30), "\n", "\n")
  }
  
  if (DEBUG) print("GoF tested")
  
  #  if (!("quasipoisson" %in% family(out)))
  if (length(intersect(c("quasipoisson", "x.quasipoisson"), family(out)))==0)
  {
    cat("Overdispersion?", "\n", "\n")
    #overdisp <- overdisp_AF(out)  #Not sure if this function works with glmer objects - need to test
    #cat("Overdispersion parameter is ", round(overdisp[1], 4), ", which yields a chi-square value of ", round(overdisp[2], 4), " and a p-value of ", round(overdisp[3], 4))
    overdisp <- overdisp_fun(out)
    cat("Overdispersion parameter is ", round(overdisp[2], 4), ", which yields a chi-squared value of ", round(overdisp[1], 4), " and a p-value of ", round(overdisp[4], 4), " at ", round(overdisp[3], 0), " degrees of freedom")
    cat("\n", rep("*", 30), "\n", "\n")
  }
  
  #Results (summary & anova tables)
  if ("glm" %in% class(out))
  {
    print(summary(out))
    if (length(intersect(c("quasipoisson", "x.quasipoisson"), family(out)))!=0) #Need to check this line
    {
      require(MuMIn)
      
      cat("\nQAICc = ", QAICc(out, chat=chat(out)), "\n")
    }
    cat("\n", rep("*", 30), "\n", "\n")
    print(Anova(out, type="II"))
    # for (i in posthocs)
    # {
    # 	curcol <- noquote(i)
    # 	if (T %in% grepl(curcol, colnames(getME(out, "X"))))
    # 	{
    # 		cat("\n")
    # 		print(summary(glht(out, linfct=mcp(curcol="Tukey"))))
    # 	} else cat(noquote(i), " is not in the model.", "\n")
    # }
  } else if (class(out) %in% c("glmerMod", "glmerModLmerTest", "glmerTest"))
  {
    print(summary(out))
    if (length(intersect(c("quasipoisson", "x.quasipoisson"), family(out)))==0)
    {
      require(MuMIn)
      
      cat("\nAICc = ", AICc(out), "\n")
    }
    cat("\n", rep("*", 30), "\n", "\n")
    print(anova(out, type="II"))
    cat("\n", rep("*", 30), "\n", "\n")
    print(drop1(out, test="Chisq"))
    # for (i in posthocs)
    # {
    # 	curcol <- noquote(i)
    # 	if (T %in% grepl(curcol, colnames(getME(out, "X"))))
    # 	{
    # 		cat("\n")
    # 		print(summary(glht(out, linfct=mcp(curcol="Tukey"))))
    # 	} else cat(noquote(i), " is not in the model.", "\n")
    # }
  }
  
  #Stop sending text output to file
  sink()
  
  #Re-set the working directory to the default
  setwd(wd0)
}

#                                                                          ----
###############################################################################
#---------------------IMPORTED DATA SETS (exc. insects)------------------- ----
###############################################################################
#                                                                          ----
#-------------------------Plot metadata----------------------------------- ----
# Notes / instructions                                                     ----
'''Remember to keep site in the active plots data set, for the random effect.
Usable plots (since some were not monitored every year) '''

# [active_plots]                                                           ----
# Read in
active_plots <- read.csv("plot_list_byseason_31July2023.csv") 

active_plots <- reclass_year(active_plots) # {custom function}
active_plots <- site.plot(active_plots) # {custom function}

active_plots$season <- as.factor(active_plots$season)

# filter & pare down data set to what is needed 
active_plots = active_plots |>
  filter(season == c("summer")) |> # remove unusable years & non-summer seasons
  filter(use == c("TRUE")) |> # remove unusable years 
  mutate(plot_notes = NULL, plot = NULL, season = NULL) |>
  droplevels() |>
  ungroup()

# [site_info]                                                              ----
 
site_info <- read.csv("site_plot_info_SRW_22July2020.csv") # Read in

site_info$treatment <- as.factor(site_info$treatment) # reclassify data
site_info = site_info |> rename(plot = plotname)  # rename "plotname" to "plot"
site_info <- site.plot(site_info) # {custom function}
site_info <- dplyr::select(site_info, c(site.plot, treatment, availwater)) # reduce dataset
site_info$availwater <- scale(site_info$availwater) # scale availwater

# [plot_dims]:Plot dimensions (Not used)                                   ----
# SKIP ALL THE BELOW - DIDN'T END UP USING ANY OF THE DATA IN THIS SPREADSHT
#plot_dims <- read.csv("SWG_plotinfo.csv", na.strings="") # plot info (location, dimensions)

#plot_dims = plot_dims %>%
#mutate(trtmt=NULL) %>%
#mutate(who_burn=NULL) %>%
#mutate(abbreviated_plot_name=NULL) %>%
#mutate(abbreviated_site_name=NULL)

#plot_dims[,c(1:2, 5)] <- lapply(plot_dims[,c(1:2, 5)], factor)
#colnames(plot_dims)[which(names(plot_dims) == "plotname")] <- "plot"
#colnames(site_info)[which(names(site_info) == "plotname")] <- "plot"
#plot_dims = plot_dims %>%
#mutate(site = recode_factor(site, "Cross Timbers" = "CrossTimbers"))

#plot_dims$site.plot = paste5(plot_dims$site, plot_dims$plot, sep=" ", na.rm=T)
#plot_dims$site.plot = as.factor(plot_dims$site.plot)

#                                                                          ----
#--------------------------Plant/ground cover data------------------------ ----
# Notes / instructions ----
''' Notes: 
- cover: There are several NAs in canopy column. I think there was no canopy 
data taken in 2017, but I need to check. DONT NEED TO INCLUDE.
- invasives: Still need to add separate BOIS/SOHA columns to the summary table, but perhaps
this is low priority.
'''
# [landscape]: Landscape cover data                                        ----

# Read in
landscape <- read.csv("SWG_landscover_18Sept2020.csv")
landscape <- landscape |> rename(plot = "plotname") 
landscape <- site.plot(landscape) # {custom function}

landscape <- dplyr::select(landscape, c(site.plot, prop_grassy_1km, prop_seminat_1km, 
                                        prop_woody_1km, prop_dev_1km, prop_wetland_1km, 
                                        prop_ag_1km, prop_water_1km))

# [cover]: Local cover data                                                ----

cover <- read.csv("cover_2017to2022.csv") # read in

cover <- cover |>
  rename(plot = plotname) |>
  filter(season == c("summer"), # check for 1 level
         year != c("2020")) |> # check for 5 levels 
  mutate(totalground = 16) # new total ground column, set all to 16 (referencing 16 subquadrats per frame)
  
cover <- site.plot(cover) # {custom function}
cover <- reclass_year(cover) # {custom function}

# drop unneeded columns
cover <- dplyr::select(cover, c(site.plot, year, herb, wood, canopy, litter, bare, 
                                rocky, totalground))

#   [cover_s]: summarized cover data                                       ----

cover_s = cover %>% group_by(year, site.plot) %>%
  summarise(
    bare_ground = sum(bare, na.rm = TRUE),
    herb = sum(herb, na.rm = TRUE),
    litter = sum(litter, na.rm = TRUE),
    canopy = sum(canopy, na.rm = TRUE),
    wood = sum(wood, na.rm = TRUE),
    rocky = sum(rocky, na.rm = TRUE),
    total_groundcover = sum(totalground, na.rm=TRUE)) %>% 
  mutate(percent_bare = bare_ground/total_groundcover, # calculate cover proportions
         percent_herb = herb/total_groundcover,
         percent_litter = litter/total_groundcover,
         percent_wood = wood/total_groundcover,
         percent_rocky = rocky/total_groundcover,
         percent_canopy = canopy/total_groundcover) |>
  ungroup() |>
  droplevels()

# Remove superfluous columns
cover_s <- dplyr::select(cover_s, c(year, site.plot, percent_bare, percent_herb, percent_litter, percent_wood, percent_rocky, percent_canopy))

#rm(cover)

# [USDA_codes]: plant codes & traits                                       ----

# Read in
USDA_codes <- read.csv("USDA Plant List 9.29.24.csv", na.strings="NA") # USDA plant codes

# drop unneeded columns
USDA_codes <- dplyr::select(USDA_codes, c(ITIS_name, USDA_name, USDA_code, Lifespan, Nativity, Seeded.Adventive))

#USDA_codes$USDA_name <- as.factor(USDA_codes$USDA_name)
USDA_codes$Nativity <- as.factor(USDA_codes$Nativity)
USDA_codes$Lifespan <- as.factor(USDA_codes$Lifespan)
USDA_codes$USDA_code <- as.factor(USDA_codes$USDA_code)
USDA_codes$Seeded.Adventive <- as.factor(USDA_codes$Seeded.Adventive)

# Remove superfluous columns - same as select if its being bitchy
USDA_codes <- USDA_codes |>
  mutate(not_in_data=NULL) |>
  mutate(X = NULL, X.1 = NULL,  X.2 = NULL,  X.3 = NULL, X.4 = NULL, X.5 = NULL,
         X.6 = NULL, X.7 = NULL, X.8 = NULL, X.9 = NULL, X.10 = NULL, X.11 = NULL,
         X.12 = NULL, X.13 = NULL, X.14 = NULL, X.15 = NULL, X.16 = NULL) |>
  mutate_if(is.factor, funs(factor(trimws(.))))  # trim white space from factor variables

# [flowers]: floral data                                                   ----

flowers <- read.csv("inflor_cleaned_17to22.csv") # read in

flowers <- flowers |>
  filter(season == c("summer"), # check for 1 level
         year != c("2020")) |> # check for 5 levels 
  rename("USDA_code" = "plant_sp") |> # rename plant_sp column
  droplevels()

# Reclassify data
flowers <- reclass_year(flowers)

flowers$USDA_code <- as.factor(flowers$USDA_code)
flowers$inflor_count = as.numeric(flowers$inflor_count)

# rename inconsistencies
flowers = flowers %>% 
  mutate(USDA_code=recode_factor(USDA_code,"COTI3"="COREO2", "COBA2"="COREO2", .default = levels(USDA_code)))%>% #coreopsis
  mutate(USDA_code=recode_factor(USDA_code,"ERMO2"="ERIGE2", "ERPH"="ERIGE2", "ERST3"="ERIGE2", "ERTE7"="ERIGE2", .default = levels(USDA_code)))%>% #erigeron
  mutate(USDA_code=recode_factor(USDA_code,"EUDA2"="EUPHO", "CHMA15"="EUPHO", "CHGL13"="EUPHO", "CHPR6"="EUPHO", "CHNU9"="EUPHO", "EUMI5"="EUPHO","EULO2"="EUPHO", .default = levels(USDA_code)))%>% #euphorbias
  mutate(USDA_code=recode_factor(USDA_code,"OESU3"="OENOT", "OESU4"="OENOT", .default = levels(USDA_code)))%>% #oenothera
  mutate(USDA_code=recode_factor(USDA_code,"SYDI2"="SYMPH4", "SYERE"="SYMPH4", "SYSU5"="SYMPH4", .default = levels(USDA_code)))%>% #symphyotrichum
  mutate(USDA_code=recode_factor(USDA_code,"AMAM3"="AMPHI8", "AMDR"="AMPHI8", "GUSA2"="AMPHI8", "GUTE2"="AMPHI8", .default = levels(USDA_code)))%>% #broomweeds
  mutate(USDA_code=recode_factor(USDA_code,"LYALL"="LYTHR", "LYCA4"="LYTHR", .default = levels(USDA_code)))%>% #lythrum
  mutate(USDA_code=recode_factor(USDA_code,"unknown"="UNKNOWN", .default = levels(USDA_code)))

flowers <- flowers |> 
  filter(USDA_code != c("UNKNOWN")) |> #filter out unknowns
  left_join(USDA_codes, by = c("USDA_code"), relationship="many-to-many") |> #merge in traits
  droplevels()

flowers <- site.plot(flowers)
# drop unneeded columns
#flowers <- dplyr::select(flowers, c(year, site.plot, samplept, USDA_code, inflor_count))

#   [flowers_s]: summarized floral data                                    ----

# create subsets
annuals <- flowers |> 
  group_by(year, site.plot) |>
  filter(Lifespan %in% grep("A", Lifespan, value = TRUE)) |>
  summarise(
    annual_rich =  n_distinct(USDA_code, na.rm = TRUE),
    annual_abun = sum(inflor_count, na.rm = TRUE)) |>
  ungroup()

perennials <- flowers |> 
  group_by(year, site.plot) |>
  filter(Lifespan %in% grep("P", Lifespan, value = TRUE)) |>
  summarise(
    per_rich =  n_distinct(USDA_code, na.rm = TRUE),
    per_abun = sum(inflor_count, na.rm = TRUE)) |>
  ungroup()

seeded <- flowers |> 
  group_by(year, site.plot) |>
  filter(Seeded.Adventive %in% grep("S", Seeded.Adventive, value = TRUE)) |>
  summarise(
    seed_rich =  n_distinct(USDA_code, na.rm = TRUE),
    seed_abun = sum(inflor_count, na.rm = TRUE)) |>
  ungroup()

adv <- flowers |>
  group_by(year, site.plot) |>
  filter(Seeded.Adventive %in% grep("A", Seeded.Adventive, value = TRUE)) |>
  summarise(
    adv_rich =  n_distinct(USDA_code, na.rm = TRUE),
    adv_abun = sum(inflor_count, na.rm = TRUE)) |>
  ungroup()

# Create final summary table - "for each year & site.plot combination, produce one value for flr_rich & flr_abun"
flowers_s = flowers |>
  group_by(year, site.plot) |>
  summarise(
    flr_rich = n_distinct(USDA_code, na.rm = TRUE),
    flr_abun = sum(inflor_count, na.rm = TRUE)) |>
  left_join(annuals, by = c("year", "site.plot"))|>
  left_join(perennials, by = c("year", "site.plot")) |>
  left_join(seeded, by = c("year", "site.plot")) |>
  left_join(adv, by = c("year", "site.plot")) |>
  ungroup()

flowers_s[,3:12][is.na(flowers_s[,3:12])] <- 0 
flowers_s[,3:12] <- scale(flowers_s[,3:12]) 

flowers_s$flr_rich[is.na(flowers_s$flr_rich)] <- 0
flowers_s$flr_abun[is.na(flowers_s$flr_abun)] <- 0
flowers_s$annual_rich[is.na(flowers_s$annual_rich)] <- 0
flowers_s$annual_abun[is.na(flowers_s$annual_abun)] <- 0
flowers_s$per_rich[is.na(flowers_s$per_rich)] <- 0
flowers_s$per_abun[is.na(flowers_s$per_abun)] <- 0
flowers_s$seed_rich[is.na(flowers_s$seed_rich)] <- 0
flowers_s$seed_abun[is.na(flowers_s$seed_abun)] <- 0
flowers_s$adv_rich[is.na(flowers_s$adv_rich)] <- 0
flowers_s$adv_abun[is.na(flowers_s$adv_abun)] <- 0

rm(adv, annuals, perennials, seeded)

# [invasives]: invasive spp. data                                           ----

invasives <- read.csv("inv_grass_complete2017to2022_23July2023.csv", na.strings = "NA")

invasives <- site.plot(invasives)
invasives <- reclass_year(invasives)

#invasives$Prop_Bothriochloa = as.numeric(invasives$Prop_Bothriochloa)
#invasives$Prop_Sorghatrum = as.numeric(invasives$Prop_Sorghatrum)

#   [invasives_s]: invasives summary table                                 ----
invasives_s = invasives |>
  group_by(year, site.plot) |>
  filter(season == "summer") |>
  summarise(
    inv = sum(grass, na.rm = TRUE))
#BOIS = sum(#equivalent of countif - count from grass if species is bois)) 
#SOHA = same as above

invasives_s$inv <- as.numeric(invasives_s$inv)
invasives_s$inv[is.na(invasives_s$inv)] <- 0
invasives_s$inv <- scale(invasives_s$inv)

remove(invasives)

#                                                                          ----
#-------------------------Climate Data------------------------------------ ----
# f. Weather data (SUB WITH PRISM PREPROCESSING DATA)                      ----
'''BRING IN WEATHER DATA FROM PRISM PREPROCESSING 5.4.25 - SKIP ALL THIS'''

# Read in
weather <- read.csv("weather_veg_samplingdate_15July2023.csv")

# Remove superfluous columns
weather = weather |>
  mutate(site = recode_factor(site, "Cross Timbers" = "CrossTimbers")) # Rename "Cross Timbers" as "CrossTimbers"

weather <- dplyr::select(weather, c(year, season, site, plot, precip90, precip60, 
                                    avgprecip90, avgprecip60, avgtemp90, avgtemp60, 
                                    maxtemp90, maxtemp60))

# Reclassify data
weather[,1:4] <- lapply(weather[,1:4], factor)
weather$year <- year(weather$year) # [lubridate]

# Create site.plot column & reclassify
weather$site.plot = as.factor(paste5(weather$site, weather$plot, sep=" ", na.rm=T))

# Filter (remove unusable years & non-summer seasons)
weather = weather |> 
  filter(season == c("summer")) |>
  mutate(season = NULL, site = NULL, plot = NULL) |>
  droplevels()

#weather_sp18 <- weather |>
filter(season == c("spring")) |>
  filter(year == c("2018")) |>
  mutate(season = NULL, year = NULL) |>
droplevels()

# rename columns in spring 2018 dataset
colnames(weather_sp18)[which(names(weather_sp18) == "precip90")] <- "precip90_sp18" 
colnames(weather_sp18)[which(names(weather_sp18) == "precip60")] <- "precip60_sp18" 
colnames(weather_sp18)[which(names(weather_sp18) == "maxtemp90")] <- "maxtemp90_sp18" 
colnames(weather_sp18)[which(names(weather_sp18) == "maxtemp60")] <- "maxtemp60_sp18" 
colnames(weather_sp18)[which(names(weather_sp18) == "avgprecip90")] <- "avgprecip90_sp18" 
colnames(weather_sp18)[which(names(weather_sp18) == "avgprecip60")] <- "avgprecip60_sp18" 
colnames(weather_sp18)[which(names(weather_sp18) == "avgtemp90")] <- "avgtemp90_sp18" 
colnames(weather_sp18)[which(names(weather_sp18) == "avgtemp60")] <- "avgtemp60_sp18" 

#                                                                          ----
#-------------------------Other Data-------------------------------------- ----
# [func_traits]: bee functional traits NOT DONE                            ----
'''Need to get project-specific dataset from Gabby / Elinor'''

# import functional traits dataset
func_traits <- read.csv("Func_Traits_simplified.csv")

# reclassify variables as factor
func_traits[,1:7] <- lapply(func_traits[,1:7], factor) 

func_traits <- func_traits |>
  mutate(Native = recode_factor(Native, "y" = "Y"))

# merge functional traits to bee data
insects <- left_join(insects, func_traits, by="Genus_species", relationship = "many-to-many") #merge!

rm(func_traits)
#                                                                          ----
###############################################################################
#----------------------- [insects]: INSECT DATA--------------------------- ----
###############################################################################
# a. read in & merge                                                       ----

# Read in (note: no 2020 b/c of COVD)
insect17 <- read.csv("SWG_summer17_insects_clean.csv")
insect18 <- read.csv("SWG_summer18_insects_clean.csv")
insect19 <- read.csv("SWG_summer19_insects_clean.csv")
insect21 <- read.csv("SWG_summer21_insects_clean.csv")
insect22 <- read.csv("SWG_summer22_insects_clean.csv")

# merge
insects = bind_rows(insect17, insect18, insect19, insect21, insect22) 
rm(insect17, insect18, insect19, insect21, insect22) # remove old data sets

# b. clean                                                                 ----
insects = insects %>%
  mutate(Method = recode_factor(Method, "AM netting" = "netting"),
         Method = recode_factor(Method, "PM netting" = "netting"),
         Method = recode_factor(Method, "AM Netting" = "netting"),
         Method = recode_factor(Method, "PM Netting" = "netting")) %>%
  #mutate(M.F = recode_factor(M.F, "determine" = "NA"),
  # M.F = recode_factor(M.F, "na" = "NA")) %>%
  mutate(Species = recode_factor(Species, "albitarsis (Cresson)" = "albitarsis"),
         Species = recode_factor(Species, "abdominalis (Fox)" = "abdominalis"),
         Species = recode_factor(Species, "aurata (Smith)" = "aurata"),
         Species = recode_factor(Species, "bardum (Cresson)" = "bardum"),
         Species = recode_factor(Species, "birkmanni (Swenk)" = "birkmanni"),
         Species = recode_factor(Species, "brevis (Say)" = "brevis"),
         Species = recode_factor(Species, "camberella" = "cambarella"),
         Species = recode_factor(Species, "coactum (Cresson)" = "coactum"),
         Species = recode_factor(Species, "communis (Cresson)" = "communis"),
         Species = recode_factor(Species, "connexum (Cresson)" = "connexum"),
         Species = recode_factor(Species, "communis?" = "communis"),
         Species = recode_factor(Species, "condignus?" = "condignus"),
         Species = recode_factor(Species, "coquilleti (Cockerell)" = "communis"),
         Species = recode_factor(Species, "coreopsis (Robertson)" = "coreopsis"),
         Species = recode_factor(Species, "cressonii (Robertson)" = "cressonii"),
         Species = recode_factor(Species, "cressoni" = "cressonii"),
         Species = recode_factor(Species, "deflexa (Cresson)" = "deflexa"),
         Species = recode_factor(Species, "disparile (Cresson)" = "disparile"),
         Species = recode_factor(Species, "exilis (Cresson)" = "exilis"),
         Species = recode_factor(Species, "enevata" = "enavata"),
         Species = recode_factor(Species, "genilis" = "gentilis"),
         Species = recode_factor(Species, "griseocolis" = "griseocollis"),
         Species = recode_factor(Species, "hudsoniella" = "hudsoniellum"),
         Species = recode_factor(Species, "inimica (Cresson)" = "inimica"),
         Species = recode_factor(Species, "infuscata (Engel & Michez)" = "infuscata"),
         Species = recode_factor(Species, "jonesi (Cockerell)" = "jonesi"),
         Species = recode_factor(Species, "lanosa (Mocsary)" = "lanosa"),
         Species = recode_factor(Species, "mendica (Cresson)" = "mendica"),
         Species = recode_factor(Species, "metallica (Fabricius)" = "metallica"),
         Species = recode_factor(Species, "mettalica" = "metallica"),
         Species = recode_factor(Species, "mellivealvis" = "melliventris"),
         Species = recode_factor(Species, "melliventris (Cresson)" = "melliventris"),
         Species = recode_factor(Species, "mucida (Cresson)" = "mucida"),
         Species = recode_factor(Species, "obliqua (Say)" = "obliqua"),
         Species = recode_factor(Species, "montivaga?" = "montivaga"),
         Species = recode_factor(Species, "parallela (Smith)" = "parallela"),
         Species = recode_factor(Species, "paralella" = "parallela"),
         Species = recode_factor(Species, "paralellus (Say)" = "parallela"),
         Species = recode_factor(Species, "petulans (Cresson)" = "petulans"),
         Species = recode_factor(Species, "petulca (Cresson)" = "petulca"),
         Species = recode_factor(Species, "pilosifrons (Cresspon)" = "pilosifrons"),
         Species = recode_factor(Species, "policaris (Say)" = "policaris"),
         Species = recode_factor(Species, "rhodognathum (Cockerell)" = "rhodognathum"),
         Species = recode_factor(Species, "rudbeckiae (Robertson)" = "rudbeckiae"),
         Species = recode_factor(Species, "rugifrons (Smith)" = "rugifrons"),
         Species = recode_factor(Species, "tegulare?" = "tegulare"),
         Species = recode_factor(Species, "texanus (Cresson)" = "texanus"),
         Species = recode_factor(Species, "verecunda (Cresson)" = "verecunda"),
         Species = recode_factor(Species, "xanthisma (Cockerell)" = "xanthisma"),
         Genus = recode_factor(Genus, "Melissodes (Melissodes)" = "Melissodes"),
         Genus = recode_factor(Genus, "Lasioglossum (Dialictus)" = "Lasioglossum"),
         Order = recode_factor(Order, "Lipidoptera" = "Lepidoptera"),
         Order = recode_factor(Order, "Diptrera" = "Diptera"),
         Family = recode_factor(Family, "Halictidae?" = "Halictidae"),
         Family = recode_factor(Family, "Dianthidium?" = "Megachilidae"),
         #USDA_name = recode_factor(USDA_name, "new croton Oka Highway fall 17" = "Croton spp."),
         USDA_code = recode_factor(USDA_code, "OkaCroton" = "CROTO"),
         Genus_species = recode_factor(Genus_species, "Epimellisodes_petulca" = "Svastra_petulca"))

#drop_na(insects$Genus)

# Reclassify & rename data
insects[,1:8] <- lapply(insects[,1:8], factor) # make all factors (exc. year)
insects[,10:11] <- lapply(insects[,10:11], factor) # make all factors (exc. plant spp.)
# insects <- insects |> mutate(season = NULL) # KEEP SEASON

insects = insects %>%
  rename(plot = plotname) # rename plotname column
  
#insects <- site.plot(insects)
insects <- reclass_year(insects)

# Create Genus_species column and make factor
insects$Genus_species = paste5(insects$Genus, insects$Species, sep="_", na.rm=T)
insects <- insects |> mutate(Genus_species = na_if(Genus_species, '#N/A_#N/A'))
insects$Genus_species = as.factor(insects$Genus_species) 


# c. insect types                                                          ----

# Specify common type name (beetle, fly, etc.) for each order
insects <- insects |> mutate(
  insect.type = case_when(Order == "Coleoptera"~"beetle", Order == "Diptera"~"fly",
                          Order == "Lepidoptera"~"lep", Order == "Hymenoptera" & 
                          Family %in% c("Andrenidae", "Melittidae", "Apidae", 
                                        "Colletidae", "Halictidae","Megachilidae")
                          ~"bee", T~"wasp"))

insects$insect.type = as.factor(insects$insect.type) # reclassify

# d. add in plant codes                                                    ----

#insects$full_plant_sp = as.character(insects$full_plant_sp) 
#USDA_codes$USDA_name = as.character(USDA_codes$USDA_name) 

# populate USDA plant codes
insects <- left_join(insects, USDA_codes, by = c("full_plant_sp"="ITIS_name"), 
                     relationship = "many-to-many")

is.nan.data.frame(insects)
#USDA_code = na_if(USDA_code, '?'),
#USDA_code = na_if(USDA_code, 'NA'),
# Genus_species = na_if(Genus_species, '#N/A_#N/A'))

# clean up NAs
insects <- insects |>
  mutate(USDA_code = na_if(USDA_code, '#N/A'),
         USDA_code = na_if(USDA_code, '?'),
         USDA_code = na_if(USDA_code, 'NA'),
         Genus_species = na_if(Genus_species, '#N/A_#N/A'))

#Check whether any rows have full plant name but no USDA code
insects_chk <- insects |>
  filter(!is.na(full_plant_sp), is.na(USDA_code)) 
unique(insects_chk1$full_plant_sp)
rm(insects_chk1)

#Check whether any rows have a USDA code but no full plant name
insects_chk2 <- insects |> 
  filter(!is.na(USDA_code),
         is.na(full_plant_sp))
unique(insects_chk2$full_plant_sp)
rm(insects_chk2)

# rm(USDA_codes)

# e. pare down                                                             ----
insects <- dplyr::select(insects, c(SWG_ID, site.plot, Method, Order, 
                                    Family, Genus, Species, site, year, season, 
                                    insect.type, Genus_species, USDA_code, Lifespan, 
                                    Nativity, Seeded.Adventive))
#                                                                          ----
###############################################################################
#------------------ [weather]: WEATHER DATA (PRISM) ---------------------- ----
###############################################################################
# [sample_dates]                                                           ----

# read in
sample_dates <- read_csv("SWG_sampling_dates.csv")

# create site.plot column
sample_dates$site.plot <- as.factor(paste5(sample_dates$site, sample_dates$plot, sep=" "))

sample_dates <- sample_dates |>
  filter(season == "summer") |> # filter by summer season (dplyr)
  mutate(Date = mdy(date)) # format date (lubridate)

sample_dates <- dplyr::select(sample_dates, c(year, Date, site.plot))

# [PRISM]                                                                  ----

# read in
#PRISM <- read_excel("PRISM data 5.5.2025.xls", sheet = "Data")
PRISM2017 <- read_csv("PRISM_2017.csv")
PRISM2018 <- read_csv("PRISM_2018.csv")
PRISM2019 <- read_csv("PRISM_2019.csv")
PRISM2020 <- read_csv("PRISM_2020.csv")
PRISM2021 <- read_csv("PRISM_2021.csv")
PRISM2022 <- read_csv("PRISM_2022.csv")
PRISM2023 <- read_csv("PRISM_2023.csv")

PRISM <- bind_rows(PRISM2017, PRISM2018, PRISM2019, PRISM2020, PRISM2021, 
                   PRISM2022, PRISM2023) #bind together all data

rm(PRISM2017, PRISM2018, PRISM2019, PRISM2020, PRISM2021, PRISM2022, PRISM2023) 

PRISM <- PRISM |>
  rename("site.plot" = "Name", "precip_mm" = "ppt (mm)",
         "min_temp_C" = "tmin (degrees C)", "mean_temp_C" = "tmean (degrees C)",
         "max_temp_C" = "tmax (degrees C)") |>
  mutate(Date = mdy(Date))  # format date 

PRISM <- dplyr::select(PRISM, c(site.plot, Date, precip_mm, min_temp_C, mean_temp_C, max_temp_C))

PRISM$site.plot <- as.factor(PRISM$site.plot)

# Clean up site names 
PRISM <- PRISM |> 
  mutate(site.plot = recode_factor(site.plot, "Davis Cedars" = "Davis Cedarsnew"),
         site.plot = recode_factor(site.plot, "Davis Middle" = "Davis Middlenew"),
         site.plot = recode_factor(site.plot, "Davis Townne" = "Davis Townnew"),
         
         site.plot = recode_factor(site.plot, "XT 50 Acres" = "CrossTimbers 50 Acres"),
         site.plot = recode_factor(site.plot, "XT East" = "CrossTimbers East"),
         site.plot = recode_factor(site.plot, "XT Highway" = "CrossTimbers Highway"),
         site.plot = recode_factor(site.plot, "XT Middle" = "CrossTimbers Middle"),
         
         site.plot = recode_factor(site.plot, "Hagerman Ben" = "Hagerman Bennett Hill"),
         site.plot = recode_factor(site.plot, "Hagerman Cro" = "Hagerman Crow Hill"),
         site.plot = recode_factor(site.plot, "Hagerman Mea" = "Hagerman Meadow Pond"),
         site.plot = recode_factor(site.plot, "Hagerman Upp" = "Hagerman Upper"),
         
         site.plot = recode_factor(site.plot, "Lewisville A" = "Lewisville Airfield"),
         site.plot = recode_factor(site.plot, "Lewisville D" = "Lewisville Dam"),
         site.plot = recode_factor(site.plot, "Lewisville O" = "Lewisville Oakland"),
         site.plot = recode_factor(site.plot, "Lewisville W" = "Lewisville Westlake"),
         
         site.plot = recode_factor(site.plot, "Stillhouse 5" = "Stillhouse 50 Acres"),
         site.plot = recode_factor(site.plot, "Stillhouse D" = "Stillhouse Dana Peak"),
         site.plot = recode_factor(site.plot, "Stillhouse U" = "Stillhouse Union Grove"),
         site.plot = recode_factor(site.plot, "Stillhouse W" = "Stillhouse WMA5"),
         
         site.plot = recode_factor(site.plot, "Waco Area 11" = "Waco Area 111"),
         
         site.plot = recode_factor(site.plot, "Wagley Airst" = "Wagley Airstrip"),
         site.plot = recode_factor(site.plot, "Wagley Water" = "Wagley Watertrap"),
         
         site.plot = recode_factor(site.plot, "Whitney 50 A" = "Whitney 50 Acres"),
         site.plot = recode_factor(site.plot, "Whitney Floo" = "Whitney Flood"),
         site.plot = recode_factor(site.plot, "Whitney Wood" = "Whitney Woods"))

# [weather]: run 90-day calculations                                       ----
'''
Explanation:
slide() applies a rolling function that returns a numeric (double) result.
.before = days(89) creates a 90-day window ending on each date.
.index = date aligns the sliding window to actual calendar time rather than row position.
This method is memory-efficient and fast, especially for time-series data grouped 
by factors like site.plot.
'''

# Make sure date is in proper format
PRISM <- PRISM %>%
  arrange(site.plot, Date)

# Apply 90-day rolling sum using slider
PRISM <- PRISM |>
  group_by(site.plot) |>
  mutate(precip90 = slide(
    .x = precip_mm,
    .f = ~ sum(.x, na.rm = TRUE),
    .before = 89,
    .complete = FALSE,
    .index = Date)) |>
  mutate(avgprecip90 = slide(
    .x = precip_mm,
    .f = ~ mean(.x, na.rm = TRUE),
    .before = 89,
    .complete = FALSE,
    .index = Date)) |>
  mutate(avgtemp90 = slide(
    .x = mean_temp_C,
    .f = ~ mean(.x, na.rm = TRUE),
    .before = 89,
    .complete = FALSE,
    .index = Date)) |>
  mutate(maxtemp90 = slide(
    .x = max_temp_C,
    .f = ~ max(.x, na.rm = TRUE),
    .before = 89,
    .complete = FALSE,
    .index = Date)) |>
  ungroup()

# extract year into its own column
PRISM$year <- year(PRISM$Date)

PRISM$precip90 <- as.numeric(PRISM$precip90)
PRISM$avgprecip90 <- as.numeric(PRISM$avgprecip90)
PRISM$avgtemp90 <- as.numeric(PRISM$avgtemp90)
PRISM$maxtemp90 <- as.numeric(PRISM$maxtemp90)

weather <- sample_dates |>
  left_join(PRISM, by = c("site.plot", "year", "Date"))

weather <- dplyr::select(weather, c(year, site.plot, precip90, avgprecip90, 
                                    avgtemp90, maxtemp90))

rm(PRISM, sample_dates)
#                                                                          ----
###############################################################################
#-----------------DERIVED DATA SETS: ALL BEES (inc. Apis)----------------- ----
###############################################################################
# Notes & instructions                                                     ----
'''Right now, this analysis JUST looks at bees. Can change it in the future to 
include all insects.

expand.grid function creates a new (empty) data frame from all combinations 
of factor variables, EVEN THOSE THAT HAVE NO DATA ASSOCIATED.
This is DIFFERENT from what is represented in the summary datasets created above - 
THEY ARE NOT INTERCHANGABLE. 

Removed Family from analysis!'''

# [bees]: filtered insects data set                                        ----
bees <- insects |> filter(insect.type == "bee")

# [bees_s]: compiled abundance & richness data                             ----

bees_s <- insects |>
  filter(insect.type == "bee") |>
  drop_na(Genus_species) |>
  group_by(year, site.plot, Method) |> # removed Family as grouping factor
  summarise(
    bee_abun = n_distinct(SWG_ID),
    bee_rich = n_distinct(Genus_species, na.rm = FALSE)) |>
  ungroup() |>
  droplevels()

#write_csv(bees, "bees.csv")

# [bees_exp]: Create empty data frame                                      ----

# Use expand.grid to create a record for each year-site.plot-taxon combination
bees_exp = expand.grid("year" = year, "site.plot" = 
                         levels(insects$site.plot), "Method" = Method, #"Family" = Family, 
                       stringsAsFactors = TRUE)

bees_exp$year <- year(as.Date(as.character(bees_exp$year), format = "%Y")) # [lubridate]

# [beesfull]: Full effect set (bees + all fixed effects)                   ----

# Add bee abundance/richness data for each combination (expect lots of NAs here!)

# create full all-bee dataframe
beesfull <- bees_exp |>
  left_join(active_plots, by = c("year", "site.plot")) |>
  left_join(site_info, by = c("site.plot")) |>
  left_join(bees_s, by = c("year","site.plot", "Method")) |>
  left_join(landscape, by = c("site.plot")) |>
  left_join(flowers_s, by = c("year","site.plot")) |>
  left_join(weather, by = c("year","site.plot"), relationship = "many-to-many") |>
  #left_join(weather_sp18, by = c("site.plot")) |>
  left_join(cover_s, by = c("year","site.plot")) |>
  left_join(invasives_s, by = c("year","site.plot"), relationship = "many-to-many") |>
  filter(use == TRUE) |>
  mutate(use = NULL) |>
  filter(site != c("Oka")) |>
  filter(site.plot != "Stillhouse Union Grove") |> 
  droplevels()

beesfull$year <- year(as.Date(as.character(beesfull$year), format = "%Y")) # [lubridate]

beesfull$bee_abun[is.na(beesfull$bee_abun)] <- 0 # Turn abundance & richness NAs to 0s
beesfull$bee_rich[is.na(beesfull$bee_rich)] <- 0
beesfull$flr_abun[is.na(beesfull$flr_abun)] <- 0 # Turn abundance & richness NAs to 0s
beesfull$flr_rich[is.na(beesfull$flr_rich)] <- 0

# [net_beesfull]: filter by netted                                         ----
net_beesfull = beesfull |>
  filter(Method == c("netting")) |>
  mutate(Method = NULL)

#                                                                          ----
###############################################################################
#-----------------DERIVED DATA SETS: WILD BEES (NO Apis)------------------ ----
###############################################################################
# [wildbees_s]: compiled abundance & richness data                         ----

wildbees_s <- insects |>
  filter(insect.type == "bee") |>
  filter(Genus_species != "Apis_mellifera") |>
  drop_na(Genus_species) |>
  group_by(year, site.plot, Method, Family) |>
  summarise(
    bee_abun = n_distinct(SWG_ID),
    bee_rich = n_distinct(Genus_species, na.rm = FALSE)) |>
  ungroup() |>
  droplevels()

#write_csv(wildbees_s, "wildbees.csv")

# [wildbees_exp]: create empty frame for non-Apis data                     ----
wildbees_exp = expand.grid("year" = year, "site.plot" = 
                             levels(insects$site.plot), "Method" = Method, #"Family" = Family,
                           stringsAsFactors = TRUE)

# [wildbeesfull]: full effect set                                          ----
wildbeesfull <- wildbees_exp |>
  left_join(site_info, by = c("site.plot")) |>
  left_join(wildbees_s, by = c("year","site.plot", "Method", "Family")) |>
  left_join(landscape, by = c("site.plot")) |>
  left_join(flowers_s, by = c("year","site.plot")) |>
  left_join(weather, by = c("year","site.plot"), relationship = "many-to-many") |>
  #left_join(weather_sp18, by = c("site.plot")) |>
  left_join(cover_s, by = c("year","site.plot")) |>
  left_join(active_plots, by = c("year","site.plot"), relationship = "many-to-many") |>
  left_join(invasives_s, by = c("year","site.plot"), relationship = "many-to-many") |>
  filter(use == TRUE) |>
  mutate(use = NULL) |>
  filter(year != c("2020")) |> # check for 5 levels 
  filter(site.plot != c("Oka Town", "Oka Highway", "Oka River", "Oka South", 
                        "Stillhouse Union Grove")) |> 
  droplevels()

  #rm(bees_exp, wildbees_exp)
  #rm(bees_s, wildbees_s)
  
wildbeesfull$abun[is.na(wildbeesfull$bee_abun)] <- 0
wildbeesfull$rich[is.na(wildbeesfull$bee_rich)] <- 0
wildbeesfull$flr_abun[is.na(wildbeesfull$flr_abun)] <- 0 # Turn abundance & richness NAs to 0s
wildbeesfull$flr_rich[is.na(wildbeesfull$flr_rich)] <- 0
# [net_wildbeesfull]: filter by netted                                     ----

net_wildbeesfull = wildbeesfull |>
  filter(Method == c("netting")) |>
  mutate(Method = NULL)

# f. check final data set                                                  ----
'''
its a good idea to write this to CSV so you can check the dataset in excel for 
any issues 
# the file should have only one value for each abundance and richness for each 
site-year combination
# Possible duplicates to check for: Crow Hill 2022, Stillhouse Dana Peak 2018, 
Waco Area 12 2022'''

write.csv(wildbeesfull,"wildbeesfull4.9.24.csv")
write.csv(beesfull,"beesfull4.9.24.csv")

# if needed, read back in, with duplicates removed
#AM_Fam_All_bees <- read.csv("all_bees_final_check.csv")

#                                                                          ----
###############################################################################
#------------------------------LAST STEPS--------------------------------- ----
###############################################################################

#                                                                          ----