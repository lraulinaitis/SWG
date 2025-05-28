
# LIVIA RAULINAITIS
# CHAPTER 1: SWG POLLINATOR DIVERSITY ANALYSES
# See SWG Data Cleaning file for data set preparation

###############################################################################
#-------------------------------TO DO LIST-------------------------------- ----
###########################################################################
'''
* 2023 data???
* Refocus on bee functional traits (nesting & body size, sometimes sociality)
* Add back in all 2017 data - DONE
  * For 2017, all seeded plants need to be recoded as adventive (since seeding was in 2018)
* Add new columns of fixed effects - prior year climate / cover 
* 6b - create trend graphs for time vs. flower lifespan, time vs. bee functional trait

'''
#                                                                          ----
###############################################################################
#----------------------------RESEARCH QUESTIONS--------------------------- ----
###############################################################################
# a. What are the drivers of bee (floral visitors only) diversity?         ----
'''
Analysis plan:
- Shannon Weiner diversity (per plot-year combination)
- GLMMs, responses = bee abundance, bee richness, bee diversity (SW)
Notes: 


'''
# b. How does bee diversity (not community) differ among treatments? before----
#    & after treatment?                                                    ----
'''
Subquestions:
- What proportion of all floral visitors are bees?

Analysis plan:

- paired t-test (2017 vs. 2022)

Notes:
Eventually need to do separate analyses - bees including A. mellifera, and bees w/o.
'''
#                                                                          ----
###############################################################################
#------------------------ANALYSIS PRELIMINARIES--------------------------- ----
###############################################################################
# [d]: Designate & prepare final dataset                                   ----
d <- net_beesfull 

# year as integer or factor (KEEP AS INTEGER IN LUBRIDATE FORMAT)
#d$year = as.numeric(unlist(d$year))
#d$year = as.factor(d$year) 
#d$inv <- as.numeric(d$inv)

# convert NAs to 0s
d[,15:28][is.na(d[,15:28])] <- 0 
#d$bee_abun[is.na(d$bee_abun)] <- 0 
#d$bee_rich[is.na(d$bee_rich)] <- 0
#d$inv[is.na(d$inv)] <- 0

# scale everything!
#d$availwater <- scale(d$availwater)
#d$inv <- scale(d$inv)
d[,15:28] <- scale(d[,15:28]) # plant abun/richness & climate measures

# set reference levels
d <- d |> mutate(treatment = relevel(treatment, ref="control"))

#                                                                          ----
###############################################################################
#------------------------EXPLORATORY DATA ANALYSIS------------------------ ----
###########################################################################
''' Guiding questions:
1. How do plant / bee community characteristics change with time since restoration?
   (characteristics: annual/perennial, abundance/richness, seeded/adv, nativity, body size)
   - Line graphs
2. How do plant/bee community characteristics change with treatment?
   - bar charts / box plots
3. How does plant diversity influence bee diversity? - Regression plots

Ask Sean for his presentations, and use those to inform questions
'''
# a. Check data distributions                                              ----
hist(d$bee_abun) # neg binom distribution (right skewed with long tail)
hist(d$bee_rich) # neg binom
hist(d$flr_abun) # neg binom
hist(d$flr_rich) # sort of neg binom?

head()
tail()
glimpse()

# b. Trend Graphs                                                          ----
  # i. Flowers by year                                                     ----
ggplot(data = d, aes(x = year, y = flr_abun, group = year)) +
  geom_boxplot() + #geom_point()+
  labs(title = "Floral Abundance by Year", x = "Year", y = "Flower Abundance (density)") +
  theme_minimal() + theme(legend.position = "bottom")

ggplot(data = d, aes(x = year, y = flr_rich, group = year, fill = treatment)) +
  geom_boxplot() +
  #geom_point()+
  labs(title = "Floral Richness by Year",
       x = "Year",
       y = "Flower Richness (n spp.)") +
  theme_minimal()

  # ii. Bees by year                                                       ----
ggplot(data = d, aes(x = year, y = abun, group = year)) +
  geom_boxplot() +
  #geom_point()+
  labs(title = "Wild Bee Abundance by Year",
       x = "Year",
       y = "Bee Abundance (density)") +
  theme_minimal()

ggplot(data = d, aes(x = year, y = rich, group = year)) +
  geom_boxplot() +
  #geom_point()+
  labs(title = "Wild Bee Richness by Year",
       x = "Year",
       y = "Bee Richness (n spp.)") +
  theme_minimal()

  # iii. Bee abundance by richness, grouped by treatment                   ----
p <- ggplot(d, aes(x = log(abun), y = flr_rich, color = treatment)) + 
  geom_point() + 
  geom_smooth(method = 'lm')

p <- ggplot(d, aes(x = flr_abun, y = flr_rich, color = treatment)) + 
  geom_point() + 
  geom_smooth(method = 'lm')
p
p +
  facet_grid(treatment~year) # rows ~ columns

  # iv. Rich / abun by treatment (bar plots)                               ----

plot(d$rich ~ d$treatment)
plot(d$abun ~ d$treatment)
plot(d$flr_rich ~ d$treatment)
plot(d$flr_abun ~ d$treatment)

  # v. ANOVA - rich/abun by treatment                                      ----

m <- aov(data = d, rich ~ treatment)
m <- aov(data = d, abun ~ treatment)
m <- aov(data = d, flr_rich ~ treatment)
m <- aov(data = d, flr_abun ~ treatment)
summary(m)
TukeyHSD(m) # must be run on aov() model

  # vi. basic predictor graphs                                             ----
ggplot(d, aes(x = inv, y = abun)) + 
  #ggplot(d, aes(x = inv, y = abun, color = treatment)) + 
  geom_point() + 
  geom_smooth(method = 'lm') + 
  xlab("n Corners w/ inv grass presence") +
  ylab("bee abundance")
# c. Check correlation matrix                                              ----

# landscape factors
corrplot(cor(d[8:11]))
  # high - woody & grassy (neg) | dev & seminat (neg)

# weather factors
corrplot(cor(c(d[14:29])))

# current weather + availwater
d_corr2 <- select(d, c(availwater, precip90, precip60, avgprecip90, avgprecip60,
                       avgtemp90, avgtemp60, maxtemp90, maxtemp60))
corrplot(cor(d_corr2))
rm(d_corr2)

# local factors
corrplot(cor(d[30:34]))

# e. Evenness                                                              ----
#                                                                          ----
###############################################################################
#-------------------------MODELS & MODEL SELECTION------------------------ ----
###############################################################################
# a. global (final) models                                                 ----
# DUMMY MODEL                                                              ----
m <- glmmTMB(data = d, bee_abun ~ year + treatment + availwater + 
               prop_wetland_1km + prop_woody_1km + prop_ag_1km + #prop_seminat_1km + prop_grassy_1km + prop_dev_1km + 
               flr_rich + per_abun + seed_abun + per_rich + #flr_abun + + seed_rich +  adv_rich + adv_abun + 
               precip90 + avgtemp90 + #maxtemp90 + avgprecip90 + 
               percent_bare + percent_wood + #percent_canopy + # percent_rocky + percent_litter + percent_herb +
               treatment:year + flr_abun:prop_woody_1km + #treatment:precip90 +
               #precip90*availwater + 
               inv + (1|site), family = nbinom2, na.action = na.fail)

  # i. bees                                                                ----

d <- d |> 
  drop_na(prop_woody_1km, availwater, precip90, avgtemp90, percent_herb, percent_wood,
          percent_canopy, percent_litter, flr_rich, flr_abun, 
          percent_bare, treatment, year, inv)

# abun - AIC = 1527.5 (note the interaction here!)
m <- glmmTMB(data = d, bee_rich ~ flr_rich + seed_rich + per_rich + flr_abun + adv_rich +treatment 
             + year + prop_woody_1km + treatment*year + inv + (1|site), family = nbinom2)

m <- glmmTMB(bee_abun ~ flr_rich + prop_woody_1km + flr_rich*prop_woody_1km 
               + percent_litter + inv + (1|site), 
             data = d, na.action = na.fail, family = nbinom2)

# rich - AIC = 1527.5
m <- glmmTMB(bee_rich ~ flr_rich + prop_woody_1km + flr_rich*prop_woody_1km 
             + percent_litter + inv + (1|site), 
             data = d, na.action = na.fail, family = nbinom2)

  # ii. flowers                                                            ----

d <- d |> 
  drop_na(prop_grassy_1km, prop_woody_1km, availwater, precip90, avgtemp90, 
          percent_canopy, percent_litter, flr_rich, flr_abun, prop_dev_1km, 
          percent_bare, treatment, year, precip90_sp18, avgtemp90_sp18)

# flr_rich: AIC = 3022.6
m <- glmmTMB(flr_rich ~ percent_bare + inv+
                prop_woody_1km + percent_litter + 
               percent_wood + percent_canopy + year*treatment + availwater*precip90+
               (1|site), data = d, family = nbinom2, na.action = na.fail)

# flr_abun: AIC = 6812.1
m <- glmmTMB(flr_abun ~  percent_bare + inv + 
               percent_wood + percent_canopy + year*treatment + availwater*precip90+
               (1|site), data = d, family = nbinom2, na.action = na.fail)

  # iii. cover                                                             ----
m <- glmmTMB(percent_canopy ~ treatment + avgtemp90 + precip90 + prop_woody_1km 
             + percent_litter + inv + year + (1|site), 
             data = d)
#                                                                          ----
# b. model selection                                                       ----
  # i. dredge                                                              ----
# option 1: dredge (MuMIN package)
dredge <- dredge(m)
dredge # output results

rm(dredge)

  # ii. stepAIC                                                            ----
# option 2: stepAIC (MASS package)
sAIC<- stepAIC(m, scope = . ~ ., direction = "both", trace = TRUE)
summary(sAIC)

rm(sAIC)

# c. model averaging                                                       ----
avg <- summary(model.avg(dredge, subset = delta < 4, fit = TRUE))
avg

#effect_sizes = Anova(get.models(m, 1)[[1]]) # use to get effect sizes

#                                                                          ----
# d. model checking                                                        ----
  # residuals / normality                                                  ----
# check residuals
hist(resid(m)) # visualize

shapiro.test(m$residuals) # shapiro-wilkes test of normality
  # if p-value is less than 0.05, resid's are NOT normally distributed 
  # ("deviation from normality")

  # overdispersion                                                         ----
# make sure we're using the right distribution for this model and data set
check_overdispersion(m) 

  # multicollinearity                                                      ----
# check for multicollinearity - only run on models with no interactions
check_collinearity(m, component = "conditional")

  # summarize output                                                       ----
# summarize model output
Anova(m, type="III")
summary(m)

#                                                                          ----
###############################################################################
###############################################################################
###############################################################################
#---------------------GRAVEYARD---------------------------------------- ####
#############Not related - get unique codes & spp names################
# Not related - get unique codes & spp names

# summarise specimen data in new sheet
flr_code_unique <- flowers %>% group_by(plant_sp) %>%
  summarise(
    abun = n_distinct(unique_ID))

flr_code_unique

# Assuming 'df' is your data frame and 'Column_Name' is the column you want to extract unique values from

# Using the unique() function to get unique values from the specified column
unique_values <- sort(unique(full_flr_data$plant_sp))

# Printing the unique values
print(unique_values)


flowers_comb_MC <- flowers |>
  left_join(site_info, by = c("site.plot")) |>
  left_join(active_plots, by = c("year","site.plot"), relationship = "many-to-many") |>
  filter(use == TRUE)

# Future automated program infrastructure----
# got a little carried away here and started to edit my original code to match this... this is a later project...
d_dict <- c(
  "YAN" = AM_Fam_Net_bees, 
  "YAP" = AM_Fam_PT_bees,
  "YAB" = AM_Fam_All_bees,
  "YUN" = c(...) "AM_Funct_Net_bees",
  "YUP" = "AM_Funct_PT_bees",
  "YUB" = "AM_Funct_All_bees",
  "NAN" = "NoAM_Fam_Net_bees",
  "NAP" = "NoAM_Fam_PT_bees",
  "NAB" = "NoAM_Fam_All_bees",
  "NUN" = "NoAM_Func_Net_bees",
  "NUP" = "NoAM_Func_PT_bees",
  "NUB" = "NoAM_Func_All_bees",
)


# by taxonomic group

# Function only summarizes by family for now, add functional groups later
# still need to adjust filter elements of function and test run it.
data_summary <- function(insects_full) {
  
  # assign insect.type to each entry
  insects_full <- insects_full |>
    mutate(insect.type=case_when(Order=="Coleoptera"~"beetle",
                                 Order=="Diptera"~"fly",
                                 Order=="Lepidoptera"~"lep",
                                 Order=="Hymenoptera" & Family %in% c("Andrenidae", "Melittidae",
                                                                      "Apidae", "Colletidae", "Halictidae","Megachilidae")~"bee", T~"wasp"))
  
  
  cat("Select a taxonomic group by which to compile richness / abundance data:", 
      paste(seq_along(levels(insects_full$insect.type)), levels(insects_full$insect.type), sep = ". "), sep = "\n")
  cat("Or enter 6 for all insects.\n")
  taxon_sel_num <- as.numeric(toupper(readline("Selection: ")))
  
  # validate input and extract associated factor
  if (taxon_sel_num != 6 && taxon_sel_num <= length(levels(insects_full$insect.type))) {
    taxon_sel <- levels(insects_full$insect.type)[taxon_sel_num] # Use # as index to get insect.type
  } else if (taxon_sel_num > 6) {
    print("Invalid.")
  }
  # use selection to determine if we need to filter dataset by insect type
  # need to merge this conditional with the AM conditional to make sure they output a single dataset
  if (taxon_sel_num != 6) {  
    insects_filtered <- insects_full |>
      filter(insect.type == taxon_sel)
  }
  
  AM <- toupper(readline("Exclude honey bees from dataset? (Y/N)"))
  if (AM[1] == "Y") {
    insects_filtered <- insects_full |>
      filter(Genus_species != "Apis_mellifera")  # Fixed: Use "Apis_mellifera" as a string
  }
  
  # Compile insect abundance/richness by FAMILY
  insects_summary <- insects_full %>% 
    # Remove the filter line if you want to consider all taxonomic groups
    group_by(year, site.plot, Order, Method, Family) %>%
    summarise(
      abun = n_distinct(SWG_ID),
      richness = n_distinct(Genus_species, na.rm = FALSE)) %>%
    ungroup() %>%
    droplevels()
  
  return(insects_summary)
  
}

# call function
data_summary <- data_summary(insects_full)

# Unused data cleaning functions ----

#Make list of taxa netted on flowers
netted <- insects17 %>% 
  filter(grepl("netting", Method)) %>% 
  distinct(taxon) %>% 
  pull(taxon)

#Look at taxa not found on flowers
notnetted <- insects17 %>% 
  filter(!(taxon %in% netted)) %>% 
  droplevels() %>% 
  distinct(taxon) %>% 
  arrange(taxon) %>% 
  pull(taxon)

# other chunks ----
# compile bee-only abundance/richness by FAMILY
bee_fam_summary <- insects_full %>% 
  filter(insect.type == c("bee")) %>%
  group_by(year, site.plot, Order, Method, Family) %>%
  summarise(
    abun = n_distinct(SWG_ID),
    richness = n_distinct(Genus_species, na.rm = FALSE)) %>%
  ungroup() %>%
  droplevels()

# compile NON-APIS bee abundance/richness by SPECIES
NoAM_bee_spp_summary <- bee_spp_summary %>% 
  filter(Genus_species != c("Apis_mellifera")) %>%
  droplevels()

# compile NON-APIS bee abundance/richness by FAMILY
NoAM_bee_fam_summary <- insects_full %>% 
  filter(Genus_species != c("Apis_mellifera")) %>%
  filter(insect.type == c("bee")) %>%
  group_by(year, site.plot, Order, Method, Family) %>%
  summarise(
    abun = n_distinct(SWG_ID),
    richness = n_distinct(Genus_species, na.rm = FALSE)) %>%
  ungroup() %>%
  droplevels()

# MOVE TO GRAVEYARD----
# issue with plot numbers for 50-acre plots. Cleaned up in excel and read back in.

# Other common misspellings to keep an eye out for:   
# Lewisville - "Lewsville"
# Meadow pond - "meadowpond" / "Meadowpond"
# Bennett Hill - "Benett Hill" / "Bennet Hill"
# Watertrap - "Water trap" / "Water Trap"

#                                                                          ----
###############################################################################
#------------------2. IMPORT DATA SETS------------------------------------ ----
###############################################################################
# a. Plot list ----
'''Remember to keep site in this data set, for the random effect.'''

active_plots <- read.csv("plot_list_byseason_31July2023.csv", na.strings = "NA") # usable plots (since some were not monitored every year)

# Reclassify data
active_plots$year <- year(as.Date(as.character(active_plots$year), format = '%Y')) # [lubridate]

active_plots$site <- as.factor(active_plots$site)
active_plots$plot <- as.factor(active_plots$plot)
active_plots$season <- as.factor(active_plots$season)

# filter & pare down data set to what is needed 
active_plots = active_plots |>
  filter(season == c("summer")) |> # remove unusable years & non-summer seasons
  filter(use == c("TRUE")) |> # remove unusable years 
  mutate(site = recode_factor(site, "Cross Timbers" = "CrossTimbers")) |> # remove unneeded columns
  mutate(site.plot = as.factor(paste5(site, plot, sep=" ", na.rm=T))) |> # create site.plot column
  mutate(plot_notes = NULL, plot = NULL, season = NULL, use = NULL) |>
  ungroup()

# b. Plot metadata ----

# Read in 
site_info <- read.csv("site_plot_info_SRW_22July2020.csv", na.strings = "NA", fileEncoding = "UTF-8-BOM") # temperature, burn info

site_info <- dplyr::select(site_info, c(site, plotname, treatment, availwater)) # pare down dataset

# Reclassify data
site_info$treatment <- as.factor(site_info$treatment)

# Rename "Cross Timbers" as "CrossTimbers"
site_info = site_info |>
  rename(plot = plotname) |> # rename "plotname" to "plot"
  mutate(site = recode_factor(site, "Cross Timbers" = "CrossTimbers")) |>
  mutate(site.plot = as.factor(paste5(site, plot, sep=" ", na.rm=T))) |> # Create site.plot columns & reclassify
  mutate(site = NULL, plot = NULL)  # remove site and plot columns


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

# c. Landscape data ----
''' notes about this dataset:
1. This is an old dataset that only extends to 2020. Need to get an updated dataset.'''

# Read in
landscape <- read.csv("SWG_landscover_18Sept2020.csv", na.strings = "", fileEncoding = "UTF-8-BOM")

landscape <- dplyr::select(landscape, c(site, plotname, prop_grassy_1km, prop_seminat_1km, prop_woody_1km, prop_dev_1km, prop_wetland_1km, prop_ag_1km, prop_water_1km))

landscape <- landscape |>
  rename(plot = "plotname") |>
  mutate(site = recode_factor(site, "Cross Timbers" = "CrossTimbers")) |>
  filter(site != c("Oka")) |> # Filter out Oka and Union Grove plots
  filter(plot != c("Union Grove")) |>
  mutate(site.plot = as.factor(paste5(site, plot, sep=" ", na.rm=T))) |>
  mutate(site=NULL, plot = NULL) # remove site and plot columns

# in case select is being bitchy
mutate(prop_seminat_500m = NULL, prop_ag_500m = NULL, prop_woody_500m = NULL, 
       prop_wetland_500m = NULL, prop_dev_500m = NULL, prop_grassy_500m = NULL,
       prop_water_500m = NULL) |>
  mutate(prop_seminat_750m = NULL, prop_ag_750m = NULL, prop_woody_750m = NULL, 
         prop_wetland_750m = NULL, prop_dev_750m = NULL, prop_grassy_750m = NULL,
         prop_water_750m = NULL) |>
  mutate(prop_seminat_1km = NULL, prop_ag_1km = NULL, 
         prop_wetland_1km = NULL, prop_dev_1km = NULL, prop_grassy_1km = NULL,
         prop_water_1km = NULL) |>
  mutate(NP_500m = NULL, TE_500m = NULL, NP_750m = NULL, TE_750m = NULL, NP_1km = NULL, TE_1km = NULL)

# d. Local cover data ----
''' Notes about this data set:
1. There are several NAs in canopy column. I think there was no canopy data taken 
in 2017, but I need to check. DONT NEED TO INCLUDE.
'''

# read in
cover <- read.csv("cover_2017to2022.csv", na.strings = "NA", fileEncoding = "UTF-8-BOM")

# drop unneeded columns
cover <- dplyr::select(cover, c(site, plotname, transectID, position, herb, wood, canopy, litter, bare, rocky, year, season))

cover <- cover |>
  rename(plot = plotname) |>
  filter(season == c("summer"), # check for 1 level
         #year != c("2020"), # check for 4 levels (2018, 2019, 2021, 2022)
         site != c("Oka"), # check for 10 levels
         plot != c("Union Grove")) |>
  mutate(site = recode_factor(site, "Cross Timbers" = "CrossTimbers")) |> # Recode "Cross Timbers" as "CrossTimbers"
  mutate(totalground = 16) |> # new total ground column, set all to 16 (referencing 16 subquadrats per frame)
  mutate_if(is.factor, funs(factor(trimws(.)))) |> # trim white space from factor variables - b/c of "Wagley" and "Wagley " 
  mutate(site.plot = as.factor(paste5(site, plot, sep=" ", na.rm=T))) |> # new site.plot column
  droplevels()

# Reclassify data
#cover[,c(8:13)] <- as.numeric(unlist(cover[,c(8:13)])) 
#cover$year <- as.factor(cover$year)

# Create summary table
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

# e. Plant codes & traits ----

# Read in
USDA_codes <- read.csv("USDA Plant List 9.29.24.csv", na.strings="NA") # USDA plant codes

# drop unneeded columns
USDA_codes <- dplyr::select(USDA_codes, c(ITIS_name, USDA_name, USDA_code, Lifespan, Nativity, Seeded.Adventive, Flower.Color))

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

# f. Floral data ----

# Read in
flowers <- read.csv("inflor_cleaned_17to22.csv", na.strings = "NA", fileEncoding = "UTF-8-BOM")

# drop unneeded columns
flowers <- dplyr::select(flowers, c(site, plot, transectID, position, plant_sp, inflor_count, year, season))

#does the same as above if select is being bitchy
flowers <- flowers |>
  mutate(unique_ID = NULL, date = NULL, observer = NULL, transectID = NULL,
         start_time = NULL, position = NULL, field_plant_sp = NULL, true_plant_sp = NULL,
         notes = NULL, samplept = NULL)

# Mutate function: create new site.plot column (combining content from site and plot columns)
flowers <- flowers |>
  filter(season == c("summer"), # check for 1 level
         #year != c("2020"), # check for 4 levels (2018, 2019, 2021, 2022)
         site != c("Oka"), # check for 10 levels
         plot != c("Union Grove")) |>
  rename("USDA_code" = "plant_sp") |>
  mutate(site = recode_factor(site, "Cross Timbers" = "CrossTimbers")) |> # Recode "Cross Timbers" as "CrossTimbers"
  mutate_if(is.factor, funs(factor(trimws(.)))) |> # trim white space from factor variables - b/c of "Wagley" and "Wagley "
  mutate(site.plot = as.factor(paste5(site, plot, sep=" ", na.rm=T))) |> # new site.plot column
  droplevels()

# Reclassify data
flowers$year = year(flowers$year) 
flowers$USDA_code <- as.factor(flowers$USDA_code)
flowers$inflor_count = as.numeric(flowers$inflor_count)

# rename inconsistencies
flowers = flowers %>% 
  mutate(USDA_code=recode_factor(USDA_code, "COTI3" = "COREO2", "COBA2" = "COREO2", .default = levels(USDA_code)))%>% #coreopsis
  mutate(USDA_code=recode_factor(USDA_code, "ERMO2" = "ERIGE2", "ERPH"="ERIGE2", "ERST3"="ERIGE2", "ERTE7"="ERIGE2", .default = levels(USDA_code)))%>% #erigeron
  mutate(USDA_code=recode_factor(USDA_code, "EUDA2" = "EUPHO", "CHMA15"="EUPHO", "CHGL13"="EUPHO", "CHPR6"="EUPHO", "CHNU9"="EUPHO", "EUMI5"="EUPHO","EULO2"="EUPHO", .default = levels(USDA_code)))%>% #euphorbias
  mutate(USDA_code=recode_factor(USDA_code, "OESU3" = "OENOT", "OESU4"="OENOT", .default = levels(USDA_code)))%>% #oenothera
  mutate(USDA_code=recode_factor(USDA_code, "SYDI2" = "SYMPH4", "SYERE"="SYMPH4", "SYSU5"="SYMPH4", .default = levels(USDA_code)))%>% #symphyotrichum
  mutate(USDA_code=recode_factor(USDA_code, "AMAM3" = "AMPHI8", "AMDR"="AMPHI8", "GUSA2"="AMPHI8", "GUTE2"="AMPHI8", .default = levels(USDA_code)))%>% #broomweeds
  mutate(USDA_code=recode_factor(USDA_code, "LYALL" = "LYTHR", "LYCA4"="LYTHR", .default = levels(USDA_code)))%>% #lythrum
  mutate(USDA_code=recode_factor(USDA_code, "unknown" = "UNKNOWN", .default = levels(USDA_code)))

#flowers_noUNK <- flowers %>% filter(! plant_sp %in% c( "UNKNOWN")) #filter out unknowns

# merge in plant traits
flowers <- flowers |>
  left_join(USDA_codes, by = c("USDA_code"), relationship = "many-to-many") 

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


# f. Weather data ----
'''BRING IN WEATHER DATA FROM PRISM PREPROCESSING 5.4.25 - SKIP ALL THIS'''
#weather$year <- as.factor(weather$year)

# Read in
weather <- read.csv("weather_veg_samplingdate_15July2023.csv", na.strings = "NA", fileEncoding = "UTF-8-BOM")

# Remove superfluous columns
weather = weather |>
  mutate(siteplot=NULL,
         date=NULL,
         site = recode_factor(site, "Cross Timbers" = "CrossTimbers")) # Rename "Cross Timbers" as "CrossTimbers"

weather <- dplyr::select(weather, c(year, season, site, plot, precip90, precip60, avgprecip90, avgprecip60,
                                    avgtemp90, avgtemp60, maxtemp90, maxtemp60))
weather <- weather |>
  mutate(precip14 = NULL, avgprecip14 = NULL, avgtemp14 = NULL, maxtemp14 = NULL, gdd14 = NULL) |>
  mutate(precip30 = NULL, avgprecip30 = NULL, avgtemp30 = NULL, maxtemp30 = NULL, gdd30 = NULL) |>
  mutate(precip60 = NULL, avgprecip60 = NULL, avgtemp60 = NULL, maxtemp60 = NULL, gdd60 = NULL) |>
  mutate(precip120 = NULL, avgprecip120 = NULL, avgtemp120 = NULL, maxtemp120 = NULL, gdd120 = NULL) |>
  mutate(avgprecip90 = NULL, maxtemp90 = NULL, gdd90 = NULL) |>
  mutate(precip150 = NULL, avgprecip150 = NULL, avgtemp150 = NULL, maxtemp150 = NULL, gdd150 = NULL) |>
  mutate(precip180 = NULL, avgprecip180 = NULL, avgtemp180 = NULL, maxtemp180 = NULL, gdd180 = NULL) |>
  mutate(precip365 = NULL, avgprecip365 = NULL, avgtemp365 = NULL, maxtemp365 = NULL, gdd365 = NULL)

# Reclassify data
weather[,1:4] <- lapply(weather[,1:4], factor)
weather$year <- year(weather$year) # [lubridate]

# Create site.plot column & reclassify
weather$site.plot = paste5(weather$site, weather$plot, sep=" ", na.rm=T)
weather$site.plot = as.factor(weather$site.plot)

weather = weather |> mutate(site=NULL, plot=NULL)

# Filter (remove unusable years & non-summer seasons)
weather = weather |> 
  filter(season == c("summer")) |>
  #filter(year != c("2017")) |>
  mutate(season = NULL) |>
  droplevels()

#weather_sp18 <- weather |>
filter(season == c("spring")) |>
  filter(year == c("2018")) |>
  mutate(season = NULL) |>
  mutate(year = NULL)
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


# g. Insect Data ----

# Read in (all separate files, need to be combined). Includes BOTH netting & pan trap data.
# Note: no insect data were collected summer 2020 b/c of covid
insect17 <- read.csv("SWG_summer17_insects_clean.csv", na.strings = "NA", fileEncoding = "UTF-8-BOM")
insect18 <- read.csv("SWG_summer18_insects_clean.csv", na.strings = "NA", fileEncoding = "UTF-8-BOM")
insect19 <- read.csv("SWG_summer19_insects_clean.csv", na.strings = "NA", fileEncoding = "UTF-8-BOM")
insect21 <- read.csv("SWG_summer21_insects_clean.csv", na.strings = "NA", fileEncoding = "UTF-8-BOM")
insect22 <- read.csv("SWG_summer22_insects_clean.csv", na.strings = "NA", fileEncoding = "UTF-8-BOM")

# Merge to one dataset
insects = bind_rows(insect17, insect18, insect19, insect21, insect22) #bind together all specimen data
rm(insect17, insect18, insect19, insect21, insect22) # remove old data sets

insects <- dplyr::select(insects, c(site.plot, Method, Order, Family, Genus, Species, 
                                    site, plotname, year, season, insect.type, full_plant_sp, SWG_ID))

# Reclassify & rename data
insects[,1:8] <- lapply(insects[,1:8], factor) # make all columns factors (exc. year)
insects[,10:11] <- lapply(insects[,10:11], factor) # make all columns factors (exc. plant spp.)
insects$year <- year(insects$year)

insects = insects %>%
  rename(plot = plotname) |> # rename plotname column
  mutate(site = recode_factor(site, "Cross Timbers" = "CrossTimbers"))  # rename all "Cross Timbers" values

# Create Genus_species column and make factor
insects$Genus_species = as.factor(paste5(insects$Genus, insects$Species, sep="_", na.rm=T)) 
#insects$year = as.factor(insects$year) 


# Specify common type name (beetle, fly, etc.) for each order
insects <- insects %>%
  mutate(insect.type=case_when(Order=="Coleoptera"~"beetle",
                               Order=="Diptera"~"fly",
                               Order=="Lepidoptera"~"lep",
                               Order=="Hymenoptera" & Family %in% c("Andrenidae", "Melittidae",
                                                                    "Apidae", "Colletidae", "Halictidae","Megachilidae")~"bee", T~"wasp"))

# Reclassify as type
#insects$insect.type = as.factor(insects$insect.type) 
# insects <- insects |> mutate(season = NULL) # KEEP SEASON

#insects$full_plant_sp = as.character(insects$full_plant_sp) 
#USDA_codes$USDA_name = as.character(USDA_codes$USDA_name) 

# populate USDA plant codes
insects <- left_join(insects, USDA_codes, by = c("full_plant_sp"="ITIS_name"), relationship = "many-to-many")

insects <- dplyr::select(insects, c(SWG_ID, site.plot, Method, Order, 
                                    Family, Genus, Species, site, year, season, 
                                    insect.type, Genus_species, USDA_code, Lifespan, 
                                    Nativity, Seeded.Adventive))


#rm(USDA_codes)

# h. Invasive spp. data ----
'''Still need to add separate BOIS/SOHA columns to the summary table, but perhaps
this is low priority.'''

invasives <- read.csv("inv_grass_complete2017to2022_23July2023.csv", na.strings = "NA")
invasives <- invasives |> 
  mutate(site = recode_factor(site, "Cross Timbers" = "CrossTimbers"))  # Recode "Cross Timbers" as "CrossTimbers"

invasives$site.plot = as.factor(paste5(invasives$site, invasives$plot, sep=" ", na.rm=T)) 
# check for 38 levels - only 34, but that's because sites were added

#invasives$year = as.factor(invasives$year)
#invasives$Prop_Bothriochloa = as.numeric(invasives$Prop_Bothriochloa)
#invasives$Prop_Sorghatrum = as.numeric(invasives$Prop_Sorghatrum)

# Create summary table
invasives_s = invasives |>
  group_by(year, site.plot) |>
  filter(season == "summer") |>
  summarise(
    inv = sum(grass, na.rm = TRUE))
#BOIS = sum(#equivalent of countif - count from grass if species is bois)) 
#SOHA = same as above

remove(invasives)

# i. Functional traits (bees) NOT DONE ----
'''Need to get project-specific dataset from Gabby / Elinor'''

# import functional traits dataset
func_traits <- read.csv("Func_Traits_simplified.csv", na.strings = "NA", fileEncoding = "UTF-8-BOM")

# reclassify variables as factor
func_traits[,1:7] <- lapply(func_traits[,1:7], factor) 

# merge functional traits to bee data
insects <- left_join(insects, func_traits, by="Genus_species", relationship = "many-to-many") #merge!

# compile NON-APIS bee abundance/richness by FAMILY - don't think I need this yet...
# NoAM_bee_nesters <- bees_full %>% 
filter(Genus_species != c("Apis_mellifera")) %>%
  group_by(year, site.plot, Order, Method, Nesting.Behavior) %>%
  summarise(
    abun = n_distinct(SWG_ID),
    richness = n_distinct(Genus_species, na.rm = FALSE)) %>%
  ungroup() %>%
  droplevels()

#                                                                          ----