# PRISM (Weather / Climate) Data Preprocessing for SWG analysis
# Livia Raulinaitis
# May 2025
#                                                                          ----
###############################################################################
#-------------------1. PACKAGES, FUNCTIONS, & GLOBAL VARIABLES------------ ----
###############################################################################
#                                                                          ----
# a. Packages ----
library(tidyverse) # includes lubridate, dplyr, purrr, readxl
library(janitor)

library(timeDate)
library(pollen)

library(data.table) # solution 2
library(slider) # solution 3
# b. Functions ----
# Paste strings & remove NAs (paste5) (DO NOT EDIT): 
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
#                                                                          ----
###############################################################################
#------------------2. IMPORT & REDUCE DATA ------------------------------- ----
###############################################################################
#                                                                          ----
'''
Issues with the PRISM dataset right now - only has years 2023 & 2024. I must have
made a mistake when downloading the dataset. Redo.
'''
# a. sample dates (Done)----

# read in
sample_dates <- read_csv("SWG_sampling_dates.csv")

# create site.plot column
sample_dates$site.plot <- as.factor(paste5(sample_dates$site, sample_dates$plot, sep=" "))

sample_dates <- sample_dates |>
  filter(season == "summer") |> # filter by summer season (dplyr)
  mutate(Date = mdy(date)) # format date (lubridate)

sample_dates <- dplyr::select(sample_dates, c(year, Date, site.plot))

# b. PRISM (Done) ----

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
# c. clean up site names (PRISM)----

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

#                                                                          ----
###############################################################################
#------------------3. RUN 90-DAY CALCULATIONS----------------------------- ----
###############################################################################
#                                                                          ----
# SOLUTION 1 - WORKS----
  # explanation ----
'''
Explanation:
- The data is grouped by site.plot to ensure calculations are plot-specific. 
- For each row, it computes the sum of precip_mm over the 90-day window ending 
on the current date.
- row_number() and sapply() are used to loop through each row within each group.

Note:
If your dataset is large, this approach may be slow. A more efficient approach 
would use a rolling join or cumulative sums, possibly with data.table or optimized 
functions in slider or zoo.
'''
  # code ----
# Arrange data by site.plot and date
PRISM <- PRISM %>%
  arrange(site.plot, Date)

#daily gdd
PRISM$daily_gdd=gdd(tmax = PRISM$max_temp_C, tmin = PRISM$min_temp_C, tbase=10, tbase_max = 90, type = "C")
# ERROR - currently doesn't reset by year or by site.plot

# Calculate the 90-day rolling sum
PRISM <- PRISM %>%
  group_by(site.plot) %>%
  mutate(
    precip90 = sapply(row_number(), function(i) {
      current_date <- Date[i]
      sum(precip_mm[Date >= current_date - days(89) & Date <= current_date], na.rm = TRUE)
    }),
    avgprecip90 = sapply(row_number(), function(i) {
      current_date <- Date[i]
      mean(precip_mm[Date >= current_date - days(89) & Date <= current_date], na.rm = TRUE)
    }),
    avgtemp90 = sapply(row_number(), function(i) {
      current_date <- Date[i]
      mean(mean_temp_C[Date >= current_date - days(89) & Date <= current_date], na.rm = TRUE)
    }),
    maxtemp90 = sapply(row_number(), function(i) {
      current_date <- Date[i]
      sum(precip_mm[Date >= current_date - days(89) & Date <= current_date], na.rm = TRUE)
    })
  ) %>%
  ungroup()

#gdd120=max(gdd(tmax=tmax..degrees.C., tmin=tmin..degrees.C., tbase=10, tbase_max=120, type = "C")),date=x)

# SOLUTION 3 - PREFERRED ----
  #explanation ----
'''
Explanation:
slide_dbl() applies a rolling function that returns a numeric (double) result.
.before = days(89) creates a 90-day window ending on each date.
.index = date aligns the sliding window to actual calendar time rather than row position.
This method is memory-efficient and fast, especially for time-series data grouped 
by factors like site.plot.
'''

  # code ----
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

weather <- dplyr::select(weather, c(year, site.plot, precip90, avgprecip90, avgtemp90, maxtemp90))

rm(PRISM, sample_dates)
# SEAN'S SOLUTION ----


func <- function(x) {

    Data14=raw_data %>% 
      filter(between(raw_data$Date, x - days(14)), x)
    
    Data30=PRISM %>% 
      filter(between(PRISM$Date, x - days(30)), x)
    
    Data60=raw_data %>% 
      select(siteplot,Date,ppt..mm.,tmin..degrees.C.,tmean..degrees.C.,tmax..degrees.C.) %>% 
      filter(between(as.Date(raw_data$Date, format='%m/%d/%Y'), 
                     (as.Date(x, format='%m/%d/%Y') - days(60)), 
                     as.Date(x, format='%m/%d/%Y')))
    
    Data90=raw_data %>% 
      select(siteplot,Date,ppt..mm.,tmin..degrees.C.,tmean..degrees.C.,tmax..degrees.C.) %>% 
      filter(between(as.Date(raw_data$Date, format='%m/%d/%Y'), 
                     (as.Date(x, format='%m/%d/%Y') - days(90)), 
                     as.Date(x, format='%m/%d/%Y')))
    
    Data120=raw_data %>% 
      select(siteplot,Date,ppt..mm.,tmin..degrees.C.,tmean..degrees.C.,tmax..degrees.C.) %>% 
      filter(between(as.Date(raw_data$Date, format='%m/%d/%Y'), 
                     (as.Date(x, format='%m/%d/%Y') - days(120)), 
                     as.Date(x, format='%m/%d/%Y')))
    
    Data150=raw_data %>% 
      select(siteplot,Date,ppt..mm.,tmin..degrees.C.,tmean..degrees.C.,tmax..degrees.C.) %>% 
      filter(between(as.Date(raw_data$Date, format='%m/%d/%Y'), 
                     (as.Date(x, format='%m/%d/%Y') - days(150)), 
                     as.Date(x, format='%m/%d/%Y')))
    
    Data180=raw_data %>% 
      select(siteplot,Date,ppt..mm.,tmin..degrees.C.,tmean..degrees.C.,tmax..degrees.C.) %>% 
      filter(between(as.Date(raw_data$Date, format='%m/%d/%Y'), 
                     (as.Date(x, format='%m/%d/%Y') - days(180)), 
                     as.Date(x, format='%m/%d/%Y')))
    
    Data365=raw_data %>% 
      select(siteplot,Date,ppt..mm.,tmin..degrees.C.,tmean..degrees.C.,tmax..degrees.C.) %>% 
      filter(between(as.Date(raw_data$Date, format='%m/%d/%Y'), 
                     (as.Date(x, format='%m/%d/%Y') - days(365)), 
                     as.Date(x, format='%m/%d/%Y')))
    
    Weather14=Data14 %>% group_by(site.plot)%>%
    summarise(precip14= sum(ppt..mm.),
              avgprecip14=mean(ppt..mm.),
              avgtemp14=mean(tmean..degrees.C.),
              maxtemp14=max(tmax..degrees.C.),
              gdd14=max(gdd(tmax=tmax..degrees.C., tmin=tmin..degrees.C., tbase=10, 
                            tbase_max=30, type = "C")), date=x)
    
    
    Weather30=Data30 %>% group_by(siteplot)%>%
      summarise(precip30= sum(ppt..mm.),
                avgprecip30=mean(ppt..mm.),
                avgtemp30=mean(tmean..degrees.C.),
                maxtemp30=max(tmax..degrees.C.),
                gdd30=max(gdd(tmax=tmax..degrees.C., tmin=tmin..degrees.C., tbase=10, tbase_max=30, type = "C")),
                date=x)
    
    
    Weather60=Data60 %>% group_by(siteplot)%>%
      summarise(precip60= sum(ppt..mm.),
                avgprecip60=mean(ppt..mm.),
                avgtemp60=mean(tmean..degrees.C.),
                maxtemp60=max(tmax..degrees.C.),
                gdd60=max(gdd(tmax=tmax..degrees.C., tmin=tmin..degrees.C., tbase=10, tbase_max=60, type = "C")),
                date=x)
    
    Weather90=Data90 %>% group_by(siteplot)%>%
      summarise(precip90= sum(ppt..mm.),
                avgprecip90=mean(ppt..mm.),
                avgtemp90=mean(tmean..degrees.C.),
                maxtemp90=max(tmax..degrees.C.),
                gdd90=max(gdd(tmax=tmax..degrees.C., tmin=tmin..degrees.C., tbase=10, tbase_max=90, type = "C")),
                date=x)
    
    Weather120=Data120 %>% group_by(siteplot)%>%
      summarise(precip120= sum(ppt..mm.),
                avgprecip120=mean(ppt..mm.),
                avgtemp120=mean(tmean..degrees.C.),
                maxtemp120=max(tmax..degrees.C.),
                gdd120=max(gdd(tmax=tmax..degrees.C., tmin=tmin..degrees.C., tbase=10, tbase_max=120, type = "C")),
                date=x)
    
    Weather150=Data150 %>% group_by(siteplot)%>%
      summarise(precip150= sum(ppt..mm.),
                avgprecip150=mean(ppt..mm.),
                avgtemp150=mean(tmean..degrees.C.),
                maxtemp150=max(tmax..degrees.C.),
                gdd150=max(gdd(tmax=tmax..degrees.C., tmin=tmin..degrees.C., tbase=10, tbase_max=150, type = "C")),
                date=x)
    
    Weather180=Data180 %>% group_by(siteplot)%>%
      summarise(precip180= sum(ppt..mm.),
                avgprecip180=mean(ppt..mm.),
                avgtemp180=mean(tmean..degrees.C.),
                maxtemp180=max(tmax..degrees.C.),
                gdd180=max(gdd(tmax=tmax..degrees.C., tmin=tmin..degrees.C., tbase=10, tbase_max=180, type = "C")),
                date=x)
    
    Weather365=Data365 %>% group_by(siteplot)%>%
      summarise(precip365= sum(ppt..mm.),
                avgprecip365=mean(ppt..mm.),
                avgtemp365=mean(tmean..degrees.C.),
                maxtemp365=max(tmax..degrees.C.),
                gdd365=max(gdd(tmax=tmax..degrees.C., tmin=tmin..degrees.C., tbase=10, tbase_max=365, type = "C")),
                date=x)
  
    Weather14=as.data.frame(Weather14)
    Weather30=as.data.frame(Weather30)
    Weather60=as.data.frame(Weather60)
    Weather90=as.data.frame(Weather90)
    Weather120=as.data.frame(Weather120)
    Weather150=as.data.frame(Weather150)
    Weather180=as.data.frame(Weather180)
    Weather365=as.data.frame(Weather365)
    
    Merge1=merge(Weather14, Weather30, by= c("siteplot","date"))
    Merge2=merge(Merge1, Weather60, by= c("siteplot","date"))
    Merge3=merge(Merge2, Weather90, by= c("siteplot","date"))
    Merge4=merge(Merge3, Weather120, by= c("siteplot","date"))
    Merge5=merge(Merge4, Weather150, by= c("siteplot","date"))
    Merge6=merge(Merge5, Weather180, by= c("siteplot","date"))
    Merge7=merge(Merge6, Weather365, by= c("siteplot","date"))
    Merge7
      
}

Data=lapply(levels(raw_data$Date), func)

Fulldata=(do.call(rbind, Data))

#setwd("D:/Dropbox/Work/UT Austin - postdoc/research projects/SWG/analysis/weather_data")

site_info= read.csv("Sampledate_SWG_2Jan2022.csv")
site_info$siteplot=trimws(site_info$siteplot)
site_info$siteplot=as.factor(site_info$siteplot)

Site_weather=left_join(site_info, Fulldata, by=c("siteplot","date"))
Site_weather$date=as.factor(Site_weather$date)

write.csv(Site_weather, file="weather_SWG_samplingdate_2Jan2022.csv")
