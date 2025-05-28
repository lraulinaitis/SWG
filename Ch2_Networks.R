
# LIVIA RAULINAITIS
# CHAPTER 2: SWG POLLINATOR NETWORKS

###############################################################################
#-------------------PACKAGES, FUNCTIONS, & GLOBAL VARIABLES---------------- ----
###############################################################################
# a. packages                                                               ----
library(tidyverse)
library(bipartite)

# b. functions                                                              ----

# c. data sets                                                              ----
  # Load in: [active_plots], [bees]

netbees <- bees |> 
  filter(Method == "netting",
         Genus != is.na(c(Genus)),
         Species != is.na(c(Species)),
         USDA_code != is.na(c(USDA_code)))


#                                                                           ----
###############################################################################
#--------------------------------PRELIMINARIES----------------------------- ----
###############################################################################
# a. prep output dataframe ----

# create dataframe to hold output values
plot_networks <- active_plots

# create new H2 column
plot_networks$H2 <- 0

# b. GRAVEYARD: check for NAs ----
#Check whether any rows have a USDA code but no full plant name
insects_chk1 <- bees |>
  filter(!is.na(full_plant_sp), 
         is.na(USDA_code)) 
unique(insects_chk1$full_plant_sp)
rm(insects_chk1)

#Check whether any rows have a USDA code but no full plant name
insects_chk2 <- bees |> 
  filter(!is.na(USDA_code),
         is.na(full_plant_sp))
unique(insects_chk2$full_plant_sp)
rm(insects_chk2)

# c. GRAVEYARD: turn data frame into matrix ----

matrix <- bees |>
  filter(year == "2022") |>
  filter(Genus_species != "NA_NA") |>
  #filter(site.plot == "Hagerman Bennett Hill") |> # control
  #filter(site.plot == "Hagerman Meadow Pond") |> # seed
  #filter(site.plot == "Hagerman Crow Hill") |> # 50
  #filter(site.plot == "Hagerman Upper") |> # 5
  #filter(Genus_species != "Apis_mellifera") |>
  group_by(USDA_code, Genus) |>
  summarise(count = n_distinct(SWG_ID)) |>
  pivot_wider(names_from = "Genus",
              values_from = "count",
              values_fill = 0) |>
  column_to_rownames("USDA_code")

# d. plot (visualize) webs ----
plotweb(matrix, labsize = .75, text.rot = 90, ybig = 1.5, col.high = "blue", col.low="green")
# arguments:
  # empty = TRUE : should empty rows be ommitted?
  # labsize = 1 : label size, default = 1
  # ybig = 1  : vert. distance b/t upper and lower
  # y.width.low = 0.1 / y.width.high = 0.1 : width of lower/upper boxes
  # low.spacing = NULL / high.spacing = NULL : distance between level boxes
  # col.interaction="grey80" : color of interaction lines
  # col.high = "grey10" / col.low="grey10" : color of upper/lower boxes 
  # high.lablength = NULL / low.lablength = NULL : max # of label characters, NULL = ALL
  # text.rot=0 = text rotation
  # text.high.col="black" / text.low.col="black": text color
  # adj.high=NULL, adj.low=NULL : see ?text

networklevel(matrix_H)
#                                                                          ----
###############################################################################
#----------------------------NETWORK ANALYSIS----------------------------- ----
###############################################################################
''' Create a dataframe with all site years. Use as counter for for loop. 
check which is upper versus lower in adjacency matrix!
'''

# Function to create interaction matrix, calculate H2, and return it
calculate_H2 <- function(filtered_data) {
  # Create an interaction matrix: rows = species, columns = other interacting species or features
  #interaction_matrix <- as.matrix(xtabs(USDA_code ~ Genus_species + year, data = filtered_data))
  interaction_matrix <- bees |>
    group_by(USDA_code, Genus_species) |>
    summarise(count = n_distinct(SWG_ID)) |>
    pivot_wider(names_from = "Genus_species",
                values_from = "count",
                values_fill = 0) |>
    column_to_rownames("USDA_code")
  
  # Calculate specialization (H2) using the bipartite package
  if (nrow(interaction_matrix) > 1 && ncol(interaction_matrix) > 1) {
    H2 <- networklevel(interaction_matrix, index = "H2")
  } else {
    H2 <- NA  # Return NA if matrix isn't sufficient for calculation
  }
  return(H2)
}

# a. calculate plot-level H2 (by plot-year combo)                          ----

# Iterate through each combination of site, treatment, and year
for (i in 1:nrow(plot_networks)) {
  
  # Extract the current site and year combination
  site.plot.counter <- as.character(plot_networks$site.plot[i])
  year.counter <- as.numeric(plot_networks$year[i])
  
  # Filter interaction dataset by the current combination
  filtered_data <- bees %>%
    filter(site.plot == site.plot.counter, year == year.counter) |>
    #filter(Method == "netting") |>
    #filter(Genus_species != "NA_NA") |>
    group_by(USDA_code, Genus) |>
    #group_by(USDA_code, Genus_species) |>
    summarise(count = n_distinct(SWG_ID)) |>
    pivot_wider(names_from = "Genus",
    # pivot_wider(names_from = "Genus_species",)
                values_from = "count",
                values_fill = 0) |>
    column_to_rownames("USDA_code")
  
  # If filtered data is not empty, calculate H2
  if (nrow(filtered_data) > 0) {
    H2 <- networklevel(filtered_data, index = "H2")
  } else {
    H2 <- NA  # No interactions found for this combination
  }
  
  # Store H2 value in the main dataset
  plot_networks$H2[i] <- H2
}


# b. run lower-level d' & add to matrix ----

##### GRAVEYARD##########

networklevel(matrix)
  # arguments: 
    # index = 
null.t.test(matrix) # explore this
plotweb(nullmodel(matrix, N = 1000)) # explore this

