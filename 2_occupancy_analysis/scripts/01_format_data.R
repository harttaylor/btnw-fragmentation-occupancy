# =============================================================================
# 01  FORMAT DATA AND COVARIATES
# =============================================================================
# Format and standardize data for occupancy analysis
#
# INPUTS  : data/2yearvisitmatrix.csv          detection history (508 sites)
#           data/visitsfordetcovs.csv          per-visit detection covariates
#           data/hyp{1,2,3}_landscape_metrics_2m_res.csv   metrics for Def A/B/C
# OUTPUTS : data/def{A,B,C}{150,500,1000}_glmmdata.csv
#
# NOTE: files are named "glmmdata" for historical reasons only; this is just a
# file-format label, not the model type. Script 03 reads them for the Bayesian
# occupancy models.
# =============================================================================

library(tidyr)
library(dplyr)
select <- dplyr::select   # avoid clashes with other packages' select()

# =============================================================================
# LOAD RAW DATA
# =============================================================================
visit_matrix          <- read.csv("data/2yearvisitmatrix.csv")
visits_for_det_covs   <- read.csv("data/visitsfordetcovs.csv")

# Landscape metrics: one file per fragmentation definition (hyp1/2/3 = A/B/C),
# each containing all three buffer sizes in a "buffer_size" column.
habitatA <- read.csv("data/hyp1_landscape_metrics_2m_res.csv")
habitatB <- read.csv("data/hyp2_landscape_metrics_2m_res.csv")
habitatC <- read.csv("data/hyp3_landscape_metrics_2m_res.csv")


# Create separate dataframes for each buffer size
habitatA_150 <- habitatA %>%
  filter(buffer_size == 150) %>%
  select(-buffer_size)

habitatA_500 <- habitatA %>%
  filter(buffer_size == 500) %>%
  select(-buffer_size)

habitatA_1000 <- habitatA %>%
  filter(buffer_size == 1000) %>%
  select(-buffer_size)

habitatB_150 <- habitatB %>%
  filter(buffer_size == 150) %>%
  select(-buffer_size)

habitatB_500 <- habitatB %>%
  filter(buffer_size == 500) %>%
  select(-buffer_size)

habitatB_1000 <- habitatB %>%
  filter(buffer_size == 1000) %>%
  select(-buffer_size)

habitatC_150 <- habitatC %>%
  filter(buffer_size == 150) %>%
  select(-buffer_size)

habitatC_500 <- habitatC %>%
  filter(buffer_size == 500) %>%
  select(-buffer_size)

habitatC_1000 <- habitatC %>%
  filter(buffer_size == 1000) %>%
  select(-buffer_size)

# =============================================================================
# DETECTION HISTORY -> LONG FORMAT, THEN JOIN DETECTION COVARIATES
# =============================================================================
# The visit matrix is wide (columns like X1_1 = year 1 visit 1). We reshape to
# one row per site-year-visit, then attach the per-visit detection covariates.
visit_matrix_long <- visit_matrix %>%
  pivot_longer(cols = starts_with("X"),  # Select columns starting with "X"
               names_to = "visit", 
               values_to = "occupancy") %>%
  mutate(
    year_id = as.numeric(substr(visit, 2, 2)),  # Extract year (1 or 2) from column name
    visit_num = as.numeric(substr(visit, 4, 4))  # Extract visit number (1, 2, or 3)
  ) %>%
  select(gisid, occupancy, year_id, visit_num)  # Keep necessary columns

# Step 2: Merge with detection covariates
visit_with_covs <- visit_matrix_long %>%
  inner_join(
    visits_for_det_covs %>%
      select(gisid, offset, year, sensor, surveyid, year_id, visit_num, duration, Easting, Northing, date_time), 
    by = c("gisid", "year_id", "visit_num")
  )
head(visit_with_covs)
str(visit_with_covs)

# Step 3: Duplicate habitat metrics for each visit within the corresponding year 
# so we can create visit specific habitat metrics for the long dataframe that icnludes entire det history
# First, create a mapping between surveyid and gisid using the visit dataframe
site_mapping <- visit_with_covs %>%
  select(surveyid, gisid, year) %>%
  distinct()


#=========================================================================================
# Defintion A
#=========================================================================================

# 150 m scale -----------------------------------------------------------------------------
# Join this mapping to the habitat metrics to add gisid
habitat_with_sites <- habitatA_150 %>%
  left_join(site_mapping, by = c("surveyid", "survey_year" = "year"))

# Now create a version of habitat metrics without the surveyid
# This prevents duplicate surveyid values in the final join
habitat_for_join <- habitat_with_sites %>%
  select(-surveyid) %>%
  distinct()

# Finally, join the visit data with habitat metrics using gisid and year
defA_150 <- visit_with_covs %>%
  left_join(habitat_for_join, by = c("gisid", "year")) #%>% 

# Check final dataset structure
print(head(defA_150))
print(dim(defA_150))
str(defA_150)
hist(defA_150$prop_habitat)
range(defA_150$prop_habitat)
hist(defA_150$edge_density)
range(defA_150$edge_density)
hist(defA_150$n_patches)
hist(defA_150$mean_patch_area)


# Standardize the dataset
defA_150_std <- defA_150 %>%
  # First, standardize the continuous habitat variables
  mutate(
    # Habitat variables
    prop_habitat_std = scale(prop_habitat, center = TRUE, scale = TRUE)[,1],
    edge_density_std = scale(edge_density, center = TRUE, scale = TRUE)[,1],
    dist_to_habitat_std = scale(dist_to_habitat, center = TRUE, scale = TRUE)[,1],
    in_habitat = ifelse(dist_to_habitat == 0, 1, 0),
    n_patches_std = scale(n_patches, center = TRUE, scale = TRUE)[,1],
    mean_patch_area_std = scale(mean_patch_area, center = TRUE, scale = TRUE)[,1],
    # shape_mn_std = scale(shape_mn, center = TRUE, scale = TRUE)[,1],
    # frac_mn_std = scale(frac_mn, center = TRUE, scale = TRUE)[,1],
    aggregation_std = scale(aggregation_index, center = TRUE, scale = TRUE)[,1],
    # para_cv_std = scale(para_cv, center = TRUE, scale = TRUE)[,1],
    # para_mn_std = scale(para_mn, center = TRUE, scale = TRUE)[,1],
    contagion_std = scale(contagion, center = TRUE, scale = TRUE)[,1],
    
    # Log-transform and then standardize certain metrics if needed
    log_mean_patch_area_std = scale(log(mean_patch_area + 1))[,1],
    log_edge_density_std = scale(log(edge_density + 1))[,1],
    
    # Detection covariates
    #julian_date_std = scale(julian_date, center = TRUE, scale = TRUE)[,1],
    #time_of_day_std = scale(time_of_day, center = TRUE, scale = TRUE)[,1],
    qpad = offset, 
    sensor = factor(sensor),
    year = factor(year),
    duration = factor(duration)
  )

hist(defA_150_std$prop_habitat_std)
hist(defA_150_std$edge_density_std)
hist(defA_150_std$para_cv_std)
hist(defA_150_std$para_mn_std)
hist(defA_150_std$n_patches_std)
hist(defA_150_std$aggregation_std)


# 500 m scale ----------------------------------------------------------------------------
# Join this mapping to the habitat metrics to add gisid
habitat_with_sites <- habitatA_500 %>%
  left_join(site_mapping, by = c("surveyid", "survey_year" = "year"))

# Now create a version of habitat metrics without the surveyid
# This prevents duplicate surveyid values in the final join
habitat_for_join <- habitat_with_sites %>%
  select(-surveyid) %>%
  distinct()

# Finally, join the visit data with habitat metrics using gisid and year
defA_500 <- visit_with_covs %>%
  left_join(habitat_for_join, by = c("gisid", "year")) 

# Check final dataset structure
print(head(defA_500))
print(dim(defA_500))
str(defA_500)
hist(defA_500$prop_habitat)
range(defA_500$prop_habitat)
hist(defA_500$edge_density)
range(defA_500$edge_density)
hist(defA_500$n_patches)
hist(defA_500$mean_patch_area)

# Standardize the dataset
defA_500_std <- defA_500 %>%
  # First, standardize the continuous habitat variables
  mutate(
    # Habitat variables
    prop_habitat_std = scale(prop_habitat, center = TRUE, scale = TRUE)[,1],
    edge_density_std = scale(edge_density, center = TRUE, scale = TRUE)[,1],
    dist_to_habitat_std = scale(dist_to_habitat, center = TRUE, scale = TRUE)[,1],
    in_habitat = ifelse(dist_to_habitat == 0, 1, 0),
    n_patches_std = scale(n_patches, center = TRUE, scale = TRUE)[,1],
    mean_patch_area_std = scale(mean_patch_area, center = TRUE, scale = TRUE)[,1],
    # shape_mn_std = scale(shape_mn, center = TRUE, scale = TRUE)[,1],
    # frac_mn_std = scale(frac_mn, center = TRUE, scale = TRUE)[,1],
    aggregation_std = scale(aggregation_index, center = TRUE, scale = TRUE)[,1],
    # para_cv_std = scale(para_cv, center = TRUE, scale = TRUE)[,1],
    # para_mn_std = scale(para_mn, center = TRUE, scale = TRUE)[,1],
    contagion_std = scale(contagion, center = TRUE, scale = TRUE)[,1],
    
    # Log-transform and then standardize certain metrics if needed
    log_mean_patch_area_std = scale(log(mean_patch_area + 1))[,1],
    log_edge_density_std = scale(log(edge_density + 1))[,1],
    
    # Detection covariates
    #julian_date_std = scale(julian_date, center = TRUE, scale = TRUE)[,1],
    #time_of_day_std = scale(time_of_day, center = TRUE, scale = TRUE)[,1],
    qpad = offset, 
    sensor = factor(sensor),
    year = factor(year),
    duration = factor(duration)
  )

hist(defA_500_std$prop_habitat_std)
hist(defA_500_std$edge_density_std)
hist(defA_500_std$para_cv_std)
hist(defA_500_std$para_mn_std)
hist(defA_500_std$n_patches_std)
hist(defA_500_std$aggregation_std)

# Check distributions 
hist(defA_500$duration)
hist(defA_500$julian_date)
hist(defA_500$time_of_day)
range(defA_500$julian_date) # ranges from 126 to 199 day of year
range(defA_500$time_of_day) # ranges from 3am to 10:30am, mainly between 4:30am and 8:30 am 

# Summary statistics for original variables
summary(defA_500$edge_density)


# 1000 m scale -----------------------------------------------------------------------
# Join this mapping to the habitat metrics to add gisid
habitat_with_sites <- habitatA_1000 %>%
  left_join(site_mapping, by = c("surveyid", "survey_year" = "year"))

# Now create a version of habitat metrics without the surveyid
# This prevents duplicate surveyid values in the final join
habitat_for_join <- habitat_with_sites %>%
  select(-surveyid) %>%
  distinct()

# Finally, join the visit data with habitat metrics using gisid and year
defA_1000 <- visit_with_covs %>%
  left_join(habitat_for_join, by = c("gisid", "year"))

# Check final dataset structure
print(head(defA_1000))
print(dim(defA_1000))
str(defA_1000)
hist(defA_1000$prop_habitat)
range(defA_1000$prop_habitat)
hist(defA_1000$edge_density)
range(defA_1000$edge_density)
hist(defA_1000$n_patches)
hist(defA_1000$mean_patch_area)

# Standardize the dataset
defA_1000_std <- defA_1000 %>%
  # First, standardize the continuous habitat variables
  mutate(
    # Habitat variables
    prop_habitat_std = scale(prop_habitat, center = TRUE, scale = TRUE)[,1],
    edge_density_std = scale(edge_density, center = TRUE, scale = TRUE)[,1],
    dist_to_habitat_std = scale(dist_to_habitat, center = TRUE, scale = TRUE)[,1],
    in_habitat = ifelse(dist_to_habitat == 0, 1, 0),
    n_patches_std = scale(n_patches, center = TRUE, scale = TRUE)[,1],
    mean_patch_area_std = scale(mean_patch_area, center = TRUE, scale = TRUE)[,1],
    # shape_mn_std = scale(shape_mn, center = TRUE, scale = TRUE)[,1],
    # frac_mn_std = scale(frac_mn, center = TRUE, scale = TRUE)[,1],
    aggregation_std = scale(aggregation_index, center = TRUE, scale = TRUE)[,1],
    # para_cv_std = scale(para_cv, center = TRUE, scale = TRUE)[,1],
    # para_mn_std = scale(para_mn, center = TRUE, scale = TRUE)[,1],
    contagion_std = scale(contagion, center = TRUE, scale = TRUE)[,1],
    
    # Log-transform and then standardize certain metrics if needed
    log_mean_patch_area_std = scale(log(mean_patch_area + 1))[,1],
    log_edge_density_std = scale(log(edge_density + 1))[,1],
    
    # Detection covariates
    #julian_date_std = scale(julian_date, center = TRUE, scale = TRUE)[,1],
    #time_of_day_std = scale(time_of_day, center = TRUE, scale = TRUE)[,1],
    qpad = offset, 
    sensor = factor(sensor),
    year = factor(year),
    duration = factor(duration)
  )


hist(defA_1000_std$prop_habitat_std)
hist(defA_1000_std$edge_density_std)
hist(defA_1000_std$para_cv_std)
hist(defA_1000_std$para_mn_std)
hist(defA_1000_std$n_patches_std)
hist(defA_1000_std$aggregation_std)

# Check distributions 
hist(defA_1000$duration)
hist(defA_1000$julian_date)
hist(defA_1000$time_of_day)
range(defA_1000$julian_date) # ranges from 126 to 199 day of year
range(defA_1000$time_of_day) # ranges from 3am to 10:30am, mainly between 4:30am and 8:30 am 

# save all those datasets
write.csv(defA_150_std, "4_occupancy_models/data/defA150_glmmdata.csv")
write.csv(defA_500_std, "4_occupancy_models/data/defA500_glmmdata.csv")
write.csv(defA_1000_std, "4_occupancy_models/data/defA1000_glmmdata.csv")


#=========================================================================================
# Defintion B 
#=========================================================================================

# 150 m scale -----------------------------------------------------------------------------
# Join this mapping to the habitat metrics to add gisid
habitat_with_sites <- habitatB_150 %>%
  left_join(site_mapping, by = c("surveyid", "survey_year" = "year"))

# Now create a version of habitat metrics without the surveyid
# This prevents duplicate surveyid values in the final join
habitat_for_join <- habitat_with_sites %>%
  select(-surveyid) %>%
  distinct()

# Finally, join the visit data with habitat metrics using gisid and year
defB_150 <- visit_with_covs %>%
  left_join(habitat_for_join, by = c("gisid", "year")) #%>% 

# Check final dataset structure
print(head(defB_150))
print(dim(defB_150))
str(defB_150)
hist(defB_150$prop_habitat)
range(defB_150$prop_habitat)
hist(defB_150$edge_density)
range(defB_150$edge_density)
hist(defB_150$n_patches)
hist(defB_150$mean_patch_area)


# Standardize the dataset
defB_150_std <- defB_150 %>%
  # First, standardize the continuous habitat variables
  mutate(
    # Habitat variables
    prop_habitat_std = scale(prop_habitat, center = TRUE, scale = TRUE)[,1],
    edge_density_std = scale(edge_density, center = TRUE, scale = TRUE)[,1],
    dist_to_habitat_std = scale(dist_to_habitat, center = TRUE, scale = TRUE)[,1],
    in_habitat = ifelse(dist_to_habitat == 0, 1, 0),
    n_patches_std = scale(n_patches, center = TRUE, scale = TRUE)[,1],
    mean_patch_area_std = scale(mean_patch_area, center = TRUE, scale = TRUE)[,1],
    # shape_mn_std = scale(shape_mn, center = TRUE, scale = TRUE)[,1],
    # frac_mn_std = scale(frac_mn, center = TRUE, scale = TRUE)[,1],
    aggregation_std = scale(aggregation_index, center = TRUE, scale = TRUE)[,1],
    # para_cv_std = scale(para_cv, center = TRUE, scale = TRUE)[,1],
    # para_mn_std = scale(para_mn, center = TRUE, scale = TRUE)[,1],
    contagion_std = scale(contagion, center = TRUE, scale = TRUE)[,1],
    
    # Log-transform and then standardize certain metrics if needed
    log_mean_patch_area_std = scale(log(mean_patch_area + 1))[,1],
    log_edge_density_std = scale(log(edge_density + 1))[,1],
    
    # Detection covariates
    #julian_date_std = scale(julian_date, center = TRUE, scale = TRUE)[,1],
    #time_of_day_std = scale(time_of_day, center = TRUE, scale = TRUE)[,1],
    qpad = offset, 
    sensor = factor(sensor),
    year = factor(year),
    duration = factor(duration)
  )

hist(defB_150_std$prop_habitat_std)
hist(defB_150_std$edge_density_std)
hist(defB_150_std$para_cv_std)
hist(defB_150_std$para_mn_std)
hist(defB_150_std$n_patches_std)
hist(defB_150_std$aggregation_std)


# 500 m scale ----------------------------------------------------------------------------
# Join this mapping to the habitat metrics to add gisid
habitat_with_sites <- habitatB_500 %>%
  left_join(site_mapping, by = c("surveyid", "survey_year" = "year"))

# Now create a version of habitat metrics without the surveyid
# This prevents duplicate surveyid values in the final join
habitat_for_join <- habitat_with_sites %>%
  select(-surveyid) %>%
  distinct()

# Finally, join the visit data with habitat metrics using gisid and year
defB_500 <- visit_with_covs %>%
  left_join(habitat_for_join, by = c("gisid", "year")) 

# Check final dataset structure
print(head(defB_500))
print(dim(defB_500))
str(defB_500)
hist(defB_500$prop_habitat)
range(defB_500$prop_habitat)
hist(defB_500$edge_density)
range(defB_500$edge_density)
hist(defB_500$n_patches)
hist(defB_500$mean_patch_area)

# Standardize the dataset
defB_500_std <- defB_500 %>%
  # First, standardize the continuous habitat variables
  mutate(
    # Habitat variables
    prop_habitat_std = scale(prop_habitat, center = TRUE, scale = TRUE)[,1],
    edge_density_std = scale(edge_density, center = TRUE, scale = TRUE)[,1],
    dist_to_habitat_std = scale(dist_to_habitat, center = TRUE, scale = TRUE)[,1],
    in_habitat = ifelse(dist_to_habitat == 0, 1, 0),
    n_patches_std = scale(n_patches, center = TRUE, scale = TRUE)[,1],
    mean_patch_area_std = scale(mean_patch_area, center = TRUE, scale = TRUE)[,1],
    # shape_mn_std = scale(shape_mn, center = TRUE, scale = TRUE)[,1],
    # frac_mn_std = scale(frac_mn, center = TRUE, scale = TRUE)[,1],
    aggregation_std = scale(aggregation_index, center = TRUE, scale = TRUE)[,1],
    # para_cv_std = scale(para_cv, center = TRUE, scale = TRUE)[,1],
    # para_mn_std = scale(para_mn, center = TRUE, scale = TRUE)[,1],
    contagion_std = scale(contagion, center = TRUE, scale = TRUE)[,1],
    
    # Log-transform and then standardize certain metrics if needed
    log_mean_patch_area_std = scale(log(mean_patch_area + 1))[,1],
    log_edge_density_std = scale(log(edge_density + 1))[,1],
    
    # Detection covariates
    #julian_date_std = scale(julian_date, center = TRUE, scale = TRUE)[,1],
    #time_of_day_std = scale(time_of_day, center = TRUE, scale = TRUE)[,1],
    qpad = offset, 
    sensor = factor(sensor),
    year = factor(year),
    duration = factor(duration)
  )

hist(defB_500_std$prop_habitat_std)
hist(defB_500_std$edge_density_std)
hist(defB_500_std$para_cv_std)
hist(defB_500_std$para_mn_std)
hist(defB_500_std$n_patches_std)
hist(defB_500_std$aggregation_std)

# Check distributions 
hist(defB_500$duration)
hist(defB_500$julian_date)
hist(defB_500$time_of_day)
range(defB_500$julian_date) # ranges from 126 to 199 day of year
range(defB_500$time_of_day) # ranges from 3am to 10:30am, mainly between 4:30am and 8:30 am 

# Summary statistics for original variables
summary(defB_500$edge_density)


# 1000 m scale -----------------------------------------------------------------------
# Join this mapping to the habitat metrics to add gisid
habitat_with_sites <- habitatB_1000 %>%
  left_join(site_mapping, by = c("surveyid", "survey_year" = "year"))

# Now create a version of habitat metrics without the surveyid
# This prevents duplicate surveyid values in the final join
habitat_for_join <- habitat_with_sites %>%
  select(-surveyid) %>%
  distinct()

# Finally, join the visit data with habitat metrics using gisid and year
defB_1000 <- visit_with_covs %>%
  left_join(habitat_for_join, by = c("gisid", "year"))

# Check final dataset structure
print(head(defB_1000))
print(dim(defB_1000))
str(defB_1000)
hist(defB_1000$prop_habitat)
range(defB_1000$prop_habitat)
hist(defB_1000$edge_density)
range(defB_1000$edge_density)
hist(defB_1000$n_patches)
hist(defB_1000$mean_patch_area)

# Standardize the dataset
defB_1000_std <- defB_1000 %>%
  # First, standardize the continuous habitat variables
  mutate(
    # Habitat variables
    prop_habitat_std = scale(prop_habitat, center = TRUE, scale = TRUE)[,1],
    edge_density_std = scale(edge_density, center = TRUE, scale = TRUE)[,1],
    dist_to_habitat_std = scale(dist_to_habitat, center = TRUE, scale = TRUE)[,1],
    in_habitat = ifelse(dist_to_habitat == 0, 1, 0),
    n_patches_std = scale(n_patches, center = TRUE, scale = TRUE)[,1],
    mean_patch_area_std = scale(mean_patch_area, center = TRUE, scale = TRUE)[,1],
    # shape_mn_std = scale(shape_mn, center = TRUE, scale = TRUE)[,1],
    # frac_mn_std = scale(frac_mn, center = TRUE, scale = TRUE)[,1],
    aggregation_std = scale(aggregation_index, center = TRUE, scale = TRUE)[,1],
    # para_cv_std = scale(para_cv, center = TRUE, scale = TRUE)[,1],
    # para_mn_std = scale(para_mn, center = TRUE, scale = TRUE)[,1],
    contagion_std = scale(contagion, center = TRUE, scale = TRUE)[,1],
    
    # Log-transform and then standardize certain metrics if needed
    log_mean_patch_area_std = scale(log(mean_patch_area + 1))[,1],
    log_edge_density_std = scale(log(edge_density + 1))[,1],
    
    # Detection covariates
    #julian_date_std = scale(julian_date, center = TRUE, scale = TRUE)[,1],
    #time_of_day_std = scale(time_of_day, center = TRUE, scale = TRUE)[,1],
    qpad = offset, 
    sensor = factor(sensor),
    year = factor(year),
    duration = factor(duration)
  )

hist(defB_1000_std$prop_habitat_std)
hist(defB_1000_std$edge_density_std)
hist(defB_1000_std$para_cv_std)
hist(defB_1000_std$para_mn_std)
hist(defB_1000_std$n_patches_std)
hist(defB_1000_std$aggregation_std)

# Check distributions 
hist(defB_1000$duration)
hist(defB_1000$julian_date)
hist(defB_1000$time_of_day)
range(defB_1000$julian_date) # ranges from 126 to 199 day of year
range(defB_1000$time_of_day) # ranges from 3am to 10:30am, mainly between 4:30am and 8:30 am 

# save all those datasets
write.csv(defB_150_std, "4_occupancy_models/data/defB150_glmmdata.csv")
write.csv(defB_500_std, "4_occupancy_models/data/defB500_glmmdata.csv")
write.csv(defB_1000_std, "4_occupancy_models/data/defB1000_glmmdata.csv")


#=========================================================================================
# Defintion C 
#=========================================================================================

# 150 m scale -----------------------------------------------------------------------------
# Join this mapping to the habitat metrics to add gisid
habitat_with_sites <- habitatC_150 %>%
  left_join(site_mapping, by = c("surveyid", "survey_year" = "year"))

# Now create a version of habitat metrics without the surveyid
# This prevents duplicate surveyid values in the final join
habitat_for_join <- habitat_with_sites %>%
  select(-surveyid) %>%
  distinct()

# Finally, join the visit data with habitat metrics using gisid and year
defC_150 <- visit_with_covs %>%
  left_join(habitat_for_join, by = c("gisid", "year")) #%>% 

# Check final dataset structure
print(head(defC_150))
print(dim(defC_150))
str(defC_150)
hist(defC_150$prop_habitat)
range(defC_150$prop_habitat)
hist(defC_150$edge_density)
range(defC_150$edge_density)
hist(defC_150$n_patches)
hist(defC_150$mean_patch_area)


# Standardize the dataset
defC_150_std <- defC_150 %>%
  # First, standardize the continuous habitat variables
  mutate(
    # Habitat variables
    prop_habitat_std = scale(prop_habitat, center = TRUE, scale = TRUE)[,1],
    edge_density_std = scale(edge_density, center = TRUE, scale = TRUE)[,1],
    dist_to_habitat_std = scale(dist_to_habitat, center = TRUE, scale = TRUE)[,1],
    in_habitat = ifelse(dist_to_habitat == 0, 1, 0),
    n_patches_std = scale(n_patches, center = TRUE, scale = TRUE)[,1],
    mean_patch_area_std = scale(mean_patch_area, center = TRUE, scale = TRUE)[,1],
    # shape_mn_std = scale(shape_mn, center = TRUE, scale = TRUE)[,1],
    # frac_mn_std = scale(frac_mn, center = TRUE, scale = TRUE)[,1],
    aggregation_std = scale(aggregation_index, center = TRUE, scale = TRUE)[,1],
    # para_cv_std = scale(para_cv, center = TRUE, scale = TRUE)[,1],
    # para_mn_std = scale(para_mn, center = TRUE, scale = TRUE)[,1],
    contagion_std = scale(contagion, center = TRUE, scale = TRUE)[,1],
    
    # Log-transform and then standardize certain metrics if needed
    log_mean_patch_area_std = scale(log(mean_patch_area + 1))[,1],
    log_edge_density_std = scale(log(edge_density + 1))[,1],
    
    # Detection covariates
    #julian_date_std = scale(julian_date, center = TRUE, scale = TRUE)[,1],
    #time_of_day_std = scale(time_of_day, center = TRUE, scale = TRUE)[,1],
    qpad = offset, 
    sensor = factor(sensor),
    year = factor(year),
    duration = factor(duration)
  )

hist(defC_150_std$prop_habitat_std)
hist(defC_150_std$edge_density_std)
hist(defC_150_std$para_cv_std)
hist(defC_150_std$para_mn_std)
hist(defC_150_std$contagion_std)
hist(defC_150_std$aggregation_std)
hist(defC_150_std$mean_patch_area_std)


# 500 m scale ----------------------------------------------------------------------------
# Join this mapping to the habitat metrics to add gisid
habitat_with_sites <- habitatC_500 %>%
  left_join(site_mapping, by = c("surveyid", "survey_year" = "year"))

# Now create a version of habitat metrics without the surveyid
# This prevents duplicate surveyid values in the final join
habitat_for_join <- habitat_with_sites %>%
  select(-surveyid) %>%
  distinct()

# Finally, join the visit data with habitat metrics using gisid and year
defC_500 <- visit_with_covs %>%
  left_join(habitat_for_join, by = c("gisid", "year")) 

# Check final dataset structure
print(head(defC_500))
print(dim(defC_500))
str(defC_500)
hist(defC_500$prop_habitat)
range(defC_500$prop_habitat)
hist(defC_500$edge_density)
range(defC_500$edge_density)
hist(defC_500$n_patches)
hist(defC_500$mean_patch_area)

# Standardize the dataset
defC_500_std <- defC_500 %>%
  # First, standardize the continuous habitat variables
  mutate(
    # Habitat variables
    prop_habitat_std = scale(prop_habitat, center = TRUE, scale = TRUE)[,1],
    edge_density_std = scale(edge_density, center = TRUE, scale = TRUE)[,1],
    dist_to_habitat_std = scale(dist_to_habitat, center = TRUE, scale = TRUE)[,1],
    in_habitat = ifelse(dist_to_habitat == 0, 1, 0),
    n_patches_std = scale(n_patches, center = TRUE, scale = TRUE)[,1],
    mean_patch_area_std = scale(mean_patch_area, center = TRUE, scale = TRUE)[,1],
    # shape_mn_std = scale(shape_mn, center = TRUE, scale = TRUE)[,1],
    # frac_mn_std = scale(frac_mn, center = TRUE, scale = TRUE)[,1],
    aggregation_std = scale(aggregation_index, center = TRUE, scale = TRUE)[,1],
    # para_cv_std = scale(para_cv, center = TRUE, scale = TRUE)[,1],
    # para_mn_std = scale(para_mn, center = TRUE, scale = TRUE)[,1],
    contagion_std = scale(contagion, center = TRUE, scale = TRUE)[,1],
    
    # Log-transform and then standardize certain metrics if needed
    log_mean_patch_area_std = scale(log(mean_patch_area + 1))[,1],
    log_edge_density_std = scale(log(edge_density + 1))[,1],
    
    # Detection covariates
    #julian_date_std = scale(julian_date, center = TRUE, scale = TRUE)[,1],
    #time_of_day_std = scale(time_of_day, center = TRUE, scale = TRUE)[,1],
    qpad = offset, 
    sensor = factor(sensor),
    year = factor(year),
    duration = factor(duration)
  )

hist(defC_500_std$prop_habitat_std)
hist(defC_500_std$edge_density_std)
hist(defC_500_std$para_cv_std)
hist(defC_500_std$para_mn_std)
hist(defC_500_std$n_patches_std)
hist(defC_500_std$aggregation_std)
hist(defC_500_std$contagion_std)
hist(defC_500_std$mean_patch_area_std)

# Check distributions 
hist(defC_500$duration)
hist(defC_500$julian_date)
hist(defC_500$time_of_day)
range(defC_500$julian_date) # ranges from 126 to 199 day of year
range(defC_500$time_of_day) # ranges from 3am to 10:30am, mainly between 4:30am and 8:30 am 

# Summary statistics for original variables
summary(defC_500$edge_density)


# 1000 m scale -----------------------------------------------------------------------
# Join this mapping to the habitat metrics to add gisid
habitat_with_sites <- habitatC_1000 %>%
  left_join(site_mapping, by = c("surveyid", "survey_year" = "year"))

# Now create a version of habitat metrics without the surveyid
# This prevents duplicate surveyid values in the final join
habitat_for_join <- habitat_with_sites %>%
  select(-surveyid) %>%
  distinct()

# Finally, join the visit data with habitat metrics using gisid and year
defC_1000 <- visit_with_covs %>%
  left_join(habitat_for_join, by = c("gisid", "year"))

# Check final dataset structure
print(head(defC_1000))
print(dim(defC_1000))
str(defC_1000)
hist(defC_1000$prop_habitat)
range(defC_1000$prop_habitat)
hist(defC_1000$edge_density)
range(defC_1000$edge_density)
hist(defC_1000$n_patches)
hist(defC_1000$mean_patch_area)

# Standardize the dataset
defC_1000_std <- defC_1000 %>%
  # First, standardize the continuous habitat variables
  mutate(
    # Habitat variables
    prop_habitat_std = scale(prop_habitat, center = TRUE, scale = TRUE)[,1],
    edge_density_std = scale(edge_density, center = TRUE, scale = TRUE)[,1],
    dist_to_habitat_std = scale(dist_to_habitat, center = TRUE, scale = TRUE)[,1],
    in_habitat = ifelse(dist_to_habitat == 0, 1, 0),
    n_patches_std = scale(n_patches, center = TRUE, scale = TRUE)[,1],
    mean_patch_area_std = scale(mean_patch_area, center = TRUE, scale = TRUE)[,1],
    # shape_mn_std = scale(shape_mn, center = TRUE, scale = TRUE)[,1],
    # frac_mn_std = scale(frac_mn, center = TRUE, scale = TRUE)[,1],
    aggregation_std = scale(aggregation_index, center = TRUE, scale = TRUE)[,1],
    # para_cv_std = scale(para_cv, center = TRUE, scale = TRUE)[,1],
    # para_mn_std = scale(para_mn, center = TRUE, scale = TRUE)[,1],
    contagion_std = scale(contagion, center = TRUE, scale = TRUE)[,1],
    
    # Log-transform and then standardize certain metrics if needed
    log_mean_patch_area_std = scale(log(mean_patch_area + 1))[,1],
    log_edge_density_std = scale(log(edge_density + 1))[,1],
    
    # Detection covariates
    #julian_date_std = scale(julian_date, center = TRUE, scale = TRUE)[,1],
    #time_of_day_std = scale(time_of_day, center = TRUE, scale = TRUE)[,1],
    qpad = offset, 
    sensor = factor(sensor),
    year = factor(year),
    duration = factor(duration)
  )

hist(defC_1000_std$prop_habitat_std)
hist(defC_1000_std$edge_density_std)
hist(defC_1000_std$para_cv_std)
hist(defC_1000_std$para_mn_std)
hist(defC_1000_std$n_patches_std)
hist(defC_1000_std$aggregation_std)
hist(defC_1000_std$contagion_std)
hist(defC_1000_std$mean_patch_area_std)

# Check distributions 
hist(defC_1000$duration)
hist(defC_1000$julian_date)
hist(defC_1000$time_of_day)
range(defC_1000$julian_date) # ranges from 126 to 199 day of year
range(defC_1000$time_of_day) # ranges from 3am to 10:30am, mainly between 4:30am and 8:30 am 

# save all those datasets
write.csv(defC_150_std, "4_occupancy_models/data/defC150_glmmdata.csv")
write.csv(defC_500_std, "4_occupancy_models/data/defC500_glmmdata.csv")
write.csv(defC_1000_std, "4_occupancy_models/data/defC1000_glmmdata.csv")


# check ranges, means and sds for all fragmentation metrics 
range(defA_150$mean_patch_area)
range(defA_150$edge_density)
range(defA_150$aggregation_index)
range(defA_150$prop_habitat)

range(defA_500$mean_patch_area)
range(defA_500$edge_density)
range(defA_500$aggregation_index)
range(defA_500$prop_habitat)

range(defA_1000$mean_patch_area)
range(defA_1000$edge_density)
range(defA_1000$aggregation_index)
range(defA_1000$prop_habitat)

range(defB_150$mean_patch_area)
range(defB_150$edge_density)
range(defB_150$aggregation_index)
range(defB_150$prop_habitat)

range(defB_500$mean_patch_area)
range(defB_500$edge_density)
range(defB_500$aggregation_index)
range(defB_500$prop_habitat)

range(defB_1000$mean_patch_area)
range(defB_1000$edge_density)
range(defB_1000$aggregation_index)
range(defB_1000$prop_habitat)

range(defC_150$mean_patch_area)
range(defC_150$edge_density)
range(defC_150$aggregation_index)
range(defC_150$prop_habitat)

range(defC_500$mean_patch_area)
range(defC_500$edge_density)
range(defC_500$aggregation_index)
range(defC_500$prop_habitat)

range(defC_1000$mean_patch_area)
range(defC_1000$edge_density)
range(defC_1000$aggregation_index)
range(defC_1000$prop_habitat)



# =============================================================================
# HELPER: build one standardised modelling dataset
# =============================================================================
# Given the landscape-metric rows for one definition at one buffer size, this
# attaches gisid, joins the metrics onto the visit-level detection data, and
# standardises the continuous covariates (mean 0, SD 1). Standardised columns
# end in "_std". (Script 03 re-standardises from the raw columns as well, so
# these are provided mainly for inspection.)

build_dataset <- function(habitat_scale) {

  # Attach gisid to the survey-keyed metrics, then drop surveyid so the final
  # join does not duplicate rows.
  habitat_for_join <- habitat_scale %>%
    left_join(site_mapping, by = c("surveyid", "survey_year" = "year")) %>%
    select(-surveyid) %>%
    distinct()

  # Join metrics onto every visit at that site-year.
  merged <- visit_with_covs %>%
    left_join(habitat_for_join, by = c("gisid", "year"))

  # Standardise.
  merged %>%
    mutate(
      prop_habitat_std        = scale(prop_habitat)[, 1],
      edge_density_std        = scale(edge_density)[, 1],
      dist_to_habitat_std     = scale(dist_to_habitat)[, 1],
      in_habitat              = ifelse(dist_to_habitat == 0, 1, 0),
      n_patches_std           = scale(n_patches)[, 1],
      mean_patch_area_std     = scale(mean_patch_area)[, 1],
      aggregation_std         = scale(aggregation_index)[, 1],
      contagion_std           = scale(contagion)[, 1],
      log_mean_patch_area_std = scale(log(mean_patch_area + 1))[, 1],
      log_edge_density_std    = scale(log(edge_density + 1))[, 1],
      qpad     = offset,
      sensor   = factor(sensor),
      year     = factor(year),
      duration = factor(duration)
    )
}

# =============================================================================
# BUILD AND SAVE ALL NINE DATASETS
# =============================================================================
# Each habitat file holds all three buffers; filter() picks one buffer, then
# build_dataset() does the join + standardisation, and we write the result.

# ---- Definition A -----------------------------------------------------------
write.csv(build_dataset(filter(habitatA, buffer_size == 150)),
          "data/defA150_glmmdata.csv",  row.names = FALSE)
write.csv(build_dataset(filter(habitatA, buffer_size == 500)),
          "data/defA500_glmmdata.csv",  row.names = FALSE)
write.csv(build_dataset(filter(habitatA, buffer_size == 1000)),
          "data/defA1000_glmmdata.csv", row.names = FALSE)

# ---- Definition B -----------------------------------------------------------
write.csv(build_dataset(filter(habitatB, buffer_size == 150)),
          "data/defB150_glmmdata.csv",  row.names = FALSE)
write.csv(build_dataset(filter(habitatB, buffer_size == 500)),
          "data/defB500_glmmdata.csv",  row.names = FALSE)
write.csv(build_dataset(filter(habitatB, buffer_size == 1000)),
          "data/defB1000_glmmdata.csv", row.names = FALSE)

# ---- Definition C -----------------------------------------------------------
write.csv(build_dataset(filter(habitatC, buffer_size == 150)),
          "data/defC150_glmmdata.csv",  row.names = FALSE)
write.csv(build_dataset(filter(habitatC, buffer_size == 500)),
          "data/defC500_glmmdata.csv",  row.names = FALSE)
write.csv(build_dataset(filter(habitatC, buffer_size == 1000)),
          "data/defC1000_glmmdata.csv", row.names = FALSE)

cat("Wrote 9 modelling datasets to data/.\n")

# ---- (optional) quick sanity checks on one dataset --------------------------
# defC150 <- read.csv("data/defC150_glmmdata.csv")
# range(defC150$prop_habitat); range(defC150$edge_density)
# hist(defC150$prop_habitat);  hist(defC150$edge_density)
