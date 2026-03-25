################################################################################
######## R Workflow to Extract Forest attributes from LiDAR Point Cloud ########
########         Lilian VALLET, Postdoctoral Researcher at CSU          ########
################################################################################
T1 <- Sys.time()
# Packages ---------------------------------------------------------------------
list.of.packages <- c("lidR","readr","dplyr","purrr","furrr","lasR","rlang","tidyverse",
                      "terra", "tidyterra", "devtools","broom","caret","tigris")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
library(devtools)
#devtools::install_github('hunter-stanke/rFIA')
#devtools::install_github("georgewoolsey/cloud2trees")
library(rlang)
library(lasR)
library(lidR)
library(terra)
library(furrr)
library(purrr)
library(cloud2trees)
library(tidyverse)
library(tidyterra)
library(tidyr)
library(dplyr)
library(readr)
library(tigris)
library(rFIA)
library(broom)
library(caret)

# Input and output Folders -----------------------------------------------------
dataP <- "D:/" #Path of input data (Containing one or several LiDAR campaign)
projectP <- "C:/Users/c837666159/Documents/Recherche/CSU/Projet_RGM/" #Project path
project_dataP <- paste0(projectP, "data/") #Auxiliary data path in project path
pixel_size <- 10
core_nb <- 30 #Number of core to use for parallel processing
window.size = function(x) {
  (x * 0.1) + 3
}#Window size custom function for tree segmentation
options(scipen = 999) 

# LiDAR Campaign ---------------------------------------------------------------
roi <- "CO_ArapahoRooseveltPikeNF_D23" #Name of the LiDAR mission
campaign <- "" #Name of the LiDAR campaign. If mission = campaign, keep ""
lidar_folder <- "LAS_Tiles"#Name of the .las folder, Usually "LAS_tiles", "Point_Cloud" or "LiDAR_Point_Cloud"
lidar_file_type <- ".laz"

# Folders and Paths Initialization ---------------------------------------------
roiP <- paste0(dataP, roi, "/") #Set a path to LiDAR mission folder
campaignP <- paste0(roiP, campaign, "/") #Set a path to LiDAR campaign folder
lasP <- paste0(campaignP, lidar_folder, "/") #Set a path to .las folder

outputP <- paste0(dataP) #Set a path to folder which will contains output
output_campaignP <- paste0(outputP, roi, "/", campaign, "/", "output_data/") #Set a path to outputs of the LiDAR campaign
output_campaign_mapP <- paste0(output_campaignP, "Map/") #Set a path to final maps of LiDAR campaign
output_campaign_tilesP <- paste0(output_campaignP, "Tiles/") #Set a path to individual tiles of LiDAR campaign
output_campaign_metadataP <- paste0(output_campaignP, "MetaData/") #Set a path to individual tiles of LiDAR campaign


dir.create(outputP, showWarnings = F) #Create above-mentioned folder
dir.create(output_campaignP, showWarnings = F) #Create above-mentioned folder
dir.create(output_campaign_mapP, showWarnings = F) #Create above-mentioned folder
dir.create(output_campaign_tilesP, showWarnings = F) #Create above-mentioned folder
dir.create(output_campaign_metadataP, showWarnings = F) #Create above-mentioned folder


# File listing and formatting
LasFiles <- list.files(lasP, pattern = paste0(lidar_file_type,"$")) #List .las files
LasFiles <- LasFiles[!grepl(
  paste0("_normalized", lidar_file_type, "$"), # e.g. _normalized\.laz$
  LasFiles
)]
TilesID <- str_sub(LasFiles, start = 1, end = -5) #List Tiles ID
LaxFiles <- list.files(lasP, pattern = ".lax$") #List .lax files
MissingLax <- LasFiles[!TilesID %in% str_sub(LaxFiles, start = 1, end =
                                              -5)] #List missing .lax files

if (length(LasFiles) != length(LaxFiles)) {
  #If any missing .lax files
  print(paste0("Writing ", length(MissingLax), " lax files")) #Print message
  lasR::exec(write_lax(),
             on = paste0(lasP, MissingLax),
             ncores = concurrent_files(30))#Write missing lax files
}

# Write pipeline ---------------------------------------------------------------
class_noise <- classify_with_ivf() #Noise classification
del_noise <- delete_noise() #Noise filtering
triang <- triangulate(filter = keep_ground()) #Ground classification through triangulation
norm_height <- transform_with(triang) #Height normalization
write <- write_las(paste0(lasP, "*_normalized.laz")) #Write normalized point cloud
chm1 <- lasR::rasterize(1, "max") #CHM at 1m
chm1_fill = pit_fill(chm1, ofile = paste0(output_campaign_tilesP, "*_chm1_fill.tif")) #Gap filling
non_ground <- lasR::rasterize(
  pixel_size,
  "count",
  filter = "Classification != 2",
  ofile = paste0(output_campaign_tilesP, "*_nongroundNB.tif")
) #Nb non ground returns
ground <- lasR::rasterize(
  pixel_size,
  "count",
  filter = "Classification == 2",
  ofile = paste0(output_campaign_tilesP, "*_groundNB.tif")
) #Nb ground returns
sup2 <- lasR::rasterize(
  pixel_size,
  "count",
  filter = "Z > 2",
  default_value = 0,
  ofile = paste0(output_campaign_tilesP, "*_sup2NB.tif")
) #Nb returns above 2m
bet12 <- lasR::rasterize(
  pixel_size,
  "count",
  filter = "Z %between% 1 2",
  default_value = 0,
  ofile = paste0(output_campaign_tilesP, "*_bet12NB.tif")
) #Nb returns between 1 and 2m
bet05 <- lasR::rasterize(
  pixel_size,
  "count",
  filter = "Z %between% -1 5",
  default_value = 0,
  ofile = paste0(output_campaign_tilesP, "*_bet05NB.tif")
) #Nb returns between 0 and 5m
bet510 <- lasR::rasterize(
  pixel_size,
  "count",
  filter = "Z %between% 5 10",
  default_value = 0,
  ofile = paste0(output_campaign_tilesP, "*_bet510NB.tif")
) #Nb returns between 5 and 10m
bet1015 <- lasR::rasterize(
  pixel_size,
  "count",
  filter = "Z %between% 10 15",
  default_value = 0,
  ofile = paste0(output_campaign_tilesP, "*_bet1015NB.tif")
) #Nb returns between 10 and 15m
bet1520 <- lasR::rasterize(
  pixel_size,
  "count",
  filter = "Z %between% 15 20",
  default_value = 0,
  ofile = paste0(output_campaign_tilesP, "*_bet1520NB.tif")
) #Nb returns between 15 and 20m
bet2025 <- lasR::rasterize(
  pixel_size,
  "count",
  filter = "Z %between% 20 25",
  default_value = 0,
  ofile = paste0(output_campaign_tilesP, "*_bet2025NB.tif")
) #Nb returns between 20 and 25m
bet2530 <- lasR::rasterize(
  pixel_size,
  "count",
  filter = "Z %between% 25 30",
  default_value = 0,
  ofile = paste0(output_campaign_tilesP, "*_bet2530NB.tif")
) #Nb returns between 25 and 30m
all <- lasR::rasterize(pixel_size, "count", ofile = paste0(output_campaign_tilesP, "*_allNB.tif")) #Nb of all

pipeline <- class_noise + del_noise + triang + norm_height + write + chm1 +
  chm1_fill + non_ground + ground + sup2  + bet05 + bet510 + bet1015 + bet1520 +
  bet2025 + bet2530 + all + bet12

#Run pipeline ------------------------------------------------------------------
TilesDone <- str_sub(list.files(output_campaign_tilesP, "allNB"),
                       end=-11) #List Tiles already processed

if (length(TilesDone) > 1) {
  #If some tiles are already processed ...
  TilesToDo <- TilesID[-c(which(TilesID %in% TilesDone))] #List Tiles to be processed
} else{
  TilesToDo <- TilesID #Otherwise, All Tiles should be processed
}

f <- c(paste0(lasP, TilesToDo, lidar_file_type),
       paste0(lasP, TilesToDo, ".lax")) #Merge las and lax file path

ans = lasR::exec(pipeline, on = f, ncores = concurrent_files(core_nb)) #Execute pipelines

# Obtain Allometry relation ----------------------------------------------------
CAMPAIGN_MAP <- map(list.files(output_campaign_tilesP, "sup2NB.tif$", full.names = T)[100:300],
            rast)#Read all tiles
#If you want to get a gpkg file of the tile grid, uncomment this section :
CAMPAIGN_TILES <- map(CAMPAIGN_MAP,function(data){
  tile_id <- str_sub(basename(sources(data)),end=-12)
  Polygon <- as.polygons(ext(data))
  Polygon$tile_id=tile_id
  return(Polygon)
})
Campaign_Tiles <- vect(svc(CAMPAIGN_TILES))
writeVector(Campaign_Tiles,paste0(output_campaign_metadataP,"Campaign_Tiles_",roi,campaign,".gpkg"),overwrite=T)
Campaign_Map <- mosaic(
  sprc(CAMPAIGN_MAP),
  filename = paste0(output_campaign_mapP, roi, "_", campaign, "_sup2NB.tif"),
  overwrite = T
)#Merge Tiles and write final map
Campaign_Map_LowRes <- aggregate(Campaign_Map,fun="any",fact=50,na.rm=T)
Campaign_Map_Border <- as.polygons(Campaign_Map_LowRes)
Campaign_Map_Buffer <- Campaign_Map_Border%>%
  project("EPSG:4326") %>%
  buffer(2000)

# Determine states
States <- vect(tigris::states(cb = TRUE, year = 2021))
States_In_Aoi <- States[Campaign_Map_Buffer, ]
NeededStates <- unique(States_In_Aoi$STUSPS)

# Acquire FIA data
Fia <- getFIA(
  states = NeededStates,
  load = TRUE,
  tables = c("PLOT", "TREE", "COND", "PLOT"))

Ref <- getFIA(states = "REF", load = TRUE)

Fia_Plot_Locations <- vect(
  Fia$PLOT,
  geom = c("LON", "LAT"),
  crs = "EPSG:4326")

Plots_In_Aoi <- Fia_Plot_Locations[Campaign_Map_Buffer, ]
PlotsDF <- as.data.frame(Plots_In_Aoi)

# Identify ForTypGrp
ForTypeGrp_LUT <- Ref$REF_FOREST_TYPE_GROUP %>% select(VALUE, MEANING)
colnames(ForTypeGrp_LUT) <- c("VALUE", "ForTypGrp")
Species_LUT <- Ref$REF_SPECIES %>% select(SPCD, COMMON_NAME)

Fia_Cond <- Fia$COND %>%
  filter(PLT_CN %in% PlotsDF$CN) %>%
  group_by(PLT_CN) %>%
  slice(which.max(INVYR)) %>%
  ungroup()
Fia_Cond$FORTYPGRP <- floor(Fia_Cond$FORTYPCD / 10) * 10
Fia_ForType <- Fia_Cond %>%
  select(PLT_CN, FORTYPCD, FORTYPGRP)%>%
  left_join(., ForTypeGrp_LUT, by = c("FORTYPGRP" = "VALUE")) %>%
  na.omit()

# Remove forest type groups with < 10 plots
Fia_ForType <- Fia_ForType %>%
  group_by(ForTypGrp) %>%
  filter(n() >= 10) %>%
  ungroup()

# Create TREE dataset
Study_Area_Tree <- Fia$TREE
Study_Area_Tree <- Study_Area_Tree %>%
  group_by(CN) %>%
  slice(which.max(INVYR)) %>%
  ungroup() %>%
  left_join(., Fia_ForType, by = "PLT_CN")%>%
  mutate(DIA = DIA *2.45, ACTUALHT = ACTUALHT *0.3048)%>%
  filter(!is.na(ForTypGrp))

# Calculate species basal area in each plot
Basal_Area_Summary <- Study_Area_Tree %>%
  mutate(BA_ft2 = (pi / 4) * (DIA^2) * TPA_UNADJ / 144) %>%
  group_by(PLT_CN, SPCD) %>%
  summarise(
    basal_area_acre = sum(BA_ft2, na.rm = TRUE),
    .groups = "drop")%>%
  left_join(., Species_LUT, by = "SPCD")

# Calculate the dominant species by basal area in each plot
Dominant_Species_Plot <- Basal_Area_Summary %>%
  group_by(PLT_CN) %>%
  slice_max(basal_area_acre, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  left_join(., Fia_ForType, by = "PLT_CN")

# Calculate the most common dominant species in each ForTypGrp across plots 
Dominant_Species_ForTypGrp <- Dominant_Species_Plot %>%
  group_by(ForTypGrp, COMMON_NAME) %>%
  summarize(count = n(), .groups = "drop_last") %>%
  group_by(ForTypGrp) %>%
  mutate(total_plots = sum(count)) %>%
  slice_max(count, n = 1, with_ties = FALSE) %>%
  ungroup()

Species_Counts <- Dominant_Species_Plot %>%
  group_by(ForTypGrp, COMMON_NAME) %>%
  summarize(count = n())

Total_Counts <- Species_Counts %>%
  group_by(ForTypGrp) %>%
  summarize(total_plots = sum(count))

# Join data
Dominant_Species_ForTypGrp <- Species_Counts %>%
  group_by(ForTypGrp) %>%
  slice_max(count, n = 1, with_ties = FALSE) %>%
  left_join(Total_Counts, by = "ForTypGrp")
colnames(Dominant_Species_ForTypGrp) <- c("ForTypGrp", "dominant_species", "n_plots_dominant", "total_plots")

Study_Area_Tree_Measured <- Study_Area_Tree%>%
  # remove dead trees, predicted Ht, predicted DBH, and suppressed crow classes 
  filter(STATUSCD == 1 & DIAHTCD == 1 & HTCD == 1 & CCLCD %in% c(1,2,3,4)) %>%
  filter(!is.na(ACTUALHT)) 
Study_Area_Tree_Measured <- left_join(Study_Area_Tree_Measured, Dominant_Species_ForTypGrp, by = "ForTypGrp")

#defining the equation as a function
expoGrowth.fun <- function(ht,init,k){
  init*((ht+1.37)^k)}

#global model with all species
Nls_All <- nls(DIA ~ expoGrowth.fun(ACTUALHT, init, k), 
               data = Study_Area_Tree_Measured,
               start = list(init = 1, k = 1.0))

#global model confidence intervals of model parameters
Nls_All_Coefs <- tidy(Nls_All, conf.int = TRUE) %>%
  mutate(ForTypGrp = "All species") %>%
  select(ForTypGrp, term, estimate, std.error, conf.low, conf.high, statistic, p.value)

#get parameters coefficients and confidence intervals for each species
Nls_By_ForTypGrp_Coefs <- Study_Area_Tree_Measured %>% 
  group_by(ForTypGrp) %>%
  do(tidy(nls(DIA ~ expoGrowth.fun(ACTUALHT, init, k),
              data = .,
              start = list(init = 1, k = 1),
              algorithm = "port"), conf.int = TRUE))%>%
  select(ForTypGrp, term, estimate, std.error, conf.low, conf.high, statistic, p.value)

# merge
Nls_Coefs <- rbind(Nls_All_Coefs, Nls_By_ForTypGrp_Coefs) %>%
  mutate_at(vars(ForTypGrp), factor)


# Remove Forest types with model parameters with p-values > 0.1
Valid_Models <- Nls_Coefs %>%
  group_by(ForTypGrp) %>%
  filter(all(p.value <= 0.1)) %>%
  ungroup()

# Restructure 
Nls_Coefs_Wide <- Valid_Models %>%
  pivot_wider(
    id_cols = ForTypGrp,
    names_from = term,
    values_from = c(estimate, std.error, conf.low, conf.high, statistic, p.value),
    names_glue = "{.value}_{term}")

# Join with dominant species
Nls_Coefs_DominantBA <- left_join(Nls_Coefs_Wide, Dominant_Species_ForTypGrp, by = "ForTypGrp")

# Write Coefficient as a file
write.csv(Nls_Coefs_DominantBA,paste0(output_campaign_metadataP,"Nls_Coefs_",roi,campaign,".csv"))

# Declare Function for forest attribute extraction -----------------------------
ForestAttributes <- c(
  "chm",
  "CC",
  "SubCC",
  "VC",
  "LadderFuelC",
  "foresttype",
  "dbhmean",
  "dbhmedian",
  "dbhsd",
  "QMD",
  "treedensity",
  "V_m3_sum",
  "AGB_kg_sum",
  "Wwood_kg_sum",
  "Wfoliage_kg_sum",
  "Carbon_kg_sum"
)
process_cbh=F
if(process_cbh==T){
  ForestAttributes <- c(
    ForestAttributes,
    "cbhlowestmean",
    "cbhlowestmedian",
    "cbhlowestmin",
    "cbhlowest20th",
    "cbhlowestsd"
  )
}

process.tile <- function(tile_id,process_cbh=F,forest_mask=T) {
  Chm1 <- rast(paste0(output_campaign_tilesP, tile_id, "_chm1_fill.tif")) #Read CHM at 1m
  Chm1 <- clamp(Chm1, lower = 0) #Assure that CHM value are not below 0
  #Process tile only if tile contains more than 5000 pixels (mandatory for tree segmentation)
  #CHM --------------------------------------------------------------------
  if(nrow(Chm1)>= pixel_size&ncol(Chm1) >= pixel_size){
  Chm <- terra::aggregate(
    Chm1,
    pixel_size,
    fun = "max",
    na.rm = T,
    filename = paste0(output_campaign_tilesP, tile_id, "_chm_",pixel_size,".tif"),
    overwrite = T
  ) #Rasterize CHM 1m to CHM
  
  #Canopy cover, Sub canopy cover and Ladder fuel cover ----------------------
  #Read required files
  allNB <- rast(paste0(output_campaign_tilesP, tile_id, "_allNB.tif"))
  sup2NB <- rast(paste0(output_campaign_tilesP, tile_id, "_sup2NB.tif"))
  bet05 <- rast(paste0(output_campaign_tilesP, tile_id, "_bet05NB.tif"))
  bet12 <- rast(paste0(output_campaign_tilesP, tile_id, "_bet12NB.tif"))
  bet25 <- (bet05 + sup2NB) - allNB
  
  VC <- (bet12 + sup2NB) / allNB#Calculate Vegetation cover
  writeRaster(VC,
              paste0(output_campaign_tilesP, tile_id, "_VC_",pixel_size,".tif"),
              overwrite = T)#Write vegetation cover
  CC <- (sup2NB - bet25) / allNB#Calculate Canopy cover
  writeRaster(CC,
              paste0(output_campaign_tilesP, tile_id, "_CC_",pixel_size,".tif"),
              overwrite = T)#Write Canopy cover
  SubCC <- bet25 / bet05#Calculate sub canopy cover
  writeRaster(SubCC,
              paste0(output_campaign_tilesP, tile_id, "_SubCC_",pixel_size,".tif"),
              overwrite = T)#Write sub canopy cover
  LadderFuelC <- bet12 / (bet05 - bet25)#Calculate ladder fuel cover
  writeRaster(
    LadderFuelC,
    paste0(output_campaign_tilesP, tile_id, "_LadderFuelC_",pixel_size,".tif"),
    overwrite = T
  )#Write ladder fuel cover
  if (global(Chm1, quantile, prob = 0.9999, na.rm = T)$X99. >= 2) {
    #Process Tile only if not only low vegetation
    #Crown polygons ----------------------------------------------------------
    CrownsSF <- cloud2trees::raster2trees(chm_rast = Chm1,
                                          ws = window.size,
                                          outfolder = tempdir())#Tree segmentation
    writeVector(
      vect(CrownsSF),
      paste0(output_campaign_tilesP, tile_id, "_crowns.gpkg"),
      overwrite = T
    )#Write tree segmentation output
    if (nrow(CrownsSF) > 5) {
      #Process tile only if at least 5 trees
      #Forest type -----------------------------------------------------------
      Trees <- cloud2trees::trees_type(tree_list = CrownsSF, max_search_dist_m = 88)#Obtain forest type
      writeRaster(
        Trees$foresttype_rast,
        paste0(output_campaign_tilesP, tile_id, "_foresttype.tif"),
        overwrite = T
      )#Write forest type raster
      
      #Tree density ----------------------------------------------------------
      Tree_Count <- terra::rasterize(Trees$tree_list,
                                     Chm,
                                     fun = "count",
                                     background = 0) #Count number of tree per pixel
      Cell_Area <- cellSize(Tree_Count) #Calculate area of each pixel
      Tree_Density <- Tree_Count / Cell_Area * 10000 #Calculate tree density as nb of tree per hectare
      writeRaster(
        Tree_Density,
        paste0(output_campaign_tilesP, tile_id, "_treedensity_",pixel_size,".tif"),
        overwrite = T
      )#Write tree density raster
      
      
      #DBH -------------------------------------------------------------------
      Trees_Dbh <- Trees$tree_list #Extract tree inventory from segmented trees
      Trees_Dbh <- left_join(Trees_Dbh,
                             Nls_Coefs_DominantBA,
                             by = c("forest_type_group" = "ForTypGrp")) #Join inventory with parameters
      Trees_Dbh %>%
        mutate(
          init = if_else(is.na(estimate_init), Nls_Coefs_DominantBA$estimate_init[1], estimate_init),
          k = if_else(is.na(estimate_k), Nls_Coefs_DominantBA$estimate_k[1], estimate_k)
        ) -> Trees_Dbh #If the specis-specific parameter is not available, get the global parameters
      Trees_Dbh %>%
        mutate(tree_dbh_cm = init * ((tree_height_m + 1.37)^k)) -> Trees_Dbh #Apply DBH model
      Trees_Dbh %>%
        mutate(tree_dbh_cm = if_else(is.na(tree_dbh_cm), 500, tree_dbh_cm)) ->
        Trees_Dbh #Set maximum value of DBH at 500cm
      
      Dbh_mean <- terra::rasterize(
        Trees_Dbh,
        Chm,
        "tree_dbh_cm",
        fun = "mean",
        na.rm = T,
        filename = paste0(output_campaign_tilesP, tile_id, "_dbhmean_",pixel_size,".tif"),
        overwrite = T
      ) #rasterize mean value of DBH
      Dbh_median <- terra::rasterize(
        Trees_Dbh,
        Chm,
        "tree_dbh_cm",
        fun = median,
        na.rm = T,
        filename = paste0(output_campaign_tilesP, tile_id, "_dbhmedian_",pixel_size,".tif"),
        overwrite = T
      ) #rasterize median value of DBH
      Dbh_sd <- terra::rasterize(
        Trees_Dbh,
        Chm,
        "tree_dbh_cm",
        fun = sd,
        na.rm = T,
        filename = paste0(output_campaign_tilesP, tile_id, "_dbhsd_",pixel_size,".tif"),
        overwrite = T
      ) #rasterize standard deviation value of DBH
      
      #QMD -------------------------------------------------------------------
      Sum_Square_DBH <- terra::rasterize(
        Trees_Dbh,
        Chm,
        field = "tree_dbh_cm",
        fun = function(x) {
          sum(x^2)
        },
        background = 0
      )
      QMD = sqrt(Sum_Square_DBH / Tree_Count)
      writeRaster(QMD,
                  paste0(output_campaign_tilesP, tile_id, "_QMD_",pixel_size,".tif"),
                  overwrite = T)#Write tree density raster
      
      
      #Biomass ---------------------------------------------------------------
      Trees_Biomass <- Trees_Dbh
      Ecoregions <- vect(paste0(project_dataP, "Biomass_data/Eco_Division_t.gpkg")) #Read Ecoregion file
      if (crs(Ecoregions) != crs(Trees_Biomass)) {
        Ecoregions <- terra::project(Ecoregions, crs(Trees_Biomass))#reproject Ecoregion if needed
      }
      Trees_Biomass <- vect(Trees_Biomass)
      Trees_Biomass <- terra::intersect(Trees_Biomass, Ecoregions)#Join ecoregion and tree inventory
      Trees_Biomass <- sf::st_as_sf(Trees_Biomass)#Change object type for NSVB processing
      Trees_Biomass %>%
        mutate(
          DIVISION = PROVINCE,
          PROVINCE = NA,
          STDORGCD = NA,
          DIA = tree_dbh_cm,
          HT = tree_height_m,
          ACTUALHT = tree_height_m,
          CR = NA,
          CULL = 0,
          DECAYCD = NA,
          STATUSCD = 1
        ) -> Trees_Biomass #Adapt naming for NSVB
      if (any(is.na(Trees_Biomass$SPCD))) {
        #Is there any tree without a species associated ?
        Trees_Biomass %>%
          filter(!is.na(SPCD)) -> Trees_Biomass#Filter it
      }
      if (nrow(Trees_Biomass) > 0) {
        #Process Biomass only if at least 1 tree is present
        source(paste0(projectP, "script/nsvb_function.R"))#Read NSVB function
        Trees <- nsvb(Trees_Biomass)#Apply NSVB
        Trees %>%
          mutate(
            V_tot_ob_Gross_m3 = V_tot_ob_Gross / 35.31469989,
            AGB_pred_kg = AGB_pred / 2.20462262185,
            W_Wood_kg = (Wood_harm + Bark_harm) / 2.20462262185,
            W_foliage_kg = W_foliage / 2.20462262185,
            C_kg = C / 2.20462262185
          ) -> Trees #Adapt units
        
        writeVector(
          vect(Trees),
          paste0(output_campaign_tilesP, tile_id, "_trees.gpkg"),
          overwrite = T
        )#Write final Tree inventory
        V_sum <- terra::rasterize(
          Trees,
          Chm,
          "V_tot_ob_Gross_m3",
          fun = sum,
          na.rm = T,
          filename = paste0(output_campaign_tilesP, tile_id, "_V_m3_sum_",pixel_size,".tif"),
          overwrite = T
        )#Write volume raster
        AGB_sum <- terra::rasterize(
          Trees,
          Chm,
          "AGB_pred_kg",
          fun = sum,
          na.rm = T,
          filename = paste0(output_campaign_tilesP, tile_id, "_AGB_kg_sum_",pixel_size,".tif"),
          overwrite = T
        )#Write biomass raster
        W_Wood_sum <- terra::rasterize(
          Trees,
          Chm,
          "W_Wood_kg",
          fun = sum,
          na.rm = T,
          filename = paste0(output_campaign_tilesP, tile_id, "_Wwood_kg_sum_",pixel_size,".tif"),
          overwrite = T
        )#Write Wood biomass raster
        W_foliage_sum <- terra::rasterize(
          Trees,
          Chm,
          "W_foliage_kg",
          fun = sum,
          na.rm = T,
          filename = paste0(
            output_campaign_tilesP,
            tile_id,
            "_Wfoliage_kg_sum_",pixel_size,".tif"
          ),
          overwrite = T
        )#Write foliage biomass raster
        C_sum <- terra::rasterize(
          Trees,
          Chm,
          "C_kg",
          fun = sum,
          na.rm = T,
          filename = paste0(
            output_campaign_tilesP,
            tile_id,
            "_Carbon_kg_sum_",pixel_size,".tif"
          ),
          overwrite = T
        )#Write carbon raster
        
        
        # if(process_cbh==T){
        #   
        # }
        
      }
    }else{
      values(Chm) <- NA
      for (attribute in ForestAttributes[6:length(ForestAttributes)]) {
        Empty_Raster <- Chm
        #Otherwise write empty rasters
        if(attribute=="foresttype"){
          Empty_Raster <- aggregate(Empty_Raster,30/pixel_size)
        }
        writeRaster(
          Empty_Raster,
          paste0(
            output_campaign_tilesP,
            tile_id,
            "_",
            attribute,"_",pixel_size,
            ".tif"
          ),
          overwrite = T
        )
      }
    }
  }else{
    values(Chm) <- NA
    for (attribute in ForestAttributes) {
      Empty_Raster <- Chm
      #Otherwise write empty rasters
      if(attribute=="foresttype"){
        Empty_Raster <- aggregate(Empty_Raster,30/pixel_size)
      }
      
      writeRaster(
        Empty_Raster,
        paste0(
          output_campaign_tilesP,
          tile_id,
          "_",
          attribute,"_",pixel_size,
          ".tif"
        ),
        overwrite = T
      )
    }
  }
  }
}

#Run extraction function -------------------------------------------------------
TilesDone <- str_sub(list.files(output_campaign_tilesP, "Carbon_kg_sum"),
                       end = -22) #List Tiles already processed

if (length(TilesDone) >= 1) {
  #If some tiles are already processed ...
  TilesToDo <- TilesID[-c(which(TilesID %in% TilesDone))] #List Tiles to be processed
} else{
  TilesToDo <- TilesID #Otherwise, All Tiles should be processed
}

plan(sequential)#Set a sequential processing plan
plan(multisession, workers = core_nb)#Set a parallel processing plan
t1 <- Sys.time() #Get starting time
TEST <- future_map(TilesToDo, process.tile)#Run extraction function
print(c("Processing time : ", Sys.time() - t1))
plan(sequential)

#Merge tiles into final campaign map -------------------------------------------  roi <- "NCentralNM_2016" #Name of the LiDAR mission
  
  CHM1 <- map(list.files(output_campaign_tilesP, "chm1_fill.tif$", full.names = T),
            rast)#Read all tiles for CHM 1m
Chm1 <- mosaic(
  sprc(CHM1),
  filename = paste0(output_campaign_mapP, roi, "_", campaign, "_chm1.tif"),
  overwrite = T
)#Merge Tiles and write final map
for (attribute in ForestAttributes) {
  TILES <- map(list.files(output_campaign_tilesP, attribute, full.names = T),
               rast)#Read all tiles
  Map <- mosaic(
    sprc(TILES),
    fun="min",
    filename = paste0(
      output_campaign_mapP,
      roi,
      "_",
      campaign,
      "_",
      attribute,
      ".tif"
    ),
    overwrite = T
  )#Merge Tiles and write final map
}

print("Processing time : ", Sys.time() - T1)
#Take 25 Minutes for Tile Processing and 5 for Map processing for 32 Tiles
#So more or less 1 minute per tile with 20 cores
