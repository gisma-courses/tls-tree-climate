# Load necessary libraries
library(ncdf4)
library(ggplot2)

system("cdo selname,T /media/creu/54D84A9D5B30CD28/proj-data/envimet-hasenkopf-simulation/core_plan_v2/core_plan_v2_output/NetCDF/core_plan_v2_2025-07-21_00.00.00.nc ~/out3.nc
")

printout= FALSE

# Path to the NetCDF file
file_path <- "/home/creu/out3.nc"

# Open the NetCDF file
#nc <- nc_open(file_path)
# Extract temperature, altitude, and time variables
#temperature <- ncvar_get(nc, "T")   # Assuming "temperature" is the variable name
altitude <- ncvar_get(nc, "GridsK")         # Assuming "altitude" is the variable name
time <- ncvar_get(nc, "Time")         # Assuming "altitude" is the variable name
# Check the dimensions of temperature (should return: [lat, lon, altitude, time])
dim(temperature)

# Initialize a matrix to store the inversion altitudes and last valid altitude for each grid cell and time step
# Dimensions: [lat, lon, time]
inversion_altitude <- array(NA, dim = c(nrow(temperature), ncol(temperature), length(time),1)) 
last_na_altitude <- array(NA, dim = c(nrow(temperature), ncol(temperature), length(time),1))
altitude_difference <- array(NA, dim = c(nrow(temperature), ncol(temperature), length(time),1))
mean_temperature_between_levels <- array(NA, dim = c(nrow(temperature), ncol(temperature), length(time),1))
inv_temperature_between_levels <- array(NA, dim = c(nrow(temperature), ncol(temperature), length(time),1))
dim(altitude_difference)

# plotting function
pp = function  (df){
  # data <- data.frame(
  #   altitude = altitude,  # Altitude (on the y-axis)
  #   temperature = temp_profile  # Temperature for the selected point (on the x-axis)
  # )
  # # Plot the temperature profile
  ggplot(df, aes(y = temperature, x = altitude)) +
    geom_line(color = "blue") +
    geom_point(color = "red") +
    geom_text(aes(label = round(temperature,5)),  # Add temperature labels (rounded to 1 decimal)
              hjust = -0.5, vjust = -0.5, size = 3, color = "black") +
    
    ggtitle("Temperature Profile at Selected Point") +
    xlab("Altitude (m)") +
    ylab("Temperature (K)") +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5),
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 14)
    ) 
}



# Function to detect temperature inversion and the last NA altitude
detect_inversion <- function(temp_profile, alt_profile) {
  inversion_level <- NA
  last_na_level <- NA
  
  # Iterate through the temperature profile to find inversion
  clean_data <- na.omit(data.frame(altitude, temp_profile))
  for (i in 2:length(temp_profile)) {
    deepest_temp_index <- which.min(clean_data$temp_profile)
    deepest_temp <- clean_data$temp_profile[deepest_temp_index]
    deepest_altitude <- clean_data$altitude[deepest_temp_index]    
    if (is.na(temp_profile[i-1]) && !is.na(temp_profile[i])) {
      last_na_level <- alt_profile[i]  # Store the last valid altitude before NA
    }
  }
  return(list(inversion_altitude = deepest_altitude, surface_altitude = last_na_level))
}
timestart = 8
timeend =8
# Iterate over each latitude, longitude, and time step, and detect inversion
for (lat in 5:195) {
  for (lon in 2:143) {
    for (t in timestart:timeend) {  # You can adjust this to loop through all times if needed
      # Extract the temperature profile and altitude profile for the specific time, latitude, and longitude
      temp_profile <- temperature[lat, lon, , t]
      alt_profile <- altitude  # Assuming altitude is consistent across latitudes for each time
      
      # Detect inversion and the last valid altitude
      result <- detect_inversion(temp_profile, alt_profile)
      
      
      # Check if both inversion_level and last_na_level are not NA
      if (!is.na(result$inversion_altitude) && !is.na(result$surface_altitude)) {
        
        #altitude_difference[lat, lon, t,] <- result$inversion_altitude - result$surface_altitude
        # Find the indices of the inversion_level and last_na_level in the altitude profile
        inversion_index <- which(alt_profile == result$inversion_altitude)
        last_na_index <- which(alt_profile == result$surface_altitude)
        
        # If inversion level and last valid level are valid and have the same altitude profile
        if (length(inversion_index) > 0 && length(last_na_index) > 0) {
          
          
          # Ensure we get the correct range of indices: from last_na_index to inversion_index
          # Sort to handle whether last_na_level is above or below the inversion level
          start_index <- min(inversion_index, last_na_index)
          end_index <- max(inversion_index, last_na_index)
          
          # Subset the temperature profile between these two levels
          # temp_sub <- temp_profile[start_index:end_index+1]
          # cold_air_mean_temp <- mean(temp_sub, na.rm = TRUE)
          
          sub_data <- na.omit(data.frame(altitude[start_index:end_index+1], temp_profile[start_index:end_index+1]))
          names(sub_data) = c("altitude","temperature")
          if (nrow(sub_data) > 1 & (sub_data$temperature[1] <  sub_data$temperature[2])){
            # Calculate altitude differences
            sub_data$altitude_diff <- c(diff(sub_data$altitude), NA)
            
            # Weighted mean calculation
            # Exclude the last row where altitude_diff is NA
            valid_data <- sub_data[!is.na(sub_data$altitude_diff), ]
            
            # Calculate weighted mean temperature
            cold_air_mean_temp <- sum(valid_data$temperature * valid_data$altitude_diff) / sum(valid_data$altitude_diff)
            
            
            
            # Compute the mean temperature between the two levels
            
            inv_temp = temp_profile[end_index+1]
            
            # Store the computed mean temperature for this grid cell and time step
            mean_temperature_between_levels[lat, lon, t,] <- cold_air_mean_temp
            inv_temperature_between_levels[lat, lon, t,] <- inv_temp
            if (printout) {            # Print results for debugging
              print(paste("mean temp at lat", lat, "lon", lon, "time", t, ":", cold_air_mean_temp))        
              print(paste("inversion temp at lat", lat, "lon", lon, "time", t, ":", inv_temp))        
              print(paste("Inversion altitude at lat", lat, "lon", lon, "time", t, ":", result$inversion_altitude))
              print(paste("Last valid altitude at lat", lat, "lon", lon, "time", t, ":", result$surface_altitude))
              df=data.frame(altitude[start_index:(end_index+3)], temp_profile[start_index:(end_index+3)])
              names(df)=c("altitude","temperature")
              gplot = pp(df)
              print(gplot)
            }
            break()} 
          
          
        }
      }
    }}}

# Close the NetCDF file
nc_close(nc)

# Convert the array for the first time step to a raster
#plot(rast(altitude_difference[, , timestart,]))
#plot(rast(inversion_altitude[, , timestart,])) 
r= rast(mean_temperature_between_levels[, , timestart,])
plot(t(flip(r, direction = "horizontal",)))
r = rast(inv_temperature_between_levels[, , timestart,])
plot(t(flip(r, direction = "horizontal",)))
r = min(rast(inv_temperature_between_levels[, , timestart,])-rast(mean_temperature_between_levels[, , timestart,]))
plot(t(flip(r, direction = "horizontal",)))
r = min(rast(inv_temperature_between_levels[, , timestart,])-rast(mean_temperature_between_levels[, , timestart,]))
plot(t(flip(r, direction = "horizontal",)))

#### tipping point

# identify_tipping_point <- function(tempfile) {
#   # Reverse the tempfile to search from the deepest point toward the surface
#   reversed_temp <- rev(temp_profile)
#   
#   # Iterate through the reversed temperature profile
#   for (i in 2:length(reversed_temp)) {
#     # Ensure both values are non-NA
#     if (!is.na(reversed_temp[i]) && !is.na(reversed_temp[i - 1])) {
#       # Check if temperature starts increasing (from bottom to top in reversed_temp)
#       if (reversed_temp[i - 1] < reversed_temp[i]) {
#         # Return the original index of the tipping point (convert from reversed index)
#         return(length(tempfile) - i + 1)
#       }
#     }
#   }
#   
#   return(NA)  # Return NA if no tipping point is found
# }
# 
# tipping_point_index <- identify_tipping_point(temp_profile)
# 
# # Output the tipping point and its temperature value
# if (!is.na(tipping_point_index)) {
#   print(paste("Tipping point index:", tipping_point_index))
#   print(paste("Temperature at tipping point:", tempfile[tipping_point_index]))
# } else {
#   print("No tipping point found.")
# }