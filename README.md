These files include the code used in the paper "Unveiling Venice's hotel competition networks from dynamic pricing digital markets" by Mirko Armillotta, Konstantinos Fokianos and Andrea Guizzardi.

**Disclaimer: Data are available only for academic purposes and can be obtained from the authors after asking for permission.**


First load the data in the R workspace called "data.Rdata". 


When you use your own data, this should be formatted as an R data frame object, called "Data", where each row is a posted offer and the columns have features in the following order: "id_hotel" "data_in" "adv.book" "Price" "latitude" "longitude" "Stars" "Restaurant" "Meeting" "Rooms" "BAR"

"id_hotel": id of the hotel

"data_in": date of tourist arrival for the posted price offer

"adv.book": integer number >= 0 indicating how many days of advance booking the posted price offer refers to

"Price": the price for the posted offer

"latitude": latitude of the hotel location

"longitude": longitude of the hotel location

"Stars": hotel's number of stars

"Restaurant": dummy equal 1 if restaurant is present, 0 otherwise

"Meeting": dummy equal 1 if meeting room is present, 0 otherwise

"Rooms": hotel's total number of rooms

"BAR": log of "Price"

The row order is instead temporal, starting from the most remote offer to the most recent one.

Then, to use the code you just need to run the main script file called "hotel.R". **When running the code it is important to keep the same file names and the same row/column order shown here.**

The following R scripts are loaded:

source("network estimation.R") # Advance booking Network Autoregressive model (NAR) estimation (Algorithm 1)

source("adja_dist.R") # Compute adjacency matrix from given probabilities

source("ols.nar.R") # Compute OLS of advance booking NAR model

source("genNAR.R") # Generate values from advance booking NAR model

source("data format.R") # Format the data for the analysis
