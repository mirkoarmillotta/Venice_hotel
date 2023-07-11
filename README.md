These files contain data and code used in the paper "Unveiling Venice's hotel competition networks from dynamic pricing digital markets" 
by Mirko Armillotta, Konstantinos Fokianos and Andrea Guizzardi.

Disclaimer: **The data are available only for research purposes.**

First load the R workspace called "data.Rdata".

To use the code you just need to run the main script file called "hotel.R". 

The following R scripts are loaded:

source("network estimation.R")       # Advance booking Network Autoregressive model (NAR) estimation (Algorithm 1)

source("adja_dist.R")                # Compute adjacency matrix from given probabilities

source("ols.nar.R")                  # Compute OLS of advance booking NAR model

source("genNAR.R")                   # Generate values from advance booking NAR model

source("data format.R")             # Format the data for the analysis
