# gtMAP

gtMAP is a decoder for group testing using quantitative PCR-test results implemented in R and STAN. The associated paper can be found at URL.

## Installation

gtMAP can be installed in RStudio directly from github using the devtools library as follows

```R
install.packages('devtools') # if devtools is not installed
library(devtools)
install_github("lahunruh/gtMAP")
```

After installation the package can be loaded using

```R
library(gtMAP)
```

To check that the installation was successful, run

```R
gtMAP::run_test()
```
which will create a group testing design as well as data and conduct the group testing and decoding. Depending on the complexity of the randomly generated test, this may take a few minutes. The output should look as follows

```
Creating M
Creating X
Creating Y
Running Decoder
Finished Running Successfully!
NULL
```

(THIS INSTALLATION TEST COULD BE MADE A BIT BETTER IN THE FUTURE)

## Usage

The main function used is `run_decoder.R`. This will run the decoder on specified input data. The function takes four parameters:
- Y: The vector of PCR group test results
- M: The design matrix (rows = tests, columns = samples)
- n: The number of samples
- l: The disjuctness of the design. If this is not known/applicable, set l = 0.

The main decoder can then be run via

```R
run_decoder(Y,M,n,l)
```

The decoder can also be run with different parameter settings. For details on this please refer to the publication and the function help.

### Generating a shifted-transversal design

gtMAP also allows the generation of a design matrix M via the function `generate_STD_testing_protocol.R` which takes two required arguments: 
- n: The number of samples
- d: The desired disjunctness of the design matrix

Please ensure that n is a square number. The resulting number of tests is equal to (d+1)*sqrt(n).

### Testing a design

If you wish to supply your own design and test it you can simulate data as follows

```R
X <- create_positives_vector(n,p=p,output=X.output) # Creates a vector of viral loads with X.output = 'VL'. n is the number of samples and p is the desired prevalence.
Y <- determine_Y(M,X,mode = Y.mode) # Creates the group testing results. Y.mode = 'Ct' returns Ct-values as test results. M is the design matrix supplied by the user.
```

## In Vitro Data

The plate per-well viral loads used for in vitro validation can be found in the data folder. The resulting PCR-test results (Ct-values) can also be found there. To extract the in vitro results the following script can be run
(FINISH!!)
```R
n <-  81
l <-  3
t <- sqrt(n) * (l + 1)
M <- generate_STD_testing_protocol(n, l, sq = c(sqrt(n),0,1,2))

p <- 1 # number of plate

assignment_all <- read.csv(paste(assignment_file,sep=""))
result_file <-  paste('plate',p,'_results_',plex,'.csv',sep='')
p1_res <- read.csv(paste(dir,result_file,sep=''))

start_indeces = c(1,1+t,97,97+t,193,193+t,289,289+t)
for (idx in 1:len(start_indeces) {
  idx_start = start_indeces[idx]
  n_plex = ncol(p1_res)
  if (n_plex > 5) {
    n_plex = 5
  }
  for (i in 2:n_plex) {
    p1_res[p1_res[,i] == "No Ct",i] = 0
  }

  names(p1_res)[1] <- 'well'
  plate1 <- assignment_all

  for (i in 2:n_plex) {
    Y <- p1_res[match(plate1[idx_start:(idx_start+t-1),2],  p1_res$well), names(p1_res)[i]]
  }
}
```

## Citation

TO FOLLOW
