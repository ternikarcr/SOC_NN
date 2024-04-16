Hello! Welcome to my code repository for using memory based learning metrices for predicting SOC using VNIR lab spectra
Please work through it and ensure you understand the content. 

Prerequisites: You need to have R and RStudio installed along with internet connection to download the dependencies. Also the datasets required need to be downloaded from the open source websites

The code is written in R programming language and was tested on the R version i.e. R 4.3.2 -- "Eye Holes"
The analyses were carried out on an Intel Xeon-5320 CPU system having 3.40 GHz processor, 26 cores, 512 GB RAM and NVIDIA GeForce RTX 4060Ti graphics card.
1) Initially load all the necessary packages. SOme may not be available from the standard repositories. For example, the NPRED package may be downloaded from "https://www.hydrology.unsw.edu.au/download/software/npred" and then installed in R environment.
2) Load all the necessary functions 
3) Setup the needed working directory as per your personal file system
4) Download the datasets from the open source websites provided in the manuscripts
5) Load the datasets and format it as per R's datastructure for suitable use
6) Further in terms of steps to be followed, the data is split in cal-val and necessary operations are carried out
7) Intervention may be needed for saving the datasets and renaming the outputs at "write.csv" option and saving the generated images
8) If one is interested in more analysis, the R workspace could be saved for quick access. 


Author: Ternikar Chirag Rajendra
Enjoy! 