# invasive-grass-disease
 Code and data for field experiment and modeling of invasive grass and disease emergence.  
 
 To run the code and create tables and figures from the manuscript, follow these steps:  
 1. Download entire repository  
 2. Open the RStudio project invasive-grass-disease.Rproj using RStudio 
 3. Run the script install_packages.R to install all packages needed  
 4. Run the scripts in code/data-processing-1  
 5. Run the scripts in code/data-processing-2  
 6. Run the scripts in code/model-fitting  
 7. Run the scripts in code/figure-prep  
 8. Run the scripts in code/figures-tables  
 9. As an alternative to steps 4-8, if you are interested in recreating a specific table or figure, you can find the script associated with that table or figure in code_data_relationships.csv using the "Outputs" columns and only run the script in the "Code" columns.  
 
 It is recommended to restart the R session (unload packages and clear environment) between scripts. You will not need to run the scripts in the directory code/dynamical-models because these are called by other scripts. Scripts in the directory code/imagej were used to quantify leaf disease severity in the software ImageJ.      

Output models, tables, and figures in this respository may vary slightly from those in the manuscript due to model estimation methods. To see models, tables, and figures exactly matching those in the manuscript as well as version history of scripts, visit https://github.com/aekendig/microstegium-bipolaris. Notes that this repository has many more scripts and is not organized and archived for reproducibility.
