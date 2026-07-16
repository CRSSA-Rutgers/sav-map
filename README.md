# sav-map

**READ ME**

**This repository contains the datasets and code for evaluating seagrass occurrence in shallow waters using PlanetScope satellite imagery: in this case it is applied to Barnegat Bay-Little Egg Harbor estuary in New Jersey, looking at Zostera & Ruppia species throughout 2023**

*GENERAL INFORMATION*

1. Title of Dataset: 
“Data for: Assessment of Seagrass Status in the BBLEH Estuary System using Satellite Imagery-Based Modeling, 2023”

2. Corresponding Author Information		
NAME: Jess M. Stitt 
EMAIL: jessmstitt [at] gmail.com / js3585 [at] crssa.rutgers.edu

3. Date of data collection: 	
Satellite data acquisition: 2023-05-26 through 2023-10-16; Field data collection: 2023-05-19 through 2023-10-27

4. Geographic location of data collection: 
Barnegat Bay-Little Egg Harbor Estuary, New Jersey, USA

*SHARING/ACCESS INFORMATION*

1. Link to Technical Report that cites or uses the data: 
   https://crssa.rutgers.edu/projects/sav/downloads.html (in progress)

2. It is strongly recommended that careful attention be paid to the contents of the SAV reports and metadata files associated with these data ('DATA & METADATA', http://crssa.rutgers.edu/projects/coastal/sav/downloads.html and 'SAV MAPPING METHODS', http://crssa.rutgers.edu/projects/coastal/sav/methods.html) . CRSSA, as well as any other contributors as listed in the metadata, shall be acknowledged as data contributors to any reports or other products derived from these data.

3. This code, as presented on this web site, is meant to generate maps to provide a regional picture of SAV distribution at various points in time and are not intended for site level permit applications or litigation purposes. This code and data, alone, are not sufficient to determine the presence or absence of SAV. Conclusive evidence concerning the presence or absence of SAV requires site inspections, preferably at several points in time during the SAV growing season.

4. Licenses/restrictions placed on the data: None

5. Links to other publicly accessible locations of the data: 
   [https://github.com/orgs/CRSSA-Rutgers](https://github.com/orgs/CRSSA-Rutgers/dashboard)

6. Links/relationships to ancillary data sets: NA

7. Was data derived from another source? No

*DATA & FILE OVERVIEW*

1. File List:
   
**1.0.  00_sav25_satmap_WORKFLOW.Rmd** 
        R markdown (.rmd) file containing all R scripts & packages used for analyses for this project, including project structure, data processing, and Random Forest (RF) modeling and evaluation.
		
**1.1.  01_clean_ref-data.R**
		    R script 
		
**1.2.  02_build_img-stack.R**
		    R script 

**1.3.  03_build_dep-stack.R**
		    R script 

**1.4.  04_build_wcc-stack.R**
		    R script 

**1.5.  05_clean_pred-var.R**
		    R script 

**1.6.  06_model_rf.R**
		    R script 

**1.7.  07_clean_comp-maps.R**
		    R script 

*METHODOLOGICAL INFORMATION*

1. Description of methods used for collection/generation of data: 
All methods used for collecting and generating data can be found in the associated technical report linked to this dataset.

2. Methods for processing the data: 
Methods for processing the data can be found in the R Markdown of File 1.0 and are described in the associated publication linked to this dataset.

3. Instrument- or software-specific information needed to interpret the data: 
Data analyses were performed using the software R (v.4.5.2); all necessary packages and code to perform analyses can be found in the R Markdown of File 1.0.
