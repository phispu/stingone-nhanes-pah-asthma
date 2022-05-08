# Urinary metabolites of polycyclic aromatic hydrocarbons and short-acting beta agonist or systemic corticosteroid asthma medication use within NHANES 

***

**Citation:** We are in the process of revising a manuscript and will include a citation here when it is available. 

## 1. Overview
We evaluated the association between 1-hydroxypyrene (a metabolite of polycyclic aromatic hydrocarbons) and 30-day short-acting beta agonist or systemic corticosteroid use (indicator for asthma symptoms).

## 2. Requirements
We used the following software, R packages, and data to complete our analysis.

### 2.1 Software and R Packages
1. Download the following software: 
- [R](https://cran.r-project.org/bin/windows/base/)
- [RStudio](https://www.rstudio.com/products/rstudio/download/#download) or another R graphical user interface
2. Run the following code in R to install the following packages:
- These required packages are needed to run our code. 
	```installation	
	install.packages(c(''rio','tidyverse','survey','arsenal','zoo','here','janitor','broom'), dependencies = TRUE)
	```
3. We used the following versions of software and packages:
- **Software**:
	- *R:* 4.1.1 ("Kick Things")
	- *RStudio:* 2021.09.0+382 ("Ghost Orchid")
- **Packages**:
	- *`rio`*: 0.5.29 
	- *`tidyverse`*: 1.3.1 
	- *`survey`*: 4.1.1 
	- *`arsenal`*: 3.6.3 
	- *`zoo`*: 1.8.9 
	- *`here`*: 1.0.1 
	- *`janitor`*: 2.1.0 
	- *`broom`*: 0.8.0 

### 2.2 Data
We used the following data from multiple waves (2005-2016) of the [U.S. National Health and Nutritional Exam Survey](https://wwwn.cdc.gov/nchs/nhanes/) (document file names in parentheses).

- **Demographics** (DEMO)
- **Questionnaire**: asthma diagnosis (MCQ), prescription medications (RXQ_RX)
- **Laboratory**: urinary polycyclic aromatic hydrocarbon (PAH), urinary creatinine (ALB_CR), serum cotinine (COT)
- **Examination**: BMI (BMX)
- **Health Insurance** (HIQ)

Our final analytical dataset `nhanes-pah-asthma-analysis.RDS` can be found in [`data/processed`](https://github.com/phispu-columbia/stingone-nhanes-pah-asthma/tree/main/data/processed).

- We have included a data dictionary `nhanes-pah-asthma-analysis-DD.xlsx` in the same folder.

## 3. Code and Instructions
We used the following code to complete our analysis: 

- `code/clean_data.R`: Merge, clean, and prepare analytical dataset.
- `code/analysis.R`: Complete analyses.

## 4. Cloning this Repository with RStudio
Below are steps to clone this repository to your local device with RStudio. Please refer to this [link](https://resources.github.com/github-and-rstudio/) for more information about using git in RStudio.

1. On top this page, click on `Code` and copy the link to this git repository (starts with https://github.com/...).
2. Open RStudio.
3. In RStudio, click on `File` &rarr; `New Project...` &rarr; `Version Control` &rarr; `Git`.
4. Under "Repository URL", paste the link of the git repository.
5. Under "Project directory name", name your project directory.
6. Under "Create project as subdirectory of:", select the folder in which you would like your project directory to be located.
7. Click on `Create Project` when you are done to clone your repository! This should take a minute or two to complete.

## 5. Grant Information
This research was supported by grants from the National Institute of Environmental Health Sciences, National Institutes of Health (Grants: #T32ES007322, #R00ES027022).