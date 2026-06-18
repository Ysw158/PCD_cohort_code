# README: PCD_cohort_code

This repository contains R scripts for the analyses presented in the associated manuscript.

## 1. Purpose
- Statistical analysis of clinical phenotypic data.
- Correlation analysis between genetic variants and clinical phenotypes.
- Association tests using Fisher's Exact Test.
- Firth's bias-reduced logistic regression for rare events or small sample sizes.

## 2. Environment and Dependencies
- **OS**: macOS (Darwin 20, aarch64)
- **R version**: 4.4.1 (2024-06-14)
- **The following R packages are required to run the code:**

- dplyr
- tidyr
- tidyverse
- tableone
- logistf
- meta
- forestplot
- epiR
- ggplot2
- ggpubr
- ggsci
- ggrepel
- ggspatial
- ggpie
- cowplot
- RColorBrewer
- scales
- grid
- openxlsx
- eoffice
- sf
- maps
- mapproj
- hchinamap
- sjmisc
- showtext

Install them using:
```r
install.packages(c("dplyr", "tidyr", "tidyverse", "tableone", "logistf", 
                   "meta", "forestplot", "epiR", "ggplot2", "ggpubr", 
                   "ggsci", "ggrepel", "ggspatial", "ggpie", "cowplot", 
                   "RColorBrewer", "scales", "grid", "openxlsx", "eoffice",
                   "sf", "maps", "mapproj", "hchinamap", "sjmisc", "showtext"))
```

## 3. Repository Structure
```
©À©¤©¤ DATA/
©¦ ©À©¤©¤ China.geojson # China map data
©¦ ©À©¤©¤ China_line.geojson # China map boundary data
©¦ ©¸©¤©¤ TEST.xlsx # Test dataset
©À©¤©¤ Rscript/
©¦ ©À©¤©¤ 1.clinic_symptoms_statistic.R # Clinical symptoms statistics
©¦ ©À©¤©¤ 2.Map.R # Map visualization
©¦ ©À©¤©¤ 3.GeneRosepie.R  # Gene rose pie chart
©¦ ©À©¤©¤ 4.Fisher_test_forestplot.R  # Fisher's test & forest plot
©¦ ©¸©¤©¤ 5.Firth_regression.R # Firth logistic regression
©À©¤©¤ PCD_cohort_code.Rproj  # RStudio project file
©À©¤©¤ README.md # This file
©¸©¤©¤ LICENSE # MIT License
```

## 4. How to Run the Code
1. Clone or download this repository to your local machine.
2. Open R or RStudio and set your working directory to the `PCD_cohort_code` folder.
3. Navigate to the `Rscript/` folder and open the main script.
4. Run the scripts in the appropriate order.
5. Outputs (figures, tables) will be saved to the specified output path within the scripts.

## 5. Input Data
- Due to privacy use restrictions, the original data cannot be shared.
- The files in the `DATA/` folder are **test data** only, provided to demonstrate the code structure. They do not represent real results.

## 6. License
This code is licensed under the  MIT License. See the `LICENSE` file for details.

## 7. Contact
For questions, please contact: 22418163@zju.edu.cn