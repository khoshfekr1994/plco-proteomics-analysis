# PLCO Comprehensive Proteomics Analysis Pipeline

A unified R pipeline for comprehensive protein and peptide-level proteomics analysis of PLCO (Prostate, Lung, Colorectal and Ovarian Cancer Screening Trial) mass spectrometry data.

## 🔬 Overview

This pipeline provides a comprehensive framework for analyzing proteomics data from the PLCO Cancer Screening Trial, supporting both protein and peptide-level analyses across IgB (Immunoglobulin B) and FT (Flow-Through) mass spectrometry fractions.

### Key Features

- **🧬 Multi-Level Analysis**: Both protein and peptide-level analysis capabilities
- **📊 Dual Data Types**: Support for IgB and FT mass spectrometry partitions
- **🔄 Flexible Missing Value Handling**: Zero-assignment and NA-assignment strategies
- **📈 Comprehensive Statistical Analysis**: ROC curves, AUC calculations, fold change analysis
- **🎯 Multiple Cancer Types**: Analysis across 6 cancer types (Bladder, Lung, Ovarian, Breast, Colorectal, Pancreatic)
- **👥 Clinical Variables**: Integration of demographic and clinical characteristics
- **🧬 Breast Cancer Subtypes**: Molecular subtype classification (Luminal, HER2, TNBC, Other)
- **⚡ Parallel Processing**: Multi-core support for large-scale analyses
- **📁 Multiple Output Formats**: Excel and RData outputs with automatic file splitting
- **🔧 Configurable Scenarios**: 10 different case-control comparison scenarios

## 🛠 Installation

### Clone the Repository

**Option 1: Using Git (Recommended)**
```bash
# Clone the repository
git clone https://github.com/khoshfekr1994/plco-proteomics-analysis.git

# Navigate to the project directory
cd plco-proteomics-analysis
```

**Option 2: Download ZIP**
1. Click the green "Code" button on the GitHub repository page
2. Select "Download ZIP"
3. Extract the ZIP file to your desired location
4. Navigate to the extracted folder

### Setup

1. **Open R or RStudio** and set your working directory:
```r
# Set working directory to the cloned repository
setwd("/path/to/plco-proteomics-analysis")

# Or in RStudio: File > Open Project > select the folder
```

2. **Install Required Packages**
```r
# Check and install required packages
required_packages <- c("dplyr", "pROC", "parallel", "tidyr", "writexl")

# Install missing packages
missing_packages <- required_packages[!required_packages %in% installed.packages()[,"Package"]]
if(length(missing_packages)) install.packages(missing_packages)

# Load packages to verify installation
lapply(required_packages, library, character.only = TRUE)
```

3. **Configure File Paths**
Open `PLCO_Proteomics_Analysis.R` and update the `CONFIG` section with your file paths:
```r
CONFIG <- list(
  # UPDATE THESE PATHS FOR YOUR SYSTEM
  data_dir = "/path/to/your/data/files/",
  output_dir = "/path/to/desired/output/directory/",
  # ... other settings
)
```

4. **Verify Installation**
```r
# Source the main script
source("PLCO_Proteomics_Analysis.R")

# You should see the startup message if everything is loaded correctly
```

## 📋 Requirements

### R Dependencies
- `dplyr` - Data manipulation
- `pROC` - ROC curve analysis
- `parallel` - Parallel processing
- `tidyr` - Data tidying
- `writexl` - Excel file output

### System Requirements
- R version ≥ 4.0.0
- Minimum 8GB RAM (16GB+ recommended for large datasets)
- Multi-core CPU for parallel processing (optional but recommended)
- Git (optional, for cloning repository)

## 📂 Data Requirements

### Input File Structure

The pipeline expects CSV files with the following structure:

**Protein-level data:**
```
Subject_ID, age, cig_stat, j_breast_er_status, j_breast_her2summ, j_breast_pr_status, 
cancer_type, Gender, Stage, [additional_metadata], Protein1, Protein2, ..., ProteinN
```

**Peptide-level data:**
```
Subject_ID, age, cig_stat, j_breast_er_status, j_breast_her2summ, j_breast_pr_status,
cancer_type, Gender, Stage, [additional_metadata], PG.ProteinNames, PG.ProteinGroups,
EG.PrecursorId, Quantity, [additional_peptide_metadata]
```

### Required Data Files

Place your data files in the configured data directory:

**For Protein Analysis:**
- `protein_IgB_PLCO.csv`
- `protein_FT_PLCO.csv`

**For Peptide Analysis:**
- `peptide_IgB_PLCO.csv`

### Clinical Variables

- `age`: Subject age
- `cig_stat`: Smoking status (0=Never, 1=Current, 2=Former)
- `j_breast_er_status`: Estrogen receptor status (for breast cancer)
- `j_breast_pr_status`: Progesterone receptor status (for breast cancer)  
- `j_breast_her2summ`: HER2 status (for breast cancer)
- `cancer_type`: Cancer diagnosis or "Non-Case"
- `Gender`: Subject gender (F/M)
- `Stage`: Cancer stage (1-4)

## 🚀 Usage

### Quick Start

```r
# Source the main script
source("PLCO_Proteomics_Analysis.R")

# Run full analysis with default settings
run_plco_analysis()

# Run protein-level analysis only
run_plco_analysis(
  analysis_level = "protein",
  data_types = c("IgB", "FT", "SUM")
)

# Run peptide-level analysis for specific cancer
run_plco_analysis(
  analysis_level = "peptide", 
  data_types = "IgB",
  cancer_type = "Breast"
)
```

### Configuration

Modify the `CONFIG` list at the top of the script to customize file paths and parameters:

```r
CONFIG <- list(
  # File paths - UPDATE THESE FOR YOUR SYSTEM
  data_dir = "/path/to/your/data/",
  output_dir = "/path/to/output/",
  
  # Analysis parameters
  cut_off_95 = 0.95,
  cut_off_86 = 0.86,
  protein_column = 15,
  
  # Cancer types to analyze
  cancer_types = c("Bladder", "Lung", "Ovarian", "Breast", "Colorectal", "Pancreatic"),
  scenarios = 1:2  # Extend to 1:10 for all scenarios
)
```

### Command Line Usage

```bash
# Analyze specific cancer type
Rscript PLCO_Proteomics_Analysis.R "Breast" "protein" "IgB"

# Arguments: cancer_type, analysis_level, data_type
```

## 📊 Analysis Scenarios

The pipeline supports 10 different case-control comparison scenarios:

| Scenario | Cases | Controls | Stage |
|----------|-------|----------|--------|
| 1 | Cancer X (all stages) | Non-Case subjects | All |
| 2 | Cancer X (all stages) | All other cancers | All |
| 3 | Cancer X Stage 1 | Non-Case subjects | 1 |
| 4 | Cancer X Stage 1 | All other cancers | 1 |
| 5 | Cancer X Stage 2 | Non-Case subjects | 2 |
| 6 | Cancer X Stage 2 | All other cancers | 2 |
| 7 | Cancer X Stage 3 | Non-Case subjects | 3 |
| 8 | Cancer X Stage 3 | All other cancers | 3 |
| 9 | Cancer X Stage 4 | Non-Case subjects | 4 |
| 10 | Cancer X Stage 4 | All other cancers | 4 |

## 📈 Output Description

### Generated Metrics

For each protein/peptide-cancer-scenario combination, the pipeline calculates:

**Performance Metrics:**
- Area Under the Curve (AUC) with confidence intervals
- Sensitivity and specificity at 95% and 86% cutoffs
- True/False positive/negative rates
- Accuracy and precision metrics

**Statistical Tests:**
- Fold change between cases and controls
- t-test p-values
- ROC curve p-values

**Demographics:**
- Age statistics for cases and controls  
- Smoking status distributions
- Breast cancer molecular subtypes (when applicable)

**Missing Data:**
- Missing value counts and percentages
- Valid sample sizes

### Output Files

Results are saved in organized directories by date:
```
output_dir/
├── protein/
│   ├── zero_assigned/
│   │   └── November_26/
│   │       ├── results_protein_IgB_zero_assigned.xlsx
│   │       ├── results_protein_FT_zero_assigned.xlsx
│   │       ├── results_protein_SUM_zero_assigned.xlsx
│   │       └── *.RData files
│   └── NA_assigned/
└── peptide/
    ├── zero_assigned/
    └── NA_assigned/
```

## 🔧 Advanced Configuration

### Missing Value Strategies

**Zero Assignment (`"zero_assigned"`):**
- Replaces NA values with 0
- Treats 0 as missing for statistical calculations
- Suitable for datasets where missing = not detected

**NA Assignment (`"NA_assigned"`):**  
- Keeps NA values as missing
- Excludes missing values from calculations
- Suitable for datasets with true missing data

### Parallel Processing

Enable parallel processing for large datasets:

```r
run_plco_analysis(
  analysis_level = "protein",
  data_types = c("IgB", "FT"), 
  use_parallel = TRUE  # Uses all available cores - 1
)
```

### Custom Analysis

For specialized analysis, you can call individual functions:

```r
# Load specific dataset
data <- load_data("path/to/file.csv", missing_strategy = "zero_assigned")

# Process single combination  
result <- process_protein_combination(data, params_row, missing_strategy)

# Calculate ROC metrics
roc_metrics <- calculate_roc_metrics(data, "protein_name")
```

## 📊 Interpretation Guide

### AUC Values
- **AUC > 0.8**: Excellent discrimination
- **AUC 0.7-0.8**: Good discrimination  
- **AUC 0.6-0.7**: Fair discrimination
- **AUC < 0.6**: Poor discrimination

### Fold Change
- **FC > 2**: Strong upregulation in cases
- **FC 1.5-2**: Moderate upregulation
- **FC 0.5-1.5**: Minimal change
- **FC < 0.5**: Downregulation in cases

### Statistical Significance
- **p < 0.001**: Highly significant
- **p < 0.01**: Very significant
- **p < 0.05**: Significant
- **p ≥ 0.05**: Not significant

## 🐛 Troubleshooting

### Common Issues

**Memory errors:**
- Reduce the number of proteins analyzed at once
- Enable parallel processing to manage memory better
- Increase system RAM

**File not found errors:**
- Verify file paths in CONFIG section
- Ensure data files are in the correct directory
- Check file naming conventions

**Missing dependencies:**
```r
# Install missing packages
if (!requireNamespace("pROC", quietly = TRUE)) {
  install.packages("pROC")
}
```

**Empty results:**
- Check that cancer types match your data
- Verify column names match expected format
- Ensure sufficient sample sizes for analysis

## 📚 Citation

If you use this pipeline in your research, please cite:

```
Khoshfekr Rudsari, H. (2024). PLCO Comprehensive Proteomics Analysis Pipeline. 
GitHub repository: https://github.com/[your-username]/plco-proteomics-analysis
```


## 👨‍💻 Author & Contact

**Hamid Khoshfekr Rudsari, Ph.D.**
- 📧 Email: hkhoshfekr@mdanderson.org
- 📧 Alternative: khoshfekr1994@gmail.com
- 🏥 Affiliation: The University of Texas MD Anderson Cancer Center


### Development Guidelines
- Follow R coding best practices
- Add appropriate documentation for new functions
- Test changes with sample data
- Update README for new features

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## 🔄 Version History

- **v1.0.0** (November 2024): Initial release with unified pipeline
  - Protein and peptide-level analysis
  - IgB and FT data support  
  - Comprehensive statistical analysis
  - Parallel processing capabilities

---
