# TNBC (Triple Negative Breast Cancer) - Repository

Triple-Negative Breast Cancer (TNBC) is a distinct subtype of breast cancer characterized by the absence of estrogen receptor (ER), progesterone receptor (PR), and human epidermal growth factor receptor 2 (HER2) expression.
This repository is dedicated to supporting computational approaches for TNBC research, including molecular data analysis, automated ADMET property prediction, and efficient data processing workflows. The core focus is on identifying potential natural drug candidates against TNBC using advanced computational methods.


**Note:** The file listing below may be incomplete. To see the full list, please visit the [GitHub UI contents page](https://github.com/Ankita-Kulkarni-134/TNBC-Triple-Negative-breast-cancer-/tree/main).

## Table of Contents

- [Overview](#overview)
- [Repository Structure](#repository-structure)
- [Getting Started](#getting-started)
- [Scripts Summary](#scripts-summary)
- [Usage](#usage)
- [License](#license)

---

## Overview

Triple Negative Breast Cancer (TNBC) is a subtype of breast cancer that lacks expression of ER, PR, and HER2. This repository aims to facilitate computational approaches to analyze molecular data, automate ADMET property calculations, and streamline data processing for TNBC research. Focuses on finding the natural drugs against the TNBC by using the computational methods 

## Repository Structure

- `4For_text_files.py` - Text file processing utilities.
- `5Admet_Automation_with_SSL_Bypass2.py` - Automates ADMET (Absorption, Distribution, Metabolism, Excretion, Toxicity) analysis with SSL bypass.
- `6Combinecsv.py` - Combines multiple CSV files for downstream processing.
- `Clean_Invalid_Molecules.py` - Filters out invalid molecules from datasets.
- `SimilarityTNBC.py` - Calculates molecular similarity relevant to TNBC.
- `Vina_res_csv.py` - Processes results from AutoDock Vina into CSV format.
- `admet_filters.py` - Additional ADMET filtering utilities.
- `impactdrugs.py` - Downloading the natural compunds from the IMPACT database
- `merge_text_files.py` - Merges multiple text files for batch processing.
- `sdftosmi.py` - Converts SDF format to SMILES.
- `supernatural.py` - Handles data from supernatural compound libraries.

> For the full and latest list of files, visit the [repository contents page](https://github.com/Ankita-Kulkarni-134/TNBC-Triple-Negative-breast-cancer-/tree/main).

## Getting Started

1. **Clone the repository:**
   ```bash
   git clone https://github.com/Ankita-Kulkarni-134/TNBC-Triple-Negative-breast-cancer-.git
   cd TNBC-Triple-Negative-breast-cancer-
   ```

2. **Set up your environment:**
   - Requires Python 3.7 or above.
   - Install dependencies as required by individual scripts (see script headers for details).

## Scripts Summary

- **Molecule Preparation:** `Clean_Invalid_Molecules.py`, `sdftosmi.py`
- **Data Processing:** `4For_text_files.py`, `merge_text_files.py`, `6Combinecsv.py`
- **ADMET Analysis:** `5Admet_Automation_with_SSL_Bypass2.py`, `admet_filters.py`
- **Similarity & Analysis:** `SimilarityTNBC.py`, `impactdrugs.py`, `supernatural.py`
- **Docking Results:** `Vina_res_csv.py`

## Usage

Each script can be run independently. Please refer to comments within each script for specific usage instructions and required input data formats.

Example:
```bash
python Clean_Invalid_Molecules.py
```

## License

Distributed under the MIT License. See [LICENSE](LICENSE) for more information.

---

**Disclaimer:** This repository is for research and educational purposes related to computational biology and cheminformatics. Results should be verified and validated before use in clinical or production settings.
