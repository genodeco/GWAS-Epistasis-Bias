# GWAS-Epistasis-Bias

## Description
This is the code repository for the manuscript "[Inflation of genome-wide association test statistics due to omitted interactions](https://www.biorxiv.org/content/10.1101/2025.11.21.689603v1)". 

## Installation

### Prerequisites
- Python 3.9
- Conda (optional, for managing the environment)

### Setting Up the Environment

#### Using Conda (Recommended)
1. **Clone the repository**:
    ```bash
    git clone https://github.com/genodeco/GWAS-Inflation
    cd GWAS-Inflation
    ```

2. **Create a Conda environment**:
    ```bash
    conda create --name my_env python=3.9
    conda activate my_env
    ```

3. **Install the dependencies**:
    ```bash
    pip install -r requirements.txt
    ```

#### Using Virtualenv
1. **Clone the repository**:
    ```bash
    git clone https://github.com/genodeco/GWAS-Inflation
    cd GWAS-Inflation
    ```

2. **Create a virtual environment**:
    ```bash
    python -m venv venv
    ```

3. **Activate the virtual environment**:
    - On Windows:
        ```bash
        venv\Scripts\activate
        ```
    - On macOS and Linux:
        ```bash
        source venv/bin/activate
        ```

4. **Install the dependencies**:
    ```bash
    pip install -r requirements.txt
    ```

### Running the Scripts
The scripts can be run independently. Please note that simulation analyses were performed with the Estonian Biobank dataset which cannot be shared publicly. Therefore, mock genotypes and covarites are generated in relevant sections. Please use your own genotype and covariate data to obtain meaningful results.


### Citation

```bibtex
@misc{yelmen_bias_2025,
	title = {Bias in genome-wide association test statistics due to omitted interactions},
	url = {https://www.biorxiv.org/content/10.1101/2025.11.21.689603v1},
    doi = {10.1101/2025.11.21.689603},
	publisher = {bioRxiv},
	author = {Yelmen, Burak and Güler, Merve Nur and Team, Estonian Biobank Research and Kollo, Tõnu and Möls, Märt and Charpiat, Guillaume and Jay, Flora},
	date = {2025-11-22},
}
