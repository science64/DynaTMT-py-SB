# DynaTMT

The **DynaTMT tool** can be used to analyze **m**ultiplexed **e**nhanced **pro**tein **d**ynamic mass spectrometry (**mePROD**) data, as well as **mePRODmt** (mitochondrial protein import proteomics) data. mePROD uses pulse SILAC combined with Tandem Mass Tag (TMT) labelling to profile newly synthesized proteins. Through a booster channel, that contains a fully heavy labelled digest, the identification rate of labelled peptides is greatly enhanced, compared to other pSILAC experiments. Through the multiplexing capacity of TMT reagents it is possible during the workflow to use the boost signal as a carrier that improves survey scan intensities, but does not interfere with quantification of the pulsed samples. This workflow makes labelling times of minutes (down to 15min in the original publication) possible.

**mePRODmt** extends this approach to specifically quantify mitochondrial protein import dynamics, allowing researchers to measure the rate of protein translocation into mitochondria under various conditions.

Additionally, mePROD utilizes a baseline channel, comprised of a non-SILAC labelled digest that serves as a proxy for isolation interference and greatly improves quantification dynamic range (note: baseline correction is **required for MS2 data** but **not necessary for MS3 data** due to reduced co-isolation interference in MS3). Quantification values of a heavy labelled peptide in that baseline channel are derived from co-fragmented heavy peptides and will be subtracted from the other quantifications.

The package can also be used to analyse any pSILAC/TMT dataset.

## Version

Current version: **2.9.4** (2026-03-19). See the [CHANGELOG](CHANGELOG.md) for details on what has changed.

## References

If you use DynaTMT in your research, please cite the following publications:

1. **Original mePROD Publication (2020):**
   > Klann K, Tascher G, Münch C. Functional Translatome Proteomics Reveal Converging and Dose-Dependent Regulation by mTORC1 and eIF2α. *Molecular Cell*. 2020;77(4):913-925.e4.
   > [DOI: 10.1016/j.molcel.2019.11.010](https://doi.org/10.1016/j.molcel.2019.11.010)

2. **Mitochondrial Protein Import Proteomics mePRODmt (2022):**
   > Schäfer JA, Bozkurt S, Michaelis JB, Klann K, Münch C. Global mitochondrial protein import proteomics reveal distinct regulation by translation and translocation machinery. *Molecular Cell*. 2022;82(2):435-446.e7.
   > [DOI: 10.1016/j.molcel.2021.11.004](https://www.cell.com/molecular-cell/fulltext/S1097-2765(21)00954-0)

3. **mePRODmt Proteomics and PBLMM Methods Paper (2024):**
   > Bozkurt S, Parmar BS, Münch C. Quantifying mitochondrial protein import by mePRODmt proteomics. *Methods in Enzymology*. 2024;706:449-474.
   > [DOI: 10.1016/bs.mie.2024.07.017](https://pubmed.ncbi.nlm.nih.gov/39455229/)

4. **Injection Time Adjustment and Targeted Mass Differences (TMD) Method:**
   > Klann K, Münch C. Instrument Logic Increases Identifications during Multiplexed Translatome Measurements. *Analytical Chemistry*. 2020.
   > [DOI: 10.1021/acs.analchem.0c01749](https://doi.org/10.1021/acs.analchem.0c01749)

---

## Table of Contents

- [Installation](#installation)
- [Quick Start](#quick-start)
- [MS2 Workflow](#ms2-workflow)
- [MS3 Workflow](#ms3-workflow)
- [Functions Reference](#functions-reference)
- [Loading Data](#loading-data)
- [Proteome Discoverer Setup](#proteome-discoverer-setup)
- [API Documentation](#api-documentation)

---

## Installation

### Requirements

- Python 3.8+
- Git (for installation from GitHub)
- JupyterLab or Jupyter Notebook (recommended for interactive analysis)

### Install from GitHub (Recommended)

```bash
pip install --upgrade git+https://github.com/science64/DynaTMT.git
```

### Install from Source

```bash
# Clone the repository
git clone https://github.com/science64/DynaTMT.git
cd DynaTMT

# Install the package
pip install .
```

### Required Dependencies

```bash
pip install pandas==2.3.3 numpy matplotlib statsmodels scipy
```

### Optional: Install PBLMM for Statistical Analysis

For peptide-based linear mixed model (PBLMM) statistical analysis:

```bash
pip install --upgrade git+https://github.com/science64/PBLMM.git
```

---

## Quick Start

```python
import DynaTMT as mePROD
import pandas as pd

# Load your PSMs file
psms = pd.read_csv("your_PSMs_file.txt", sep='\t', header=0)

# Remove booster channel (adjust column name as needed)
data = psms.drop('Abundance: 131C', axis=True)

# Initialize the processor
process = mePROD.PD_input(data)

# Filter and process
filtered = process.filter_PSMs(data)
normalized = process.total_intensity_normalisation(filtered)
heavy = process.extract_heavy(normalized)

# Continue with MS2 or MS3 specific workflow (see below)
```

---

## MS2 Workflow

**Use this workflow for MS2-based TMT quantitation data.**

MS2-based quantification requires:
- **Injection time (IT) adjustment** - Required to correct for variable ion injection times
- **Baseline correction** - Required to correct for co-isolation interference

### Complete MS2 Example

```python
from datetime import date
import pandas as pd
import DynaTMT as mePROD
import PBLMM

# Configuration
wd = "Example data/MS2_data"
dataName = "20200724_SB_CCCP_ISRIB_Import_PSMs.txt"
conditions = ['Light', 'DMSO', 'DMSO', 'DMSO', 'CCCP', 'CCCP', 'CCCP', 'CCCP_ISRIB', 'CCCP_ISRIB', 'CCCP_ISRIB']
pairs = [['CCCP', 'DMSO'], ['CCCP_ISRIB', 'DMSO'], ['CCCP_ISRIB', 'CCCP']]

# Load and prepare data
psms = pd.read_csv(f'{wd}/{dataName}', sep='\t', header=0)
booster_removed = psms.drop('Abundance: 131C', axis=True)  # Remove booster channel

# Initialize processor
process = mePROD.PD_input(booster_removed)

# Step 1: Filter PSMs (removes contaminants, shared peptides, isolation interference >50%)
filter_data = process.filter_PSMs(booster_removed)

# Step 2: IT adjustment (REQUIRED for MS2)
IT_adjusted = process.IT_adjustment(filter_data)

# Step 3: Normalization
sumNorm = process.total_intensity_normalisation(IT_adjusted)

# Step 4: Extract heavy peptides
heavy = process.extract_heavy(sumNorm)

# Step 5: Baseline correction (REQUIRED for MS2)
peptide_data = process.baseline_correction(heavy, threshold=15, i_baseline=0, random=True)

# Step 6: Statistical analysis with PBLMM
hypothesis_testing = PBLMM.HypothesisTesting()
resultFinal = hypothesis_testing.peptide_based_lmm(peptide_data, conditions=conditions, pairs=pairs)
resultFinal.reset_index(inplace=True)
resultFinal.rename(columns={'index': 'Accession'}, inplace=True)

# Export results
resultFinal.to_excel(f'{wd}/results_MS2_{date.today().strftime("%d.%m.%Y")}.xlsx', index=False)
```

### MS2 Workflow Key Points

| Step | Function | Required | Notes |
|------|----------|----------|-------|
| 1 | `filter_PSMs()` | Yes | Removes contaminants, shared peptides, isolation interference >50% |
| 2 | `IT_adjustment()` | **Yes** | Corrects for variable ion injection times |
| 3 | `total_intensity_normalisation()` | Yes | Or use `Median_normalisation()` or `TMM()` |
| 4 | `extract_heavy()` | Yes | Extracts heavy-labeled peptides |
| 5 | `baseline_correction()` | **Yes** | Subtracts baseline channel; converts PSMs to peptides |
| 6 | PBLMM analysis | Recommended | Peptide-based linear mixed model |

---

## MS3 Workflow

**Use this workflow for MS3-based TMT quantitation data.**

MS3-based quantification has reduced co-isolation interference, therefore:
- **Injection time adjustment** - NOT required (skip this step)
- **Baseline correction** - NOT required (optional); remove baseline channel instead

### Complete MS3 Example

```python
from datetime import date
import pandas as pd
import DynaTMT as mePROD
import PBLMM

# Configuration
wd = "Example data/MS3_data"
dataName = "20240109_LU_LC2_MAA_DB_mePROD_MS3_PSMs.txt"
conditions = ['Cont', 'Cont', 'Cont', 'USP39', 'USP39', 'USP39']
pairs = [['USP39', 'Cont']]

# Load and prepare data
psms = pd.read_csv(f'{wd}/{dataName}', sep='\t', header=0)
booster_removed = psms.drop('Abundance 131C', axis=True)  # Remove booster channel
baseline_removed = booster_removed.drop('Abundance 126', axis=True)  # Remove baseline channel

# Initialize processor
process = mePROD.PD_input(baseline_removed)

# Step 1: Filter PSMs
filter_data = process.filter_PSMs(baseline_removed)

# Step 2: Normalization (NO IT adjustment for MS3!)
sumNorm = process.total_intensity_normalisation(filter_data)

# Step 3: Extract heavy peptides
heavy = process.extract_heavy(sumNorm)

# Step 4: Convert PSMs to Peptides (NO baseline correction for MS3!)
peptide_data = process.PSMs_to_Peptide(heavy)

# Step 5: Statistical analysis with PBLMM
hypothesis_testing = PBLMM.HypothesisTesting()
resultFinal = hypothesis_testing.peptide_based_lmm(peptide_data, conditions=conditions, pairs=pairs)
resultFinal.reset_index(inplace=True)
resultFinal.rename(columns={'index': 'Accession'}, inplace=True)

# Export results
resultFinal.to_excel(f'{wd}/results_MS3_{date.today().strftime("%d.%m.%Y")}.xlsx', index=False)
```

### MS3 Workflow Key Points

| Step | Function | Required | Notes |
|------|----------|----------|-------|
| 1 | `filter_PSMs()` | Yes | Removes contaminants, shared peptides, isolation interference >50% |
| 2 | `IT_adjustment()` | **No** | Skip for MS3 data |
| 3 | `total_intensity_normalisation()` | Yes | Or use `Median_normalisation()` or `TMM()` |
| 4 | `extract_heavy()` | Yes | Extracts heavy-labeled peptides |
| 5 | `PSMs_to_Peptide()` | Yes | Combines PSMs into peptides |
| 6 | `baseline_correction()` | **No** | Skip for MS3 data; remove baseline channel instead |
| 7 | PBLMM analysis | Recommended | Peptide-based linear mixed model |

---

## MS2 vs MS3 Comparison

| Feature | MS2 | MS3 |
|---------|-----|-----|
| IT Adjustment | **Required** | Not needed |
| Baseline Correction | **Required** | Not needed (remove baseline channel) |
| Co-isolation Interference | Higher | Lower |
| Signal Intensity | Higher | Lower |
| PSM to Peptide | Via `baseline_correction()` | Via `PSMs_to_Peptide()` |

---

## Functions Reference

### PD_input Class Functions

| Function | Description |
|----------|-------------|
| `get_channels(input)` | Returns list of unnormalized abundance column names |
| `filter_peptides(input)` | Filters peptide files: removes shared peptides, contaminants, NA values |
| `filter_PSMs(input)` | Filters PSM files: removes shared peptides, contaminants, isolation interference >50%, NA values |
| `IT_adjustment(input)` | Adjusts abundances for ion injection times (**MS2 only**) |
| `total_intensity_normalisation(input)` | Normalizes to total intensity across TMT channels |
| `Median_normalisation(input)` | Normalizes to median of TMT channels |
| `TMM(input)` | Trimmed Mean of M-values normalization |
| `extract_heavy(input)` | Extracts heavy-labeled (SILAC) peptides |
| `extract_light(input)` | Extracts light-labeled peptides |
| `baseline_correction(input, threshold, i_baseline, random)` | Subtracts baseline channel, filters by threshold, converts PSMs to peptides |
| `PSMs_to_Peptide(input)` | Combines PSMs into peptides without baseline correction |
| `protein_rollup(input, method)` | Aggregates peptides to protein level ('sum', 'mean', 'median') |
| `log2(input)` | Log2 transforms all TMT intensities |

---

## Loading Data

### ProteomeDiscoverer Output (Default)

DynaTMT uses ProteomeDiscoverer PSM or Peptide file outputs in tab-delimited text format:

```python
import pandas as pd

# Load PSMs file
psms = pd.read_csv("PSMs.txt", sep='\t', header=0)

# Or load Peptides file
peptides = pd.read_csv("Peptides.txt", sep='\t', header=0)
```

### Plain Text Input

For custom file formats, use the `plain_text_input` class. Column order:
1. Protein Accession
2. Ion Injection Time (optional)
3. Modifications
4. TMT Abundances (all remaining columns)

```python
from DynaTMT.DynaTMT import plain_text_input

df = pd.read_csv("custom_data.txt", sep='\t', header=0)
processor = plain_text_input(df, it_adj=True)  # Set it_adj=False if no injection time column
```

---

## Proteome Discoverer Setup

### Custom Modifications for Heavy Lysine + TMT

Since search engines cannot handle two modifications on the same residue, create custom modifications combining TMT and heavy lysine:

| Modification Name | Description | ΔM (monoisotopic) | ΔM (average) |
|-------------------|-------------|-------------------|--------------|
| `Label:13C(6)15N(4)` | Heavy Arginine (PD default) | +10.0083 | - |
| `TMTK8` | TMT + Heavy Lysine (K8) | +237.1771 | +237.2061 |
| `TMTproK8` | TMTpro + Heavy Lysine (K8) | +312.2213 | +312.2554 |
| `TMTK4` | TMT + Heavy Lysine (K4) | - | - |
| `TMTK6` | TMT + Heavy Lysine (K6) | - | - |
| `TMTproK4` | TMTpro + Heavy Lysine (K4) | - | - |
| `TMTproK6` | TMTpro + Heavy Lysine (K6) | - | - |

### Workflow Setup in Proteome Discoverer

1. Create an MS2 or MS3 reporter ion-based workflow depending on your acquisition method
2. Use SequenceHT node for database searches
3. Set enzyme to trypsin
4. Add TMT as static modification at N-terminus and lysines
5. Add TMTK8/TMTproK8 and Label:13C(6)15N(4) (Arg10) as dynamic modifications
6. Use 1% FDR threshold at peptide and protein levels
7. Use default Percolator settings
8. Export PSMs file in text format

---

## API Documentation

### class PD_input

Class for analyzing ProteomeDiscoverer peptide/PSM output with default column names.

#### `__init__(self, input)`

Initializes the PD_input class with input DataFrame. Automatically extracts TMT abundance channels.

**Parameters:**
- `input` (DataFrame): ProteomeDiscoverer PSM or Peptide output

---

#### `get_channels(self, input)`

Identifies and returns abundance column names from the input DataFrame.

**Returns:**
- List of column names containing abundance data (excludes normalized columns)

---

#### `filter_peptides(self, filtered_input)`

Filters peptide-level data by removing:
- Shared peptides (non-unique quantification)
- Contaminant proteins
- Rows with NA values in abundance channels

**Parameters:**
- `filtered_input` (DataFrame): Peptide data

**Returns:**
- DataFrame: Filtered peptide data

---

#### `filter_PSMs(self, filtered_input)`

Filters PSM-level data by removing:
- Shared peptides (multiple protein accessions)
- Empty protein accessions
- Contaminant proteins
- PSMs with isolation interference >50%
- Rows with NA values in abundance channels

**Parameters:**
- `filtered_input` (DataFrame): PSM data

**Returns:**
- DataFrame: Filtered PSM data

---

#### `IT_adjustment(self, input)`

Adjusts TMT abundances for ion injection times. **Required for MS2 data only.**

Formula: `Abundance_adjusted = (Abundance / Injection_Time) * 1000`

**Parameters:**
- `input` (DataFrame): Data with "Ion Inject Time" column

**Returns:**
- DataFrame: IT-adjusted data

---

#### `total_intensity_normalisation(self, input)`

Normalizes abundances to the total intensity across all TMT channels. The channel with the lowest total intensity is used as reference.

**Parameters:**
- `input` (DataFrame): Data to normalize

**Returns:**
- DataFrame: Normalized data

---

#### `Median_normalisation(self, input)`

Normalizes abundances to the median of each TMT channel. The channel with the lowest median is used as reference.

**Parameters:**
- `input` (DataFrame): Data to normalize

**Returns:**
- DataFrame: Normalized data

---

#### `TMM(self, input)`

Implements Trimmed Mean of M-values (TMM) normalization (Robinson & Oshlack, 2010, Genome Biology).

**Parameters:**
- `input` (DataFrame): Data to normalize

**Returns:**
- DataFrame: TMM-normalized data

---

#### `extract_heavy(self, input)`

Extracts heavy-labeled peptides based on modification strings. Searches for: `TMTK8`, `Label`, `TMTproK8`, `TMTK4`, `TMTK6`, `TMTproK4`, `TMTproK6`

**Parameters:**
- `input` (DataFrame): Data containing Modifications column

**Returns:**
- DataFrame: Heavy-labeled peptides only

---

#### `extract_light(self, input)`

Extracts light-labeled peptides (peptides WITHOUT heavy modifications).

**Parameters:**
- `input` (DataFrame): Data containing Modifications column

**Returns:**
- DataFrame: Light-labeled peptides only

---

#### `PSMs_to_Peptide(self, input)`

Combines PSMs into peptides by grouping on Annotated Sequence, Modifications, and Master Protein Accessions, then summing abundances.

**Parameters:**
- `input` (DataFrame): PSM-level data

**Returns:**
- DataFrame: Peptide-level data

---

#### `baseline_correction(self, input_file, threshold=5, i_baseline=0, random=True)`

Subtracts baseline channel from all other channels and converts PSMs to peptides. **Required for MS2 data.**

**Parameters:**
- `input_file` (DataFrame): Heavy peptide data
- `threshold` (float): Minimum mean signal after baseline subtraction (default: 5)
- `i_baseline` (int): Index of baseline channel (default: 0, typically TMT 126)
- `random` (bool): If True, replaces zero/negative values with random values 0-1 (recommended for PBLMM)

**Returns:**
- DataFrame: Baseline-corrected peptide data

**Notes:**
- Negative values after subtraction are set to zero
- If `random=True`, zero values are replaced with random values (0-1) to avoid zero-inflation in statistical models
- Mean signal threshold filters low-quality peptides
- Automatically converts PSMs to peptides

---

#### `protein_rollup(self, input_file, method='sum')`

Aggregates peptide-level data to protein-level by grouping on Master Protein Accessions.

**Parameters:**
- `input_file` (DataFrame): Peptide-level data
- `method` (str): Aggregation method - 'sum', 'mean', or 'median' (default: 'sum')

**Returns:**
- DataFrame: Protein-level data (proteins as index)

---

#### `log2(self, input)`

Log2 transforms all TMT abundance values.

**Parameters:**
- `input` (DataFrame): Data to transform

**Returns:**
- DataFrame: Log2-transformed data

---

### class plain_text_input

Class for analyzing pSILAC data from custom text files where column identity is determined by position rather than name.

**Column Order:**
1. Protein Accession
2. Ion Injection Time (optional)
3. Modifications
4. All subsequent columns: TMT Abundances

#### `__init__(self, input, it_adj=True)`

**Parameters:**
- `input` (DataFrame): Input data
- `it_adj` (bool): If True, expects injection time in column 2. If False, modifications are in column 2.

**All other methods function identically to PD_input class.**

---

## Jupyter Notebook Templates

The repository includes ready-to-use Jupyter notebooks:

| Notebook | Description |
|----------|-------------|
| `mePROD LMM simplified_MS2.ipynb` | Simplified MS2 workflow with PBLMM |
| `mePROD LMM simplified_MS3.ipynb` | Simplified MS3 workflow with PBLMM |
| `mePROD LMM template.ipynb` | Full template with visualization and QC |

---

## Example Data

Example data is provided in the `Example data/` folder:

- `MS2_data/` - Example MS2 PSM data with conditions and pairs files
- `MS3_data/` - Example MS3 PSM data with conditions and pairs files

---

## Troubleshooting

### Common Issues

1. **Column not found error**: Ensure your data uses standard ProteomeDiscoverer column names, or use `plain_text_input` class.

2. **Empty result after filtering**: Check that your data contains valid values and proper modification annotations.

3. **Baseline correction produces all zeros**: Adjust the threshold parameter to a lower value.

4. **PBLMM installation fails**: Ensure Git is installed and accessible from command line.

---

## License

This project is licensed under the GNU General Public License v3.0 - see the [LICENSE](LICENSE) file for details.

---

## Authors

- **Kevin Klann** - Original development
- **Süleyman Bozkurt** - Maintenance and updates

**Contact:** sbozkurt.mbg@gmail.com

---

## Acknowledgments

This work was developed in the [Münch Lab](https://molsysmed.de/) at Goethe University Frankfurt.
