# DeepCRISPR-Typer

DeepCRISPR-Typer is a comprehensive CRISPR-Cas system classification tool that combines deep learning models with traditional bioinformatics approaches for accurate CRISPR array detection and Cas protein classification.

## Features

- **CRISPR Array Detection**: Uses MinCED for precise CRISPR array identification
- **Repeat Sequence Classification**: Deep learning-based repeat subtype classification using CNN-Attention architecture
- **Cas Protein Classification**: ESM-2 based protein family classification with LoRA fine-tuning
- **HMM Scanning**: Family-specific HMM scanning for detailed protein domain analysis
- **Rule-based Integration**: Expert rules for final system classification combining multiple evidence sources
- **Parallel Processing**: Optimized pipeline with concurrent execution for faster processing

## Installation

### Prerequisites

Before installing DeepCRISPR-Typer, ensure you have the following tools installed on your system:

1. **MinCED**: CRISPR array detection tool
   - Download from: https://github.com/ctSkennerton/minced
   - Required file: `minced.jar`

2. **Prodigal**: Protein coding gene prediction
   - Download from: http://prodigal.ornl.gov/
   - Should be available in your system PATH or specify custom path

3. **HMMER**: Profile HMM scanning
   - Download from: http://hmmer.org/
   - Should be available in your system PATH or specify custom path

4. **Python**: Version 3.7 or higher
5. **Java**: Required for running MinCED

### Python Dependencies

Install the required Python packages:

```bash
pip install -r requirements.txt
```

### Model and Data Setup

DeepCRISPR-Typer requires pre-trained models and reference databases. You need to obtain these from the following sources:

#### Pre-trained Models

1. **CNN-Attention Model** (`cnn_att_large.pth`)
   - Source: CRISPRclassify-CNN-Att (https://github.com/Xingyu-Liao/CRISPRclassify-CNN-Att)
   - Place in: `data/models/cnn_att_large.pth`

2. **TEMC-Cas Model** (`temc_cas_best_model.pth`) - Optional
   - Source: TEMC-Cas (https://github.com/Xingyu-Liao/TEMC-Cas)  
   - Place in: `data/models/temc_cas_best_model.pth`

#### CRISPRclassify-CNN-Att Source Code

DeepCRISPR-Typer requires access to the CRISPRclassify-CNN-Att source code to load the CNN-Attention model architecture. 

- **Source**: Download the CRISPRclassify-CNN-Att repository from its official source
- **Setup**: Add the `source` directory from CRISPRclassify-CNN-Att to your Python path, or place it in a location where DeepCRISPR-Typer can access it
- **Alternative**: The pipeline will attempt to dynamically load the required modules if the source code is available in the expected locations

### Reference Databases

The following files should be obtained from CRISPRCasTyper project:

1. **HMM Profiles Database**
   - Source: https://github.com/Russel88/CRISPRCasTyper
   - Place entire `Profiles` directory in: `data/Profiles/`

2. **Adaptation Rules** (`adaptation.json`)
   - Source: CRISPRCasTyper project
   - Place in: `data/adaptation.json`

3. **Interference Rules** (`interference.json`)
   - Source: CRISPRCasTyper project  
   - Place in: `data/interference.json`

4. **Cas Scoring Matrix** (`CasScoring.csv`)
   - This file should be included in the repository root
   - If missing, please contact the authors

### Directory Structure

After setup, your directory structure should look like this:

```
DeepCRISPR-Typer/
├── config/
│   ├── model_config.py
│   └── paths_config.py
├── data/
│   ├── input/          # Input FASTA files go here
│   ├── models/
│   │   ├── cnn_att_large.pth
│   │   └── temc_cas_best_model.pth (optional)
│   ├── output/         # Results will be saved here
│   ├── Profiles/       # HMM database
│   ├── adaptation.json
│   └── interference.json
├── src/
├── main.py
├── requirements.txt
├── README.md
└── CasScoring.csv
```

## Configuration

Edit the configuration file `config/paths_config.py` to match your system setup:

```python
# External tool paths - update these to match your installation
MINCED_PATH = "/path/to/your/minced.jar"        # Path to minced.jar
PRODIGAL_PATH = "/path/to/prodigal"             # Path to prodigal executable
HMMSCAN_PATH = "/path/to/hmmscan"               # Path to hmmscan executable

# Model paths - these should point to your downloaded models
CNN_ATT_MODEL_PATH = "data/models/cnn_att_large.pth"
TEMC_CAS_MODEL_PATH = "data/models/temc_cas_best_model.pth"

# Data paths - these should point to your downloaded reference data
HMM_DB_PATH = "data/Profiles"
ADAPTATION_RULES_PATH = "data/adaptation.json" 
INTERFERENCE_RULES_PATH = "data/interference.json"
```

**Note**: The default configuration uses relative paths assuming you follow the directory structure above.

## Usage

### Environment Validation

Before running the pipeline, validate your environment configuration:

```bash
python main.py --validate
```

This will check if all required tools, models, and data files are properly configured.

### Basic Usage

Run DeepCRISPR-Typer on a genome FASTA file:

```bash
python main.py -i input_genome.fasta
```

### Custom Output

Specify a custom output prefix:

```bash
python main.py -i input_genome.fasta -o my_analysis
```

### Example Workflow

1. Place your input FASTA file in `data/input/` or specify the full path
2. Run the validation command to ensure everything is set up correctly
3. Execute the main pipeline with your input file
4. Results will be saved as CSV files in `data/output/<output_prefix>/`

## Output Format

The pipeline generates CSV files containing:

- Array ID and genomic coordinates
- Repeat sequences and predicted subtypes
- Cas protein predictions and families
- HMM scan results and domain matches
- Final classification with confidence scores
- Adaptation and interference scores

## Troubleshooting

### Common Issues

1. **"Missing required paths" error**
   - Ensure all external tools are installed and paths are correctly configured
   - Verify that pre-trained models and reference data are in the correct locations

2. **MinCED not found**
   - Make sure Java is installed and `minced.jar` path is correct
   - Test MinCED independently: `java -jar minced.jar -h`

3. **Model loading errors**
   - Verify that PyTorch version is compatible (>=1.12.0)
   - Ensure model files are not corrupted

4. **Memory issues with large genomes**
   - The pipeline includes parallel processing optimizations
   - For very large datasets, consider splitting the input FASTA file

### Validation Steps

If you encounter issues, run the validation command first:

```bash
python main.py --validate
```

This will provide specific information about missing components.

## Citation

If you use DeepCRISPR-Typer in your research, please cite the relevant publications for:

1. **CNN-Attention Model**: X. Liao; Y. Li; Y. Wu; X. Li; X. Shang. Deep Learning-Based Classification of CRISPR Loci Using Repeat Sequences. ACS Synth. Biol. 2025, 14 (5), 1813–1828. 
2. **TEMC-Cas Model**: X. Liao; Y. Li; Y. Wu; L. Wen; M. Jing; B. Chen; X. Li; X. Shang. TEMC-Cas: Accurate Cas Protein Classification via Combined Contrastive Learning and Protein Language Models. ACS Synth. Biol. 2025, 14 (11), 4586–4596.   
3. **CRISPRCasTyper**: Russel, J.; Pinilla-Redondo, R.; Mayo-Muñoz, D.; Shah, S. A.; Sørensen, S. J. CRISPRCasTyper: automated identification, annotation, and classification of CRISPR-Cas loci. CRISPR J. 2020, 3, 462– 469,  DOI: 10.1089/crispr.2020.0059
4. **ESM-2**: Zeming Lin, Halil Akin, Roshan Rao, et al. Language models of protein sequences at the scale of evolution enable accurate structure prediction. bioRxiv, 2022.

---

**Note**: This tool integrates multiple third-party components. Please ensure you comply with their respective licenses and terms of use.