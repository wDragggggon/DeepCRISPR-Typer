# DeepCRISPR-Typer Setup Guide

This guide provides step-by-step instructions for setting up DeepCRISPR-Typer on your system.

## Step 1: Install External Tools

### MinCED Installation
```bash
# Download MinCED
wget https://github.com/ctSkennerton/minced/releases/download/0.4.2/minced-0.4.2.tar.gz
tar -xzf minced-0.4.2.tar.gz
# The minced.jar file will be in the extracted directory
```

### Prodigal Installation
```bash
# On Ubuntu/Debian
sudo apt-get install prodigal

# On CentOS/RHEL  
sudo yum install prodigal

# Or compile from source
wget http://prodigal.ornl.gov/prodigal.v2.6.3.tar.gz
tar -xzf prodigal.v2.6.3.tar.gz
cd prodigal.v2.6.3
make install
```

### HMMER Installation
```bash
# On Ubuntu/Debian
sudo apt-get install hmmer

# On CentOS/RHEL
sudo yum install hmmer

# Or compile from source
wget http://eddylab.org/software/hmmer/hmmer-3.3.2.tar.gz
tar -xzf hmmer-3.3.2.tar.gz
cd hmmer-3.3.2
./configure
make
sudo make install
```

## Step 2: Python Environment Setup

```bash
# Create a virtual environment (recommended)
python -m venv deepcrispr_env
source deepcrispr_env/bin/activate  # On Windows: deepcrispr_env\Scripts\activate

# Install Python dependencies
pip install -r requirements.txt
```

## Step 3: Download Required Data and Models

### CRISPRCasTyper Data Files
```bash
# Clone CRISPRCasTyper repository to get reference data
git clone https://github.com/Russel88/CRISPRCasTyper.git

# Copy required files to your DeepCRISPR-Typer data directory
cp CRISPRCasTyper/data/adaptation.json data/
cp CRISPRCasTyper/data/interference.json data/
cp -r CRISPRCasTyper/data/Profiles data/
```

### Pre-trained Models
You need to obtain the following pre-trained models from their respective publications:

1. **CNN-Attention Model** (`cnn_att_large.pth`)
   - [check repository](https://github.com/Xingyu-Liao/CRISPRclassify-CNN-Att)
   - Place in: `data/models/cnn_att_large.pth`

2. **CRISPRclassify-CNN-Att Source Code**
   - Download the complete CRISPRclassify-CNN-Att repository
   - The pipeline needs access to the `source` directory containing the CNNClassifier implementation
   - You can either:
     - Add the source directory to your Python PATH, or
     - Place it in a location where the pipeline can dynamically import it

3. **TEMC-Cas Model** (`temc_cas_best_model.pth`) - Optional
   - [check repository](https://github.com/Xingyu-Liao/TEMC-Cas)  
   - Place in: `data/models/temc_cas_best_model.pth`

## Step 4: Configure Paths

Edit `config/paths_config.py` to match your installation:

```python
# Update these paths according to your setup
MINCED_PATH = "/path/to/minced.jar"        # Full path to minced.jar
PRODIGAL_PATH = "prodigal"                 # Usually just "prodigal" if in PATH
HMMSCAN_PATH = "hmmscan"                   # Usually just "hmmscan" if in PATH

# Model paths (relative paths work if you follow the directory structure)
CNN_ATT_MODEL_PATH = "data/models/cnn_att_large.pth"
TEMC_CAS_MODEL_PATH = "data/models/temc_cas_best_model.pth"

# Data paths (relative paths)
HMM_DB_PATH = "data/Profiles"
ADAPTATION_RULES_PATH = "data/adaptation.json"
INTERFERENCE_RULES_PATH = "data/interference.json"
```

## Step 5: Validate Installation

Run the validation command to ensure everything is set up correctly:

```bash
python main.py --validate
```

You should see:
```
✅ All required paths are valid
✅ Environment configuration validation passed!
```

## Step 6: Test with Example Data

Place a test genome FASTA file in `data/input/` or use any FASTA file:

```bash
python main.py -i your_genome.fasta
```

Results will be saved in `data/output/<genome_name>/`.

## Troubleshooting

### Common Path Issues

If you get "command not found" errors for prodigal or hmmscan:

1. Check if they are in your PATH: `which prodigal`
2. If not, provide full paths in `paths_config.py`
3. Ensure executable permissions: `chmod +x /path/to/prodigal`

### Java Issues with MinCED

MinCED requires Java. Ensure Java is installed:

```bash
java -version
```

If Java is not available, install OpenJDK:

```bash
# Ubuntu/Debian
sudo apt-get install openjdk-11-jre

# CentOS/RHEL
sudo yum install java-11-openjdk
```

### Model File Permissions

Ensure model files are readable:

```bash
chmod 644 data/models/*.pth
```

## Directory Structure Verification

Your final directory structure should look like:

```
DeepCRISPR-Typer/
├── config/
├── data/
│   ├── input/          # Your input FASTA files
│   ├── models/         # Pre-trained models (.pth files)
│   ├── output/         # Results will be generated here
│   ├── Profiles/       # HMM database from CRISPRCasTyper
│   ├── adaptation.json # From CRISPRCasTyper
│   └── interference.json # From CRISPRCasTyper
├── src/
├── main.py
├── requirements.txt
├── README.md
└── CasScoring.csv
```

Once everything is set up correctly, you're ready to use DeepCRISPR-Typer for CRISPR-Cas system classification!