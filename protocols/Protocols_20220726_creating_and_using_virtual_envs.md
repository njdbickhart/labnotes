# Installing and using virtual environments
---

## General usage

Virtual environments are a way to completely customize your $PATH variable without needing to rewrite the variable each time. By using a virtual environment, you will use a combination of your current $PATH and a new path that leads to a subdirectory. There are several steps for using virtual environments:

1. Create the environment (only has to be done once)
2. Activate the environment (swaps your $PATH)
3. Run your intended program
4. Deactivate the environment (returns your $PATH to its original configuration)

## Specific installations of virtual environments

#### Installation of required modules for Phase Genomics Hic_qc

```bash
# Must be done only once! Not needed as it is already created in your folder
conda create -p /home/tim.smith2/environments/hic_qc python=3.6 numpy pysam matplotlib pdfkit markdown scipy

### Essential steps
# Activate the environment
conda activate /home/tim.smith2/environments/hic_qc

# Run the script
python ~/phase_genomics/hic_qc/hic_qc.py

# Deactivate the environment
conda deactivate
```

#### Installation and use of Flye

```bash
# Must be done only once! Not needed as it is already created in your folder
conda create -p /home/tim.smith2/environments/flye flye

### Essential steps
# Activate the environment
conda activate /home/tim.smith2/environments/flye

# Run the script
flye -h

# Deactivate the environment
conda deactivate
```