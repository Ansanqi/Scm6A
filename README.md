# Scm6A
Scm6A is a machine learning tool written in Python 3.8.10 to predict the level of m6A modification in single cells.

# Prerequisites
Before running the script, ensure you have the following libraries installed:

pandas
numpy
pickle
joblib
argparse
scikit-learn

You can install these libraries using pip:

pip install pandas numpy joblib scikit-learn

# Usage
The script is designed to be run from the command line with specified input and output paths.

Command Line Arguments
-I, --input: The absolute path to the input file.
-O, --output: The absolute path to the directory where the output file will be saved.

Example:
python Scm6A.py -I /path/to/input.csv -O /path/to/output

# Other Details

1. Input File
The input file needs to be single-cell gene expression data, with one gene symbol per row and one cell per column.
The input file can be in various formats such as CSV, XLS, XLSX, TSV, or TXT.

2.Output File
The final output is a CSV file containing the predictions along with relevant genomic information.

