# Drug Data Processing Pipeline

![Python](https://img.shields.io/badge/Python-3.9%2B-blue.svg)
![License](https://img.shields.io/badge/License-MIT-green.svg)

This project provides a Python pipeline to process a CSV file of drug names and enrich it with chemical information.  

The pipeline includes:
1. Fetching **SMILES** notation for each drug from the PubChem API  
2. Calculating selected **molecular descriptors** with RDKit  
3. Computing **topological indices** from the molecular graph with NetworkX  
4. Exporting all results into a single CSV file  

---

## ⚙️ Installation

Clone the repository and install dependencies:

```bash
git clone https://github.com/aliberkanbektas/psychiatric-qspr-toolkit.git
cd drug-pipeline

python -m venv venv
source venv/bin/activate   # Linux/Mac
venv\Scripts\activate    # Windows

pip install -r requirements.txt
```

---

## 🚀 Usage

Prepare a CSV file `drugs.csv` with the following format:

```csv
Drug
Aspirin
Paracetamol
Ibuprofen
Lorazepam
```

Run the pipeline:

```bash
python process_drugs.py --input drugs.csv --output final_data.csv
```

Options:
- `--input, -i` : Input CSV file (default: `drugs.csv`)  
- `--output, -o` : Output CSV file (default: `final_data.csv`)  
- `--sleep` : Delay between API requests in seconds (default: `0.25`)  

---

## 📂 Project Structure

```
.
├── process_drugs.py        # Main pipeline script
├── drugs.csv               # Example input file
├── final_data.csv          # Output file (generated)
├── requirements.txt        # Dependencies
├── LICENSE                 # MIT license
└── README.md               # Documentation
```

---

## 🧪 Dependencies

- pandas  
- requests  
- rdkit-pypi  
- networkx  

---

## 📄 License

This project is licensed under the [MIT License](https://github.com/aliberkanbektas/psychiatric-qspr-toolkit/blob/main/LICENSE.txt).

## Citation

If you use this project, please cite it using the following DOI: 
[![DOI](https://zenodo.org/badge/1068629608.svg)](https://doi.org/10.5281/zenodo.17252907)
