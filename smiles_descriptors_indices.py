# -*- coding: utf-8 -*-
"""
Drug Data Processing Pipeline.

This script orchestrates a pipeline for processing a CSV file of drug names.
The process includes the following steps:
1. Fetches SMILES (Simplified Molecular Input Line Entry System) notation for
   each drug from the PubChem API.
2. Calculates a specific subset of molecular descriptors using the RDKit library.
3. Computes a specific subset of topological indices from the molecular graph
   representation using the NetworkX library.
4. Combines all generated data into a single, comprehensive CSV output file.

Dependencies:
    - pandas
    - requests
    - rdkit-pypi
    - networkx

Usage:
    python process_drugs.py --input <input_file.csv> --output <output_file.csv>
"""
import argparse
import logging
import os
import time
from math import sqrt
from typing import Dict, List, Tuple, Optional
from urllib.parse import quote

import pandas as pd
import requests

# Check for and import required scientific libraries with user-friendly error messages
try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors
except ImportError as e:
    raise ImportError(
        "RDKit library is required. Please install it via 'pip install rdkit-pypi'."
    ) from e

try:
    import networkx as nx
except ImportError as e:
    raise ImportError(
        "NetworkX library is required. Please install it via 'pip install networkx'."
    ) from e

# --- Constants and Configuration ---

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)

# Configuration for PubChem API requests
PUBCHEM_BASE_URL = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{}/property/SMILES/TXT"
REQUEST_TIMEOUT = 6.0
REQUEST_RETRIES = 2
REQUEST_BACKOFF = 0.5

# Definition of topological indices to be calculated
REQUIRED_TOPOLOGICAL_INDICES = {
    "mM2", "H", "SC", "Mod. Sombor", "Red. Mod. Sombor",
    "BSO1", "GA", "RR", "ABC", "N"
}


# --- Stage 1: SMILES Fetching ---

def get_smiles_from_pubchem(drug_name: str) -> Optional[str]:
    """Queries the PubChem API to get the canonical SMILES for a drug.

    This function includes retry logic with a backoff delay to handle
    transient network issues or API rate limiting.

    Args:
        drug_name (str): The name of the drug to query.

    Returns:
        Optional[str]: The SMILES string if found and valid, otherwise None.
    """
    url_encoded_name = quote(drug_name, safe='')
    url = PUBCHEM_BASE_URL.format(url_encoded_name)

    for attempt in range(REQUEST_RETRIES + 1):
        try:
            response = requests.get(url, timeout=REQUEST_TIMEOUT)
            if response.status_code == 200:
                smiles = response.text.strip()
                # Ensure the response is a plausible SMILES string
                return smiles if smiles and any(c.isalpha() for c in smiles) else None
            elif response.status_code == 404:
                logging.debug("'%s' not found on PubChem (404).", drug_name)
                return None
            else:
                logging.warning(
                    "PubChem returned status %d for '%s'. Retrying...",
                    response.status_code, drug_name
                )
        except requests.RequestException as e:
            logging.warning("Request failed for '%s': %s. Retrying...", drug_name, e)

        if attempt < REQUEST_RETRIES:
            time.sleep(REQUEST_BACKOFF)
    return None


def fetch_smiles_for_drugs(df: pd.DataFrame, sleep_per_req: float) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Fetches SMILES from PubChem for each drug in the input DataFrame.

    Iterates through a DataFrame, queries the PubChem API for each drug,
    and segregates the results into found and missing records.

    Args:
        df (pd.DataFrame): DataFrame with a 'Drug' column.
        sleep_per_req (float): Delay in seconds between API calls to be polite.

    Returns:
        Tuple[pd.DataFrame, pd.DataFrame]: A tuple containing two DataFrames:
            - The first with drugs and their found SMILES.
            - The second with drugs for which SMILES were not found.
    """
    found_records, missing_records = [], []
    total = len(df)
    logging.info("Fetching SMILES for %d drugs from PubChem...", total)

    for i, row in df.iterrows():
        drug_name = str(row["Drug"]).strip()
        logging.info("[%d/%d] Querying: %s", i + 1, total, drug_name)
        smiles = get_smiles_from_pubchem(drug_name)
        time.sleep(sleep_per_req) 

        if smiles:
            found_records.append({"Drug": drug_name, "SMILES": smiles})
        else:
            missing_records.append({"Drug": drug_name})

    df_found = pd.DataFrame(found_records)
    df_missing = pd.DataFrame(missing_records)
    logging.info("SMILES search complete. Found: %d, Missing: %d", len(df_found), len(df_missing))
    return df_found, df_missing


# --- Stage 2: RDKit Descriptor Calculation ---

def compute_rdkit_descriptors(df: pd.DataFrame) -> pd.DataFrame:
    """Calculates a selected set of RDKit molecular descriptors.

    Args:
        df (pd.DataFrame): DataFrame must contain a 'SMILES' column.

    Returns:
        pd.DataFrame: A new DataFrame with added descriptor columns.
    """
    df_copy = df.drop_duplicates(subset=["Drug", "SMILES"]).copy()
    logging.info("Calculating RDKit descriptors for %d unique molecules...", len(df_copy))
    
    descriptor_data = []
    for smiles in df_copy["SMILES"]:
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            logging.warning("RDKit could not parse SMILES: %s. Descriptors will be NaN.", smiles)
            descriptor_data.append({})
            continue
        try:
            descriptors = {
                "MolWt": Descriptors.MolWt(mol),
                "HAC": Descriptors.HeavyAtomCount(mol),
            }
            descriptor_data.append(descriptors)
        except Exception as e:
            logging.error("Descriptor calculation failed for '%s': %s", smiles, e)
            descriptor_data.append({})

    df_descriptors = pd.DataFrame(descriptor_data, index=df_copy.index)
    return pd.concat([df_copy, df_descriptors], axis=1)


# --- Stage 3: Topological Index Calculation ---

def calculate_topological_indices(smiles: str) -> Optional[Dict[str, float]]:
    """Creates a molecular graph and calculates topological indices.

    Args:
        smiles (str): A valid SMILES string.

    Returns:
        Optional[Dict[str, float]]: A dictionary of calculated indices, or None if
                                    the SMILES string is invalid.
    """
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return None

    edges = [(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()) for bond in mol.GetBonds()]
    graph = nx.Graph(edges)
    indices = {key: 0.0 for key in REQUIRED_TOPOLOGICAL_INDICES}

    for u, v in graph.edges():
        du, dv = graph.degree[u], graph.degree[v]
        
        # Avoid division by zero
        if du == 0 or dv == 0:
            continue
            
        # --- Index Calculations ---
        indices["mM2"] += 1.0 / (du * dv)
        indices["H"] += 2.0 / (du + dv)
        indices["SC"] += 1.0 / sqrt(du + dv)
        indices["RR"] += 1.0 / sqrt(du * dv)
        indices["GA"] += (2.0 * sqrt(du * dv)) / (du + dv)
        indices["ABC"] += sqrt((du + dv - 2) / (du * dv)) if (du + dv - 2) > 0 else 0.0
        indices["BSO1"] += sqrt(du**2 + dv**2) / (du * dv)
        indices["Mod. Sombor"] += 1.0 / sqrt(du**2 + dv**2)
        indices["Red. Mod. Sombor"] += 1.0 / sqrt((du - 1)**2 + (dv - 1)**2) if (du > 1 or dv > 1) else 0.0
        indices["N"] += sqrt(du + dv)
        
    return indices


def compute_topological_indices(df: pd.DataFrame) -> pd.DataFrame:
    """Applies topological index calculation to a DataFrame.

    Args:
        df (pd.DataFrame): DataFrame containing a 'SMILES' column.

    Returns:
        pd.DataFrame: A DataFrame with topological indices added, with failed
                      rows removed.
    """
    logging.info("Calculating topological indices for %d molecules...", len(df))
    topo_data = []
    failed_indices = []

    for index, row in df.iterrows():
        indices = calculate_topological_indices(row["SMILES"])
        if indices:
            topo_data.append(indices)
        else:
            logging.warning(
                "Skipping topological calculation for Drug '%s' due to invalid SMILES.",
                row.get("Drug", "N/A")
            )
            failed_indices.append(index)
            
    # Align data by dropping failed rows before concatenating
    df_clean = df.drop(index=failed_indices)
    df_topo = pd.DataFrame(topo_data, index=df_clean.index)
    
    return pd.concat([df_clean, df_topo], axis=1)


# --- Main Execution ---

def prepare_input_dataframe(filepath: str) -> Optional[pd.DataFrame]:
    """Loads and validates the input CSV file."""
    if not os.path.exists(filepath):
        logging.error("Input file not found: '%s'", filepath)
        return None
    
    logging.info("Reading input file: %s", filepath)
    df = pd.read_csv(filepath)

    if "Drug" in df.columns:
        return df
    
    # Attempt to find and rename a likely 'Drug' column
    for col in df.columns:
        if col.lower() in ("drug", "name", "compound"):
            df.rename(columns={col: "Drug"}, inplace=True)
            logging.info("Renamed column '%s' to 'Drug' for processing.", col)
            return df
            
    logging.error("A 'Drug' column is required but was not found.")
    return None


def main():
    """Main function to parse arguments and run the full data pipeline."""
    parser = argparse.ArgumentParser(
        description="A pipeline to fetch SMILES and calculate molecular properties for drugs.",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument(
        "--input", "-i",
        default="drugs.csv",
        help="Path to the input CSV file. Must contain a column with drug names (default: 'Drug')."
    )
    parser.add_argument(
        "--output", "-o",
        default="final_data.csv",
        help="Name of the final output CSV file."
    )
    parser.add_argument(
        "--sleep",
        type=float,
        default=0.25,
        help="Seconds to wait between PubChem API requests to avoid rate-limiting."
    )
    args = parser.parse_args()

    try:
        # Step 1: Load and Prepare Data
        df_input = prepare_input_dataframe(args.input)
        if df_input is None:
            return

        # Step 2: Fetch SMILES
        df_smiles, _ = fetch_smiles_for_drugs(df_input, sleep_per_req=args.sleep)
        if df_smiles.empty:
            logging.warning("No SMILES could be found. Terminating pipeline.")
            return

        # Step 3: Calculate RDKit Descriptors
        df_with_descriptors = compute_rdkit_descriptors(df_smiles)

        # Step 4: Calculate Topological Indices
        df_final = compute_topological_indices(df_with_descriptors)

        # Step 5: Save Output
        df_final.to_csv(args.output, index=False)
        logging.info("✅ Pipeline finished successfully. Results saved to '%s'.", args.output)

    except Exception as e:
        logging.critical("An unhandled error occurred during the pipeline execution: %s", e, exc_info=True)


if __name__ == "__main__":
    main()