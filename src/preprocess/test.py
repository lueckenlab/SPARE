import numpy as np
import pandas as pd

import pytest
import scanpy as sc
from subprocess import run, check_output
import os
import sys

## VIASH START
meta = {
    'executable': '../../target/preprocess',
    # 'executable': 'target/preprocess',
    'resources_dir': '../../data/'
}
## VIASH END


import sys
sys.path.append(meta['resources_dir'])

input_file = f"{meta['resources_dir']}/synthetic.h5ad"
output_file = f"{meta['resources_dir']}/synthetic_processed.h5ad"

print(f"meta['resources_dir'] {meta['resources_dir']}")
print(f"meta['executable'] { meta['executable']}")

def create_synthetic_data():
    n_genes = 100  # Number of genes
    n_cells_per_patient = 10  # number of cells per patient
    n_patients = 10  # Number of patients
    cell_types = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J'] * n_cells_per_patient

    # Creating a synthetic dataset with distinct clusters
    data = np.zeros((n_cells_per_patient * n_patients, n_genes))
    for i in range(n_patients):
        if i < 3:  # First group of patients (dead)
            data[i * n_cells_per_patient:(i + 1) * n_cells_per_patient, :] = np.random.rand(n_cells_per_patient, n_genes) + i * 0.1
        elif i < 7:  # Second group of patients (sick)
            data[i * n_cells_per_patient:(i + 1) * n_cells_per_patient, :] = np.random.rand(n_cells_per_patient, n_genes) + 1.0 + (i - 3) * 0.1
        else:  # Third group of patients (healthy)
            data[i * n_cells_per_patient:(i + 1) * n_cells_per_patient, :] = np.random.rand(n_cells_per_patient, n_genes) + 2.0 + (i - 7) * 0.1

    obs = pd.DataFrame({
        'patient': np.repeat([f'P{i+1}' for i in range(n_patients)], n_cells_per_patient),
        'cell_type': cell_types,
        'outcome': np.repeat([1]*3 + [2]*4 + [3]*3, n_cells_per_patient)
    })

    adata = sc.AnnData(X=data, obs=obs)
    adata.layers["X_raw_counts"] = adata.X.copy()
    return adata



@pytest.fixture
def synthetic_data():
    print("<<<<<<<<<<<<<<<<<<<<< create data >>>>>>>>>>>>>>>>>>>")
    adata = create_synthetic_data()
    adata.write(input_file)
    return adata

def test_preprocess_script(synthetic_data):
    # Run the preprocess script
    print(">>> Run executable")
    cmd_args = [
        meta["executable"],
        "--input", input_file,
        "--output", output_file,
        "--cell_type_key", "cell_type",
        "--batch_covariates", "patient",
        "--batch_key", "patient",
    ]
    result = run(cmd_args, check=True)

    assert result.returncode == 0
    # Load the processed data
    adata_processed = sc.read(output_file)
    
    # Check the structure and content of the processed data
    assert 'X_pca' in adata_processed.obsm_keys()  
    assert 'X_harmony' in adata_processed.obsm_keys() 
    assert 'X_scVI_patient' in adata_processed.obsm_keys()  
    assert 'X_scANVI_patient' in adata_processed.obsm_keys()  
    assert adata_processed.shape[1] == 100  # Ensure the number of genes remains the same
    print("++++++++++++++++++++++++++++++++++++++++++")


if __name__ == '__main__':
    sys.exit(pytest.main([__file__]))
