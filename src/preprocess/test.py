import numpy as np
import pandas as pd

import pytest
import scanpy as sc
from subprocess import run
import os

# def create_synthetic_data():
#     n_genes = 100  # Number of genes
#     n_cells = 12  # 3 cells per patient
#     n_patients = 4
#     cell_types = ['A', 'B', 'C'] * n_patients

#     # Creating a synthetic dataset with the specified relationships
#     data = np.random.rand(n_cells, n_genes)
#     data[3:6] = data[0:3]  # Patients 1 and 2 are identical
#     data[9:12] = data[0:3] + 0.1  # Patient 4 is close to patients 1 and 2

#     obs = pd.DataFrame({
#         'patient': np.repeat(['P1', 'P2', 'P3', 'P4'], 3),
#         'cell_type': cell_types
#     })

#     adata = sc.AnnData(X=data, obs=obs)
#     adata.layers["X_raw_counts"] = adata.X.copy()
#     return adata

def create_synthetic_data():
    n_genes = 100  # Number of genes
    n_cells_per_patient = 10  # Increase the number of cells per patient
    n_patients = 10  # Increase the number of patients
    cell_types = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J'] * n_cells_per_patient

    # Creating a synthetic dataset with distinct clusters
    data = np.zeros((n_cells_per_patient * n_patients, n_genes))
    for i in range(n_patients):
        if i < 5:  # First half of patients have similar expression patterns
            data[i * n_cells_per_patient:(i + 1) * n_cells_per_patient, :] = np.random.rand(n_cells_per_patient, n_genes) + i * 0.1
        else:  # Second half of patients have distinct but close expression patterns
            data[i * n_cells_per_patient:(i + 1) * n_cells_per_patient, :] = np.random.rand(n_cells_per_patient, n_genes) + 1.5 + i * 0.1

    obs = pd.DataFrame({
        'patient': np.repeat([f'P{i+1}' for i in range(n_patients)], n_cells_per_patient),
        'cell_type': cell_types
    })

    adata = sc.AnnData(X=data, obs=obs)
    adata.layers["X_raw_counts"] = adata.X.copy()
    return adata

@pytest.fixture
def synthetic_data():
    adata = create_synthetic_data()
    adata.write("synthetic.h5ad")
    return adata

def test_preprocess_script(synthetic_data):
    # Run the preprocess script
    result = run(["viash","run","config.vsh.yaml","-p" ,"native", "--", "--input", "synthetic.h5ad", "--output", "synthetic_processed.h5ad", "--cell_type_key", "cell_type", "--batch_covariates=patient", "--batch_effect_covariate", "patient"], check=True)
    assert result.returncode == 0

    # Load the processed data
    adata_processed = sc.read("synthetic_processed.h5ad")
    
    # Check the structure and content of the processed data
    assert 'X_pca' in adata_processed.obsm_keys()  
    # assert 'X_harmony' in adata_processed.obsm_keys() 
    assert 'X_scVI_patient' in adata_processed.obsm_keys()  
    assert 'X_scANVI_patient' in adata_processed.obsm_keys()  
    assert adata_processed.shape[1] == 100  # Ensure the number of genes remains the same


# Run the test with pytest
pytest.main(["-v", "test.py"])