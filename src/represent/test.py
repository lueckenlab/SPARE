import pytest
import scanpy as sc
from subprocess import run
import numpy as np
import sys
import patient_representation as pr

## VIASH START
meta = {
    'executable': '../../target/represent',
    'resources_dir': '../../data/'
}
## VIASH END
sys.path.append(meta['resources_dir'])

input_file = f"{meta['resources_dir']}/synthetic_processed.h5ad"
output_file = f"{meta['resources_dir']}/synthetic_represent.h5ad"
print(input_file)


# @pytest.fixture
def synthetic_processed_data():
    # Assuming the preprocess test has already been run and synthetic_processed.h5ad exists
    input_file = f"{meta['resources_dir']}/synthetic_processed.h5ad"
    return sc.read(input_file)

# def test_represent_script(synthetic_processed_data):
def test_represent_script():
    # Run the represent script
    print(">>> Run executable")
    cmd_args = [
        meta["executable"],
        "--input", input_file,
        "--output", output_file,
        "--cell_type_key", "cell_type",
        "--batch_covariates=patient",
        "--sample_key", "patient"
    ]
    result = run(cmd_args, check=True)

    assert result.returncode == 0


    # Load the represented data
    adata_represent = sc.read(output_file)
    print("adata_represent: ")
    print(adata_represent)

    # Check if the representations were stored correctly
    assert 'pseudobulk_pca_distances' in adata_represent.uns_keys()
    assert 'pseudobulk_pca_UMAP' in adata_represent.uns_keys()
    # assert 'ct_pseudobulk_pca_distances' in adata_represent.uns_keys()
    # assert 'ct_pseudobulk_pca_UMAP' in adata_represent.uns_keys()

    distances_pseudobulk = adata_represent.uns['pseudobulk_pca_distances']
    # distances_ct_pseudobulk = adata_represent.uns['ct_pseudobulk_X_pca_distances']
    assert distances_pseudobulk.shape == (10, 10) 
    # assert distances_ct_pseudobulk.shape == (10, 10) 


synthetic_processed_data()
test_represent_script()

# pytest.main(["-v", "test.py"])