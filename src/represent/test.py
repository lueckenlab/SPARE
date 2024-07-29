import pytest
import scanpy as sc
from subprocess import run

@pytest.fixture
def synthetic_processed_data():
    # Assuming the preprocess test has already been run and synthetic_processed.h5ad exists
    return sc.read("synthetic_processed.h5ad")

def test_represent_script(synthetic_processed_data):
    # Run the represent script
    result = run(["viash","run","config.vsh.yaml","--", "--input", "synthetic_processed.h5ad", "--output", "synthetic_represent.h5ad", "--cell_type_key", "cell_type","--batch_covariates=patient", "--sample_key", "patient"], check=True)
    assert result.returncode == 0

    # Load the represented data
    adata_represent = sc.read("synthetic_represent.h5ad")

    # Check if the representations were stored correctly
    assert 'pseudobulk_pca_distances' in adata_represent.uns_keys()
    assert 'pseudobulk_pca_UMAP' in adata_represent.uns_keys()
    # assert 'ct_pseudobulk_pca_distances' in adata_represent.uns_keys()
    # assert 'ct_pseudobulk_pca_UMAP' in adata_represent.uns_keys()

    distances_pseudobulk = adata_represent.uns['pseudobulk_pca_distances']
    # distances_ct_pseudobulk = adata_represent.uns['ct_pseudobulk_X_pca_distances']
    assert distances_pseudobulk.shape == (10, 10) 
    # assert distances_ct_pseudobulk.shape == (10, 10) 

pytest.main(["-v", "test.py"])
