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

    print("evaluate:")
    pat_instance = pr.tl.TotalPseudobulk(sample_key="patient", cells_type_key="cell_type")
    pat_instance.prepare_anndata(adata_represent, sample_size_threshold=0, cluster_size_threshold=0)
    metadata = pat_instance._extract_metadata(["outcome","patient"])
    outcome = metadata["outcome"]
    patients = metadata["patient"].astype('category').cat.codes
    # outcome = metadata["outcome"].reset_index(drop=True)
    # patients = metadata["patient"].reset_index(drop=True)

    # outcome = pat_instance._extract_metadata(["outcome"])["outcome"].reset_index(drop=True).astype('category').cat.codes
    # patients = pat_instance._extract_metadata(["patient"])["patient"].reset_index(drop=True).astype('category').cat.codes
    print(outcome)
    print(patients)
    # outcome = adata_represent.obs['outcome'].copy()
    # patients = adata_represent.obs['patient'].copy()

    #distances for eval
    distances = [
        # "random_vec_distances",
        "pseudobulk_pca_distances",
       "pseudobulk_scVI_patient_distances", 
        "wasserstein_scANVI_patient_distances"
    ]
    targets = [outcome, patients]
    tasks = ["classification", "regression"]
    for task in tasks:
        print(f"Starting task: {task}")
        for distance in distances:
            for target in targets:
                # print(f"{distance}: {pr.tl.evaluate_representation(distance, target,num_donors_subset=None, proportion_donors_subset=None, method='knn', n_neighbors=5, task=task)}")
                target_name = "outcome" if target is outcome else "patient"
                try:
                    print(f"dist: {distance},task: {task},target: {target_name}: {pr.tl.evaluate_representation(adata_represent.uns[distance], target,num_donors_subset=None, proportion_donors_subset=None, method='knn', n_neighbors=5, task=task)}")
                except Exception as e:
                    print("failed: ", distance, " ", task, "\n", e)
  
    print("The End!")

synthetic_processed_data()
test_represent_script()

# pytest.main(["-v", "test.py"])
# if __name__ == '__main__':
#     sys.exit(pytest.main([__file__]))