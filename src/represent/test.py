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



@pytest.fixture(scope="module")
def represented_data():
    # obtaining representations
    print(">>> Run executable")
    cmd_args = [
        meta["executable"],
        "--input", input_file,
        "--output", output_file,
        "--cell_type_key", "cell_type",
        "--batch_covariates=patient",
        "--sample_key", "patient",
        #setting sample size threshold to 0 for testing the synthetic data to avoid removing all samples
        "--celltype_pseudobulk_sample_size_threshold", "0",
    ]
    result = run(cmd_args, check=True)
    assert result.returncode == 0
    return sc.read(output_file)

@pytest.fixture(scope="module")
def prepare_data(represented_data):
    adata_represent = represented_data

    # Prepare data
    pat_instance = pr.tl.TotalPseudobulk(sample_key="patient", cells_type_key="cell_type")
    pat_instance.prepare_anndata(adata_represent, sample_size_threshold=0, cluster_size_threshold=0)
    metadata = pat_instance._extract_metadata(["outcome", "patient"])
    
    return metadata

def test_representation_keys(represented_data):
    adata_represent = represented_data
    print(adata_represent)

    # Check if the representations were stored correctly
    assert 'pseudobulk_pca_distances' in adata_represent.uns_keys()
    assert 'pseudobulk_pca_UMAP' in adata_represent.uns_keys()
    assert 'ct_pseudobulk_pca_distances' in adata_represent.uns_keys()
    assert 'ct_pseudobulk_pca_UMAP' in adata_represent.uns_keys()

    distances_pseudobulk = adata_represent.uns['pseudobulk_pca_distances']
    distances_ct_pseudobulk = adata_represent.uns['ct_pseudobulk_pca_distances']
    assert distances_pseudobulk.shape == (10, 10)
    assert distances_ct_pseudobulk.shape == (10, 10)

# we expect score 0 for the following, because each patient has unique identifier and can not be classified by other's patient keys
def test_classification_patient(represented_data):
    adata_represent = represented_data

    # Prepare data
    pat_instance = pr.tl.TotalPseudobulk(sample_key="patient", cells_type_key="cell_type")
    pat_instance.prepare_anndata(adata_represent, sample_size_threshold=0, cluster_size_threshold=0)
    metadata = pat_instance._extract_metadata(["outcome", "patient"])
    patients = metadata["patient"].astype('category').cat.codes

    # Perform classification on patient target
    distances = [
        "pseudobulk_pca_distances",
        "pseudobulk_scVI_patient_distances",
        "wasserstein_scANVI_patient_distances"
    ]
    
    for distance in distances:
        score = pr.tl.evaluate_representation(
            adata_represent.uns[distance], patients,
            num_donors_subset=None, proportion_donors_subset=None,
            method='knn', n_neighbors=5, task="classification"
        )['score']
        assert score == 0.0, f"Expected 0 score, but got {score} for {distance} with patient target for the classification."

# def test_represent_script():
#     # Run the represent script
#     print(">>> Run executable")
#     cmd_args = [
#         meta["executable"],
#         "--input", input_file,
#         "--output", output_file,
#         "--cell_type_key", "cell_type",
#         "--batch_covariates=patient",
#         "--sample_key", "patient",
#         #setting sample size threshold to 0 for testing the synthetic data to avoid removing all samples
#         "--celltype_pseudobulk_sample_size_threshold","0",
#     ]
#     result = run(cmd_args, check=True)

#     assert result.returncode == 0


#     # Load the represented data
#     adata_represent = sc.read(output_file)
#     print("adata_represent: ")
#     print(adata_represent)

#     # Check if the representations were stored correctly
#     assert 'pseudobulk_pca_distances' in adata_represent.uns_keys()
#     assert 'pseudobulk_pca_UMAP' in adata_represent.uns_keys()
#     assert 'ct_pseudobulk_pca_distances' in adata_represent.uns_keys()
#     assert 'ct_pseudobulk_pca_UMAP' in adata_represent.uns_keys()

#     distances_pseudobulk = adata_represent.uns['pseudobulk_pca_distances']
#     distances_ct_pseudobulk = adata_represent.uns['ct_pseudobulk_pca_distances']
#     assert distances_pseudobulk.shape == (10, 10) 
#     assert distances_ct_pseudobulk.shape == (10, 10) 

#     print("evaluate:")
#     #### add coments
#     pat_instance = pr.tl.TotalPseudobulk(sample_key="patient", cells_type_key="cell_type")
#     pat_instance.prepare_anndata(adata_represent, sample_size_threshold=0, cluster_size_threshold=0)
#     metadata = pat_instance._extract_metadata(["outcome","patient"])
#     outcome = metadata["outcome"]
#     patients = metadata["patient"].astype('category').cat.codes
 
#     print(outcome)
#     print(patients)

#     #distances for eval
#     distances = [
#         # "random_vec_distances",
#         "pseudobulk_pca_distances",
#        "pseudobulk_scVI_patient_distances", 
#         "wasserstein_scANVI_patient_distances"
#     ]
#     targets = [outcome, patients]
#     tasks = ["classification", "regression"]
#     for task in tasks:
#         print(f"Starting task: {task}")
#         for distance in distances:
#             for target in targets:
#                 target_name = "outcome" if target is outcome else "patient"
#                 try:
#                     print(f"dist: {distance},task: {task},target: {target_name}: {pr.tl.evaluate_representation(adata_represent.uns[distance], target,num_donors_subset=None, proportion_donors_subset=None, method='knn', n_neighbors=5, task=task)}")
#                 except Exception as e:
#                     print("failed: ", distance, " ", task, "\n", e)
  
#     print("The End!")



if __name__ == '__main__':
    sys.exit(pytest.main([__file__]))
