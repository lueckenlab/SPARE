import cellxgene_census

CELLXGENE_ID = "dde06e0f-ab3b-46be-96a2-a8082383c4a1"

cellxgene_census.download_source_h5ad(
    CELLXGENE_ID, to_path="data/onek1k.h5ad"
)
