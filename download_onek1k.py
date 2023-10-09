import cellxgene_census

CELLXGENE_ID = "3faad104-2ab8-4434-816d-474d8d2641db"

cellxgene_census.download_source_h5ad(
    CELLXGENE_ID, to_path="../data/onek1k.h5ad"
)
