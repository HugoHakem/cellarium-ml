import warnings
from typing import Any, Optional

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

from .utils_drug import check_missing_pairs, normalize_drug_fields, standardize_drug_dose, validate_all_doses
from .utils_perturbation import assign_perturbation_id, infer_perturbation_type

# === SECTION: Custom Type ===

class MetadataDict(dict):
    REQUIRED_KEYS = {
        "tissue": "cell line",
        "cell_type": "lymphoblast",
        "cell_line": "K562",
        "disease": "chronic myelogenous leukemia",
        "assay": "10x 3' v1",
        "organism": "Homo sapiens",
        "gene_perturbed": None,
        "gene_perturbation_type": None, # CRISPRa, CRISPRi, CRISPRko
        "drug": None,
        "drug_canonical_smiles": None,
        "drug_dose": None,
        "drug_dose_unit": None,
        "drug_targets": None,
        "drug_moa_broad": None,
        "drug_moa_fine": None,
        "drug_pubchem_cid": None
    }

    _TYPES = {
        "tissue": Optional[str],
        "cell_type": Optional[str],
        "cell_line": Optional[str],
        "disease": Optional[str],
        "assay": Optional[str],
        "organism": Optional[str],
        "gene_perturbed": Optional[str],
        "gene_perturbation_type": Optional[str],
        "drug": Optional[str],
        "drug_canonical_smiles": Optional[str],
        "drug_dose": Optional[float],
        "drug_dose_unit": Optional[str],
        "drug_targets": Optional[str],
        "drug_moa_broad": Optional[str],
        "drug_moa_fine": Optional[str],
        "drug_pubchem_cid": Optional[str],

    }

    def __init__(self, **kwargs):
        # Validate keys and types at initialization
        for key in kwargs:
            if key not in self._TYPES:
                raise KeyError(f"Unexpected metadata key: '{key}'")
            self._validate_type(key, kwargs[key])
        # Merge defaults and user values
        data = self.REQUIRED_KEYS.copy()
        data.update(kwargs)
        super().__init__(data)

    def __setitem__(self, key, value):
        if key not in self._TYPES:
            raise KeyError(f"Unexpected metadata key: '{key}'")
        self._validate_type(key, value)
        super().__setitem__(key, value)

    def __delitem__(self, key):
        if key in self._TYPES:
            raise KeyError(f"Cannot delete required metadata key: '{key}'")
        super().__delitem__(key)

    def _validate_type(self, key, value):
        expected = self._TYPES[key]
        if expected is Optional[float]:
            if value is not None and not isinstance(value, float):
                raise TypeError(f"Key '{key}' must be a float or None, got {type(value).__name__}")
        elif expected is Optional[str]:
            if value is not None and not isinstance(value, str):
                raise TypeError(f"Key '{key}' must be a str or None, got {type(value).__name__}")
        else:
            raise NotImplementedError(f"Unsupported type specification for key '{key}'")

    def __repr__(self):
        return f"MetadataDict({super().__repr__()})"


# === SECTION: Format Obs ===

OBS_DTYPE_MAP = {
    'cell_barcode': str,                 # Unique identifier per cell

    # Perturbation metadata
    'perturbation_type': 'category',     # e.g., "CRISPRi", "drug", "CRISPRa+drug"
    'gene_perturbed': 'category',        # Name of gene, can be NaN for drug-only
    'perturbation_id': 'category',       # Unified ID string, e.g. "TP53_CRISPRi+bortezomib_1(uM)"
    'is_control': bool,                  # True for negative controls (e.g., DMSO or NT guides)

    # Drug metadata
    'drug': 'category',                  # Drug name (e.g., "bortezomib")
    'drug_dose': [str, 'category'],      # List of tuples: [("drug_name", float_dose, "unit")] (Must be cast to str)
    'drug_canonical_smiles': 'category', # Molecular structure (SMILES)
    'drug_targets': 'category',          # List[str] or ";"-delimited string of targets
    'drug_moa_broad': 'category',        # e.g., "Proteasome inhibitor"
    'drug_moa_fine': 'category',         # List[str] or string (if hierarchical MOA)
    'drug_pubchem_cid': 'Int64',         # Nullable integer type

    # Sample/biological metadata
    'tissue': 'category',                # Tissue of origin
    'cell_type': 'category',             # Broad annotation (e.g., "epithelial")
    'cell_line': 'category',             # Specific model (e.g., "MCF7")
    'disease': 'category',               # Disease label (e.g., "breast cancer")
    'organism': 'category',              # e.g., "human"
    'assay': 'category'                  # e.g., "scRNA-seq", "Cell Painting"
}

def cast_obs_types(obs: pd.DataFrame, obs_dtype_map=dict):
    def _cast(ser: pd.Series, dt: Any):
        try:
            if dt == "category":
                # Ensure str conversion doesn't turn NaN into "nan"
                # ser = ser.astype(object)  # avoids auto-casting
                ser = ser.where(ser.notna())  # keep missing values
                ser = ser.astype("category")
            elif dt is str:
                # Avoid casting None to "None" and turn resulting 'nan' to NaN
                ser = ser.where(ser.notna())
                ser = ser.astype(str)
                ser = ser.where(~(ser == "nan"))
            else:
                ser = ser.astype(dt)
        except Exception as e:
            warnings.warn(f"Warning: Could not cast column {ser.name} to {dtype}: {e}", UserWarning)
        return ser

    for col, dtype in obs_dtype_map.items():
        if col in obs:
            if isinstance(dtype, list):
                for d in dtype:
                    obs[col] = _cast(obs[col], d)
            else:
                obs[col] = _cast(obs[col], dtype)


def format_obs(
        adata: ad.AnnData,
        default_metadata: MetadataDict,
        extra_dtype_map: Optional[dict]=None
    ) -> None:
    """
    Extends the .obs DataFrame of an AnnData object with standardized metadata fields in place.

    This function:
      - Adds or updates standard metadata columns (e.g., tissue, cell_type, drug, drug_dose, etc.) in adata.obs.
      - Normalizes drug and drug_dose fields, checks for missing drug/dose/SMILES pairs, and standardizes dose format.
      - Infers and updates the 'perturbation_type' column.
      - Assigns a unique 'perturbation_id' for each row.
      - Adds an 'is_control' boolean column indicating control cells.
      - Adds a 'cell_barcode' column if not present, using the index.

    Args:
        adata (anndata.AnnData): AnnData object whose .obs will be updated.
        default_metadata (MetadataDict): Metadata dictionary with standard fields to add or update.
        extra_dtype_map (dict): additional map of obs's `column` to `dtype` to cast.

    Returns:
        None: The function modifies adata.obs in place.
    """
    default_metadata = MetadataDict(**default_metadata)
    obs = adata.obs.copy()

    for key, val in default_metadata.items():
        if key not in obs.columns and key not in ["gene_perturbation_type", "drug_dose_unit"]:
            obs[key] = val

    drug_mask = ~obs["drug"].isna()
    if drug_mask.sum():

        # check for uncomplete tuple of drug, drug_dose, drug_canonical_smiles
        check_missing_pairs(obs)

        # canonicalize drug_canonical_smiles, split drug_dose at every `+`, turn drug_dose in lists if not already
        obs = normalize_drug_fields(obs)

        # check if every present drug_dose can be standardized
        validate_all_doses(obs)

        # standardize drug_dose
        obs.loc[drug_mask, "drug_dose"] = standardize_drug_dose(
            obs,
            unit=default_metadata["drug_dose_unit"]
        )

        # when normalizing with drug fields, drug is split. Rejoin drugs with `+`
        obs.loc[drug_mask, "drug"] = obs.loc[drug_mask, "drug"].apply(
            lambda x: "+".join(x)
        )

    infer_perturbation_type(
        obs,
        gene_pert_type=default_metadata["gene_perturbation_type"])
    assign_perturbation_id(obs, inplace=True)

    obs["is_control"] = obs["perturbation_id"] == "control"

    # Save original barcode
    if "cell_barcode" not in obs.columns:
        obs = obs.reset_index().rename(columns={"index": "cell_barcode"})

    # Cast obs column into new dtype
    cast_obs_types(obs, OBS_DTYPE_MAP)
    if extra_dtype_map:
        cast_obs_types(obs, extra_dtype_map)

    adata.obs = obs


def format_var(
    adata: ad.AnnData,
    gene_name: pd.DataFrame,
    key_symbol: str="symbol"
) -> ad.AnnData:
    """
    Extends the .var DataFrame of an AnnData object with standardized gene metadata.

    Adds or updates the following columns:
        - UMI_count: Total UMI counts per gene.
        - measured: Boolean indicating if the gene was measured in the experiment.
        - symbol: columns from the provided gene_name DataFrame, joined on Ensembl ID.

    Args:
        adata (anndata.AnnData): AnnData object to extend.
        gene_name (pd.DataFrame): DataFrame with gene metadata, indexed by Ensembl ID.
        key_symbol (str): the expected key to retrieve symbol in gene_name

    Returns:
        anndata.AnnData: AnnData object with extended .var DataFrame.
    """

    adata.var["UMI_count"] = pd.Series(
        np.asarray(adata.X.sum(axis=0)).ravel(),
        index=adata.var_names
    )
    adata.var["measured"] = True

    missing_genes = gene_name.index.difference(adata.var_names)
    X_missing = sp.csr_matrix(
        (
            adata.n_obs,
            len(missing_genes)
        )
    )

    _adata = ad.AnnData(
        X = sp.hstack([adata.X, X_missing], format="csr"),
        obs = adata.obs,
        var = pd.concat(
            [adata.var,
            pd.DataFrame({"UMI_count": 0, "measured": False}, index=missing_genes)],
            axis=0
        )
    )
    adata = _adata[:,
        (_adata.var["UMI_count"] > 0) |
        _adata.var_names.isin(gene_name.index)
    ].copy()
    adata.var = (
        adata.var
        .reset_index()
        .rename(columns={"index": "ensembl_id"})
        .join(gene_name, on="ensembl_id")
        .set_index("ensembl_id")
    )
    adata.var[key_symbol] = adata.var[key_symbol].where(
        ~adata.var[key_symbol].isna(),
        adata.var_names
    )
    return adata
