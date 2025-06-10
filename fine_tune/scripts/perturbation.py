import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp


def _build_perturbation_id(row: pd.Series) -> str:
    """
    Generate a unique perturbation identifier string from a row of perturbation metadata.

    This function supports:
      - CRISPR perturbations: requires 'perturbation_type' and 'gene_perturbed'
      - Compound perturbations: requires 'perturbation_type', 'drug', and 'drug_dose'
        'drug_dose' should be a list of tuple [(drug, integer (the dose itself), unit dose (often uM))]
      - Both CRISPR and drug perturbations in the same row (mixed perturbations)

    It supports multiple perturbations separated by '+'. For CRISPR-based perturbations,
    it combines the gene and perturbation type (e.g., "GENE_CRISPRi"). For drug
    perturbations, it combines the drug and dose (e.g., "drug_dose(dose_unit)"). If no valid
    perturbation is found, returns "control".

    Args:
        row (pd.Series): A row from a DataFrame containing perturbation metadata.

    Returns:
        str: A string representing the unique perturbation identifier, or "control" if no perturbation is found.

    Raises:
        ValueError: If required columns for perturbation identification are missing or inconsistent.
    """
    pert_types = row.get("perturbation_type", None)
    if pert_types is None or not isinstance(pert_types, str) or not pert_types.strip():
        return "control"
    pert_types = pert_types.split("+")

    gene_perts = iter(row["gene_perturbed"].split("+") if isinstance(row.get("gene_perturbed", None), str) else [])
    comp_perts = iter(row["drug"].split("+") if isinstance(row.get("drug", None), str) else [])
    comp_doses = iter(row["drug_dose"] if isinstance(row.get("drug_dose", None), list) else [])

    pert_ids = []
    for pert in pert_types:
        if "CRISPR" in pert:
            if "gene_perturbed" not in row or row["gene_perturbed"] is None:
                raise ValueError("Missing 'gene_perturbed' for CRISPR perturbation.")
            gene = next(gene_perts, None)
            if gene and gene.upper() != "NT":
                pert_ids.append(f"{gene}_{pert}")
        elif pert == "drug":
            if (("drug" not in row) or
                ("drug_dose" not in row) or
                (row["drug"] is None)):
                raise ValueError("Missing 'drug' or 'drug_dose' for drug perturbation.")
            comp = next(comp_perts, None)
            dose = next(comp_doses, None)
            if comp and comp.upper() != "DMSO":
                if not dose:
                    raise ValueError(
                        "The drug is not `DMSO` and no dosage has been provided. "
                        f"Compound found: {comp}"
                    )
                pert_ids.append(f"{comp}_{str(dose[1]) + f'({dose[2]})'}")
        else:
            raise ValueError(
                "`perturbation_type` should be a combination of CRISPRi, CRISPRa, CRISPRko, or drug. "
                f"Found unknown type: {pert} in the full perturbation_type: {'+'.join(pert_types)}"
            )
    return "+".join(pert_ids) if pert_ids else "control"


def assign_perturbation_id(df: pd.DataFrame, inplace: bool = True) -> pd.DataFrame:
    """
    Add a 'perturbation_id' column to a DataFrame based on perturbation metadata.

    Supports:
      - CRISPR perturbations: requires 'perturbation_type' and 'gene_perturbed'
      - Compound perturbations: requires 'perturbation_type', 'drug', and 'drug_dose'
      - Both CRISPR and drug perturbations in the same row

    For each row, generates a unique identifier string using _build_perturbation_id.
    The identifier encodes the type and target of the perturbation(s) in the row.

    Example:
        .. code-block:: python

            import pandas as pd
            from fine_tune.scripts.perturbation import assign_perturbation_id

            df = pd.DataFrame({
                "perturbation_type": ["CRISPRi", "drug", "CRISPRi+drug"],
                "gene_perturbed": ["TP53", None, "BRCA1"],
                "drug": [None, "Bortezomib", "Rapamycin"],
                "drug_dose": [None, [("Bortezomib", 10, "uM")], [("Rapamycin", 5, "uM")]],
            })

            assign_perturbation_id(df, inplace=True)
            print(df["perturbation_id"])
            # Output:
            # 0           TP53_CRISPRi
            # 1       Bortezomib_10(uM)
            # 2    BRCA1_CRISPRi+Rapamycin_5(uM)

    Args:
        df (pd.DataFrame): Input DataFrame with columns:
            - 'perturbation_type'
            - 'gene_perturbed' (optional, required for CRISPR)
            - 'drug' (optional, required for drug)
            - 'drug_dose' (optional, required for drug)
        inplace (bool, optional): If True, modifies the input DataFrame in place.
                                  If False, returns a modified copy.

    Returns:
        pd.DataFrame: DataFrame with an added 'perturbation_id' column.

    Raises:
        ValueError: If required columns for perturbation identification are missing or inconsistent.
    """
    _df = df if inplace else df.copy()
    _df["perturbation_id"] = _df.apply(_build_perturbation_id, axis=1)
    if not inplace:
        return _df


def extend_var(
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
