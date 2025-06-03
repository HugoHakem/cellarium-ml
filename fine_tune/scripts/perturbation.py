import pandas as pd


def _build_perturbation_id(row: pd.Series) -> str:
    """
    Generate a unique perturbation identifier string from a row of perturbation metadata.

    This function supports:
      - CRISPR perturbations: requires 'perturbation_type' and 'gene_perturbed'
      - Compound perturbations: requires 'perturbation_type', 'compound', and 'compound_dose'
      - Both CRISPR and compound perturbations in the same row (mixed perturbations)

    It supports multiple perturbations separated by '+'. For CRISPR-based perturbations,
    it combines the gene and perturbation type (e.g., "GENE_CRISPRi"). For compound
    perturbations, it combines the compound and dose (e.g., "compound_dose"). If no valid
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
    comp_perts = iter(row["compound"].split("+") if isinstance(row.get("compound", None), str) else [])
    comp_doses = iter(row["compound_dose"].split("+") if isinstance(row.get("compound_dose", None), str) else [])

    pert_ids = []
    for pert in pert_types:
        if "CRISPR" in pert:
            if "gene_perturbed" not in row or row["gene_perturbed"] is None:
                raise ValueError("Missing 'gene_perturbed' for CRISPR perturbation.")
            gene = next(gene_perts, None)
            if gene and gene.upper() != "NT":
                pert_ids.append(f"{gene}_{pert}")
        elif pert == "compound":
            if (("compound" not in row) or
                ("compound_dose" not in row) or
                (row["compound"] is None)):
                raise ValueError("Missing 'compound' or 'compound_dose' for compound perturbation.")
            comp = next(comp_perts, None)
            dose = next(comp_doses, None)
            if comp and comp.upper() != "DMSO":
                if not dose:
                    raise ValueError(
                        "The compound is not `DMSO` and no dosage has been provided. "
                        f"Compound found: {comp}"
                    )
                pert_ids.append(f"{comp}_{dose}")
        else:
            raise ValueError(
                "`perturbation_type` should be a combination of CRISPRi, CRISPRa, CRISPRko, or compound. "
                f"Found unknown type: {pert} in the full perturbation_type: {'+'.join(pert_types)}"
            )
    return "+".join(pert_ids) if pert_ids else "control"


def assign_perturbation_id(df: pd.DataFrame, inplace: bool = True) -> pd.DataFrame:
    """
    Add a 'perturbation_id' column to a DataFrame based on perturbation metadata.

    Supports:
      - CRISPR perturbations: requires 'perturbation_type' and 'gene_perturbed'
      - Compound perturbations: requires 'perturbation_type', 'compound', and 'compound_dose'
      - Both CRISPR and compound perturbations in the same row

    For each row, generates a unique identifier string using _build_perturbation_id.
    The identifier encodes the type and target of the perturbation(s) in the row.

    Args:
        df (pd.DataFrame): Input DataFrame with columns:
            - 'perturbation_type'
            - 'gene_perturbed' (optional, required for CRISPR)
            - 'compound' (optional, required for compound)
            - 'compound_dose' (optional, required for compound)
        inplace (bool, optional): If True, modifies the input DataFrame in place.
                                  If False, returns a modified copy.

    Returns:
        pd.DataFrame: DataFrame with an added 'perturbation_id' column.

    Raises:
        ValueError: If required columns for perturbation identification are missing or inconsistent.
    """
    _df = df if inplace else df.copy()
    _df["perturbation_id"] = _df.apply(_build_perturbation_id, axis=1)
    return _df
