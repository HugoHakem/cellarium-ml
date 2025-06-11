from typing import Optional

import pandas as pd

VALID_PERT_TYPES = {"CRISPRi", "CRISPRa", "CRISPRko", "drug"}


# === SECTION: Perturbation type ===


def infer_perturbation_type(obs: pd.DataFrame, gene_pert_type: Optional[str]) -> pd.DataFrame:
    """
    Infers and standardizes the 'perturbation_type' column in a DataFrame of cell metadata.

    For each row, determines the perturbation type(s) based on the presence and count of gene and drug perturbations.
    If the 'perturbation_type' is missing or incomplete, it fills in the appropriate type(s) using the provided
    default gene perturbation type and the presence of drugs.

    Args:
        obs (pd.DataFrame): DataFrame containing at least 'gene_perturbed' and 'drug' columns.
        gene_pert_type (Optional[str]): Default perturbation type to use for gene perturbations
            (e.g., "CRISPRi", "CRISPRa", or "CRISPRko"). Must be in VALID_PERT_TYPES.

    Returns:
        pd.DataFrame: DataFrame with the 'perturbation_type' column updated or added.

    Raises:
        ValueError: If the default gene perturbation type is invalid or if the perturbation type specification
            is inconsistent with the number of gene/drug perturbations.
    """
    gene_counts = (
        obs["gene_perturbed"].astype(str).str.count(r"\+").fillna(0) +
        obs["gene_perturbed"].notna().astype(int)
    )

    drug_counts = (
        obs["drug"].astype(str).str.count(r"\+").fillna(0) +
        obs["drug"].notna().astype(int)
    )


    def resolve_pert_type(default, g_count, d_count):
        total_count = g_count + d_count

        if pd.isna(default):
            if g_count > 0 and gene_pert_type not in VALID_PERT_TYPES:
                raise ValueError(
                    f"Invalid `gene_pert_type`: {repr(gene_pert_type)}. "
                    f"Must be one of {VALID_PERT_TYPES} when gene perturbation is needed."
                )
            return "+".join(
                [gene_pert_type] * g_count +
                ["drug"] * d_count
            ) if total_count > 0 else None

        parts = default.split("+")
        if not all(p in VALID_PERT_TYPES for p in parts):
            raise ValueError(f"Invalid perturbation types: {parts}")

        # Count actual gene/drug parts
        gene_parts = [p for p in parts if p != "drug"]
        drug_parts = [p for p in parts if p == "drug"]

        # Check if counts match
        if len(gene_parts) == g_count and len(drug_parts) == d_count:
            return default  # Already valid

        # If gene count is incorrect
        if len(gene_parts) not in {0, g_count}:
            raise ValueError(
                f"Invalid `gene_pert_type`: {repr(gene_pert_type)}. "
                f"Must be one of {VALID_PERT_TYPES} when gene perturbation is needed."
            )

        # Fix missing gene part (if completely missing)
        if len(gene_parts) == 0 and g_count > 0:
            if gene_pert_type not in VALID_PERT_TYPES:
                raise ValueError(
                    "`gene_pert_type` must be in " \
                    f"{VALID_PERT_TYPES} when gene perturbation is needed."
                )
            gene_parts = [gene_pert_type] * g_count

        # Fix missing drug part
        if len(drug_parts) != d_count:
            drug_parts = ["drug"] * d_count

        return "+".join(gene_parts + drug_parts) if gene_parts or drug_parts else None

    obs["perturbation_type"] = [
        resolve_pert_type(default, g, d)
        for default, g, d in zip(obs.get("perturbation_type", [None] * len(obs)), gene_counts, drug_counts)
    ]

    return obs


# === SECTION: Perturbation ID ===


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

