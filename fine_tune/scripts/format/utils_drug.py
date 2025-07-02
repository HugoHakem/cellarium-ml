from typing import Any, Optional

import pandas as pd
from rdkit import Chem


def _canonicalize_smiles(smiles: str) -> str:
    """
    Canonicalize a SMILES string.

    Args:
        smiles (str): The SMILES representation of the molecule.

    Returns:
        str: The canonicalized SMILES string.
    """
    return Chem.CanonSmiles(smiles, useChiral=1)


def normalize_drug_fields(df: pd.DataFrame) -> pd.DataFrame:
    """
    Normalize drug and drug_dose fields to lists
    Canonicalize SMILES
    """
    df["drug_dose"] = df["drug_dose"].astype(object)
    drug_mask = df["drug"].notna()

    # Split drugs and ensure drug_dose is a list
    df.loc[drug_mask, "drug"] = df.loc[drug_mask, "drug"].str.split("+")
    df.loc[drug_mask, "drug_dose"] = df.loc[drug_mask, "drug_dose"].apply(
        lambda x: [x] if not isinstance(x, list) else x
    )
    df.loc[drug_mask, "drug_canonical_smiles"] = df.loc[drug_mask, "drug_canonical_smiles"].apply(
        lambda x: _canonicalize_smiles(x)
    )
    return df


def check_missing_pairs(df: pd.DataFrame) -> None:
    """Ensure that no row has a drug without a dose or SMILES"""
    drug_mask = df["drug"].notna()
    incomplete = df.loc[drug_mask, ["drug", "drug_dose", "drug_canonical_smiles"]].isna().any(axis=1)

    if incomplete.any():
        raise ValueError(
            f"There are {incomplete.sum()} incomplete drug entries:\n\n"
            f"{df.loc[drug_mask].loc[incomplete, ['drug', 'drug_dose', 'drug_canonical_smiles']]}"
        )


def _is_valid_dose_entry(dose: Any, drug: str) -> bool:
    """Check if a single dose entry is valid"""
    if isinstance(dose, float):
        return True
    if isinstance(dose, tuple):
        if len(dose) == 2:
            return isinstance(dose[0], float) and isinstance(dose[1], str)
        if len(dose) == 3:
            return dose[0] == drug and isinstance(dose[1], float) and isinstance(dose[2], str)
    return False


def _validate_dose_row(drugs: list[str], doses: list[Any]) -> bool:
    """Validate a row's drug/dose entries"""
    if len(drugs) != len(doses):
        return False
    return all(_is_valid_dose_entry(dose, drug) for drug, dose in zip(drugs, doses))


def validate_all_doses(df: pd.DataFrame) -> None:
    """Check all drug/dose rows for format issues"""
    drug_mask = df["drug"].notna()
    invalid_rows = df.loc[drug_mask].apply(
        lambda row: not _validate_dose_row(row["drug"], row["drug_dose"]),
        axis=1
    )

    if invalid_rows.any():
        raise ValueError(
            "Some `drug_dose` entries are not in a valid format.\n"
            "You must provide either:\n"
            "- a float (interpreted as dose),\n"
            "- a list of floats,\n"
            "- or a list of tuples (drug, dose, unit) or (dose, unit).\n\n"
            f"Invalid entries:\n\n"
            f"{df.loc[drug_mask].loc[invalid_rows, ['drug', 'drug_dose']]}"
        )


def _format_dose_row(drugs: list[str], doses: list[Any], unit: Optional[str]) -> list[tuple[str, float, str]]:
    """Standardize a row of drug/dose pairs"""
    result = []
    for drug, dose in zip(drugs, doses):
        if isinstance(dose, tuple):
            if len(dose) == 2:
                result.append((drug, dose[0], dose[1]))
            elif len(dose) == 3:
                result.append(dose)
        else:
            if not isinstance(unit, str):
                raise ValueError(
                    f"Invalid `unit`: {unit}. "
                    f"If None, must be specified."
                )
            result.append((drug, dose, unit))
    return result


def standardize_drug_dose(df: pd.DataFrame, unit: Optional[str]) -> pd.Series:
    """Returns a new column with standardized (drug, dose, unit) tuples"""
    drug_mask = df["drug"].notna()
    return df.loc[drug_mask].apply(
        lambda row: _format_dose_row(row["drug"], row["drug_dose"], unit),
        axis=1
    )
