# Copyright Contributors to the Cellarium project.
# SPDX-License-Identifier: BSD-3-Clause

"""
Testing utilities
-----------------

This module contains helper functions for testing.
"""

from argparse import Namespace
from collections.abc import Callable
from typing import Any

import lightning.pytorch as pl
import numpy as np
import pandas as pd
import torch
from lightning.pytorch.loggers.logger import Logger
from lightning.pytorch.utilities import rank_zero_only
from pytest import approx
from scipy.stats import linregress


def assert_positive(name: str, number: float) -> None:
    """
    Assert that a number is positive.

    Args:
        name: The name of the number.
        number: The number to check.

    Raises:
        ValueError: If the number is not positive.
    """
    if number <= 0:
        raise ValueError(f"`{name}` must be positive. Got {number}")


def assert_nonnegative(name: str, number: float) -> None:
    """
    Assert that a number is non-negative.

    Args:
        name: The name of the number.
        number: The number to check.

    Raises:
        ValueError: If the number is negative.
    """
    if number < 0:
        raise ValueError(f"`{name}` must be non-negative. Got {number}")


def assert_columns_and_array_lengths_equal(
    matrix_name: str,
    matrix: np.ndarray | torch.Tensor,
    array_name: str,
    array: np.ndarray | torch.Tensor,
) -> None:
    """
    Assert that the number of columns in a matrix matches the length of an array.

    Args:
        matrix_name: The name of the matrix.
        matrix: The matrix.
        array_name: The name of the array.
        array: The array.

    Raises:
        ValueError: If the number of columns in the matrix does not match the length of the array.
    """
    if matrix.shape[1] != len(array):
        raise ValueError(
            f"The number of `{matrix_name}` columns must match the `{array_name}` length. "
            f"Got {matrix.shape[1]} != {len(array)}"
        )


def assert_arrays_equal(
    a1_name: str,
    a1: np.ndarray,
    a2_name: str,
    a2: np.ndarray,
) -> None:
    """
    Assert that two arrays are equal.

    Args:
        a1_name: The name of the first array.
        a1: The first array.
        a2_name: The name of the second array.
        a2: The second array.

    Raises:
        ValueError: If the arrays are not equal.
    """
    if not np.array_equal(a1, a2):
        raise ValueError(f"`{a1_name}` must match `{a2_name}`. Got {a1} != {a2}")


def assert_slope_equals(data: pd.Series, slope: float, loglog: bool = False, atol: float = 1e-4) -> None:
    """
    Assert that the slope of a series is equal to a given value.

    Args:
        data:
            The :class:`pandas.Series` object to check.
        slope:
            Expected slope.
        loglog:
            Whether to use log-log scale.
        atol:
            The absolute tolerance.

    Raises:
        ValueError: If the slope is not equal to the given value.
    """
    x = np.log(data.index) if loglog else data.index
    y = np.log(data) if loglog else data
    actual_slope = linregress(x, y).slope
    if actual_slope != approx(slope, abs=atol):
        raise ValueError(f"Slope of data must be {slope}. Got {actual_slope}")


"""
Helper functions for coordinate data collection.
"""


def get_coord_data(
    models: dict[int, Callable[[], pl.LightningModule]],
    layer_name_to_multiplier_name: dict[str, str],
    train_loader: torch.utils.data.DataLoader,
    nsteps: int,
    nseeds: int,
) -> pd.DataFrame:
    """
    Get coordinate data for a model.

    Args:
        models:
            A dictionary mapping width to a function that returns a model.
        layer_name_to_multiplier_name:
            A dictionary mapping layer names to their corresponding multipliers.
        train_loader:
            The training data loader.
        nsteps:
            The number of steps to train for.
        nseeds:
            The number of seeds to use.

    Returns:
        A :class:`pandas.DataFrame` containing the coordinate data.
    """
    from cellarium.ml.callbacks import GetCoordData

    df = pd.DataFrame()
    for i in range(nseeds):
        torch.manual_seed(i)
        for width, lazy_model in models.items():
            logger = PandasLogger()
            trainer = pl.Trainer(
                devices=1,
                max_steps=nsteps,
                logger=logger,
                callbacks=[GetCoordData(layer_name_to_multiplier_name)],
                enable_checkpointing=False,
                enable_model_summary=False,
                enable_progress_bar=False,
            )
            model = lazy_model()
            trainer.fit(model, train_loader)
            width_df = logger.df
            width_df["seed"] = i
            width_df["width"] = width
            df = pd.concat([df, width_df], ignore_index=True)

    return df


class PandasLogger(Logger):
    def __init__(self) -> None:
        self.df = pd.DataFrame()

    @property
    def name(self):
        return "PandasLogger"

    @property
    def version(self):
        return ""

    @rank_zero_only
    def log_metrics(self, metrics: dict[str, float], step: int | None = None) -> None:
        metric_df = pd.DataFrame(metrics, index=[step])
        self.df = pd.concat([self.df, metric_df])

    @rank_zero_only
    def log_hyperparams(self, params: dict[str, Any] | Namespace, *args: Any, **kwargs: Any) -> None:
        pass


def coord_check_MLP(df: pd.DataFrame, nsteps: int, bias: bool, optimizer: str, atol: float) -> None:
    for t in range(nsteps):
        df_t = df.loc[df.t == t]
        for i in range(1, 4):  # layers
            df_i = df_t.loc[df_t.module.str.contains(str(i))]
            # output scaling: Θ(1)
            df_out = df_i.loc[df.type == "out"]
            out = df_out.groupby("width").l1.mean()
            if any(out):
                assert_slope_equals(out, 0, loglog=True, atol=atol)

            # W
            df_param = df_i.loc[df.type == "param"]
            W = df_param[df_param.module.str.contains("weight")].groupby("width").l1.mean()
            if i == 1:
                # input layer scaling: Θ(1)
                assert_slope_equals(W, 0, loglog=True, atol=atol)
            elif i == 2:
                # hidden layer scaling: Θ(1/sqrt(n))
                assert_slope_equals(W, -0.5, loglog=True, atol=atol)
            elif i == 3:
                # output layer scaling: Θ(1/n)
                if any(W):
                    assert_slope_equals(W, -1, loglog=True, atol=atol)

            # b scaling: Θ(1)
            if bias:
                b = df_param[df_param.module.str.contains("bias")].groupby("width").l1.mean()
                if any(b):
                    assert_slope_equals(b, 0, loglog=True, atol=atol)

            # ∆W
            df_delta = df_i.loc[df.type == "delta"]
            dW = df_delta[df_delta.module.str.contains("weight")].groupby("width").l1.mean()
            if i == 1:
                # input layer scaling: Θ(1)
                if any(dW):
                    assert_slope_equals(dW, 0, loglog=True, atol=atol)
            elif i == 2:
                # hidden layer scaling: Θ(1/n)
                if any(dW):
                    if optimizer == "adamw":
                        # additional scaling due to weight decay: Θ(1/sqrt(n))
                        assert_slope_equals(dW, -0.75, loglog=True, atol=0.27)
                    else:
                        assert_slope_equals(dW, -1, loglog=True, atol=atol)
            elif i == 3:
                # output layer scaling: Θ(1/n)
                assert_slope_equals(dW, -1, loglog=True, atol=atol)

            # ∆b scaling: Θ(1)
            if bias:
                db = df_delta[df_delta.module.str.contains("bias")].groupby("width").l1.mean()
                if any(db):
                    assert_slope_equals(db, 0, loglog=True, atol=atol)
