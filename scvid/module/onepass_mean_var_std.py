# Copyright Contributors to the Cellarium project.
# SPDX-License-Identifier: BSD-3-Clause

import torch
import torch.nn as nn

from scvid.data.util import get_rank_and_num_replicas

from .base_module import BaseModule
from .gather import GatherLayer


class OnePassMeanVarStd(BaseModule):
    """
    Calculate the mean, variance, and standard deviation of the data in one pass (epoch)
    using running sums and running squared sums.

    Args:
        g_genes: Number of genes.
        transform: If not ``None`` is used to transform the input data.
    """

    def __init__(self, g_genes: int, transform: nn.Module | None = None) -> None:
        super().__init__()
        self.g_genes = g_genes
        self.transform = transform
        self.x_sums: torch.Tensor
        self.x_squared_sums: torch.Tensor
        self.x_size: torch.Tensor
        self.register_buffer("x_sums", torch.zeros(g_genes))
        self.register_buffer("x_squared_sums", torch.zeros(g_genes))
        self.register_buffer("x_size", torch.tensor(0))
        self._dummy_param = torch.nn.Parameter(torch.tensor(0.0))

    @staticmethod
    def _get_fn_args_from_batch(
        tensor_dict: dict[str, torch.Tensor]
    ) -> tuple[tuple, dict]:
        x = tensor_dict["X"]
        return (x,), {}

    def forward(self, x_ng: torch.Tensor) -> None:
        if self.transform is not None:
            x_ng = self.transform(x_ng)
        _, num_replicas = get_rank_and_num_replicas()
        if num_replicas > 1:
            x_ng = torch.cat(GatherLayer.apply(x_ng), dim=0)
        self.x_sums = self.x_sums + x_ng.sum(dim=0)
        self.x_squared_sums = self.x_squared_sums + (x_ng**2).sum(dim=0)
        self.x_size = self.x_size + x_ng.shape[0]

    @property
    def mean_g(self) -> torch.Tensor:
        return self.x_sums / self.x_size

    @property
    def var_g(self) -> torch.Tensor:
        return self.x_squared_sums / self.x_size - self.mean_g**2

    @property
    def std_g(self) -> torch.Tensor:
        return torch.sqrt(self.var_g)
