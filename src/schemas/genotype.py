from dataclasses import dataclass, field
from typing import List

import numpy as np


@dataclass
class Genotype:
    group1: List[int]
    group2: List[int] = field(default_factory=list)
    predictions: List[int] = field(default_factory=list)

    @property
    def is_hetero(self) -> bool:
        return len(self.group2) > 0

    @property
    def first_allele(self):
        return int(find_nearest(self.group1, np.median(self.group1)))

    @property
    def second_allele(self):
        return int(find_nearest(self.group2, np.median(self.group2)) if self.is_hetero else '-')

    @property
    def first_allele_sz(self):
        return len(self.group1)

    @property
    def second_allele_sz(self):
        return len(self.group2) if self.is_hetero else '-'

    @property
    def alleles(self):
        return (self.first_allele, self.second_allele)


def find_nearest(array: List[int], value: float):
    return array[(np.abs(np.asarray(array) - value)).argmin()]
