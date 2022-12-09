import os
import re
from typing import Dict, List, Tuple

import numpy as np
from pandas import read_csv

import src.templates as tmpl
from src.config import main_config
from src.schemas.fast5 import normalize_signal_mad


class PoreModel():

    def __init__(self, pore_model_path: str) -> None:
        """
        Reads config info from the given yaml config file
        :param pore_model_path: str - path to pore model file
        :return: df - dataframe with kmers and its expected values
        :return: int - size of k-mer in pore model table
        """
        if not os.path.exists(pore_model_path):
            raise FileNotFoundError(
                'Not found pore model table at path', pore_model_path)

        pore_model = read_csv(pore_model_path, sep='\t', header=0)
        if 'kmer' not in pore_model and 'level_norm' not in pore_model:
            raise ValueError(
                'Pore model table do not contains "kmer" and "level_mean" columns')

        self.kmersize = len(pore_model['kmer'][0])
        pore_model['level_norm'] = normalize_signal_mad(pore_model['level_mean'])
        self.table = pore_model[['kmer', 'level_norm']]

    def _get_consecutive_diff(self, pattern: str) -> Tuple[float, float]:
        """
        Get mean and median differences between expected signals
        """
        rep = pattern*self.kmersize
        kmers = [rep[i:self.kmersize+i] for i in range(len(pattern)+1)]
        pore_values = np.abs(np.diff(
            [self.get_value(kmer) for kmer in kmers]))
        return np.mean(pore_values), np.median(pore_values)

    def get_value(self, kmer: str) -> float:
        """ Get state value of a kmer from pore model table """
        return self.table[self.table['kmer'] == kmer]['level_norm'].values[0]

    def get_diffs_for_all(self, sequence: str) -> Dict[str, Tuple[float, float]]:
        """
        Gets differences between expected signals for each unit
        """
        diffs: Dict[str, Tuple[float, float]] = {}
        brackets = ['(', ')', '{', '}']
        res = re.findall(r'[\(\{].*?[\)\}]', sequence)

        for r in res:
            base_pattern = ''.join([c for c in r if c not in brackets])
            pattern = ['']
            for char in base_pattern:
                if char not in tmpl.DNA_DICT:
                    for idx, p in enumerate(pattern):
                        pattern[idx] = p + char
                else:
                    new: List[str] = []
                    for iupac in tmpl.DNA_DICT[char]:
                        for p in pattern:
                            new.append(p+iupac)
                    pattern = new
            for p in pattern:
                diffs[p] = pore_model._get_consecutive_diff(p)

        return diffs


pore_model = PoreModel(main_config.pore_model_path)
