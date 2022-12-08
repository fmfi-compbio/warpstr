import os

from pandas import read_csv

from src.input_handler.fast5 import normalize_signal_mad
from src.input_handler.input import main_config


class PoreModel():

    def __init__(self, pore_model_path: str) -> None:
        """
        Reads config info from the given yaml config file
        :param pore_model_path: str - path to pore model file
        :return: df - dataframe with kmers and its expected values
        :return: int - size of k-mer in pore model table
        """
        print('nacitavam poremodel')
        if not os.path.exists(pore_model_path):
            raise FileNotFoundError(
                'Not found pore model table at path', pore_model_path)

        pore_model = read_csv(pore_model_path, sep='\t', header=0)
        if 'kmer' not in pore_model and 'level_norm' not in pore_model:
            raise ValueError(
                'Pore model table do not contains "kmer" and "level_mean" columns')

        self.kmersize = len(pore_model['kmer'][0])
        pore_model['level_norm'] = normalize_signal_mad(pore_model['level_mean'])
        self.pore_model = pore_model[['kmer', 'level_norm']]

    def get_value(self, kmer: str) -> float:
        """ Get state value of a kmer from pore model table """
        return self.pore_model[self.pore_model['kmer'] == kmer]['level_norm'].values[0]


pore_model = PoreModel(main_config.pore_model_path)
