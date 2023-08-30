import os
from typing import List, Tuple

import numpy as np
from matplotlib import pyplot as plt

import src.templates as tmpl
from src.schemas import Locus

from .dna_sequence import get_sequences
from .pore_model import pore_model


class Squiggler:
    """ Producer of expected signals from poremodel"""

    def __init__(self, reference_path: str):
        self.ref_path = reference_path

    def _generate_signal(self, sequence: str):
        """
        Generates expected signal for input sequence
        :param sequence: string for which to generate signal values
        :return gen_sig: numpy array of normalized expected signal values
        """
        kmers = [sequence[i:i+pore_model.kmersize] for i in range(len(sequence)-pore_model.kmersize+1)]
        gen_sig = [pore_model.get_value(kmer) for kmer in kmers]
        return np.array(gen_sig, dtype=np.float64)

    def process_locus(self, locus: Locus):
        output_path = os.path.join(locus.path, tmpl.LOCUS_INFO_SUBDIR)

        (lf_t, rf_t, lf_r, rf_r, tmpp, revp) = get_sequences(
            locus.coord, self.ref_path, locus.flank_length)

        self._save_all(output_path, (lf_t, rf_t, lf_r, rf_r, tmpp, revp))
        self._save_gen_signals(
            output_path, (lf_t, rf_t, lf_r, rf_r, tmpp, revp))

    def _save_gen_signals(self, output_path: str, sequences: Tuple[str, str, str, str, str, str]):
        """
        Store expected signals
        """
        gen_signals: List[np.ndarray] = []
        for i, seq in zip(tmpl.LOCUS_NAMES, sequences):
            result = self._generate_signal(seq)
            saving_path = os.path.join(output_path, i+'.txt')
            np.savetxt(saving_path, result, fmt='%f')
            gen_signals.append(result)

        self._save_all_images(gen_signals, output_path)

    def _save_all_images(self, data: List[np.ndarray], upper_path: str):
        """
        Store all Squiggler generated images
        """
        out_path = os.path.join(upper_path, tmpl.LOCUS_INFO_SUBDIR+'.png')
        total_images = len(data)
        _, axs = plt.subplots(total_images, figsize=(
            16, 3*total_images), sharey=True)
        for idx, i in enumerate(data):
            axs[idx].plot(i, '-o', label='expected signal')
            axs[idx].set_title(tmpl.LOCUS_NAMES[idx])
            axs[idx].set_ylabel('Normalized current')
            axs[idx].legend()
        plt.savefig(out_path, bbox_inches='tight')
        plt.close()

    def _save_all(self, output_path: str, sequences: Tuple[str, str, str, str, str, str]):
        """ store all squiggler generated sequences """
        fasta_path = os.path.join(output_path, tmpl.LOCUS_FLANKS)
        with open(fasta_path, 'w') as file:
            file.write('type,sequence\n')
            for idx, i in enumerate(tmpl.LOCUS_NAMES):
                file.write(i+','+sequences[idx].upper()+'\n')
