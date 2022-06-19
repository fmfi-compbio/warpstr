import os
import numpy as np
from matplotlib import pyplot as plt

import src.templates as tmpl
from .dna_sequence import get_sequences
from .pore_model import load_pore_model


class Squiggler:
    """
    This class is for generators of expected signals
    """

    def __init__(self, pore_model_path, reference_path):
        self.pore_model, self.kmersize = load_pore_model(pore_model_path)
        self.ref_path = reference_path

    def generate_signal(self, sequence):
        """
        Generates expected signal for input sequence
        :param sequence: string for which to generate signal values
        :return gen_sig: numpy array of normalized expected signal values
        """
        gen_sig = []
        for i in range(len(sequence)-self.kmersize+1):
            kmer = sequence[i:i+self.kmersize]
            gen_sig.append(
                self.pore_model[self.pore_model['kmer'] == kmer]['level_norm'].values[0])

        return np.array(gen_sig)

    def process_locus(self, locus_path, flank_length, coord):
        """
        Process locus
        """
        output_path = os.path.join(locus_path, tmpl.LOCUS_INFO_SUBDIR)

        (lf_t, rf_t, lf_r, rf_r, tmpp, revp) = get_sequences(
            coord, self.ref_path, flank_length)

        self.save_all(output_path, (lf_t, rf_t, lf_r, rf_r, tmpp, revp))
        self.save_gen_signals(
            output_path, (lf_t, rf_t, lf_r, rf_r, tmpp, revp))

    def save_gen_signals(self, output_path, sequences):
        """
        Store expected signals
        """
        gen_signals = []
        for idx, i in enumerate(tmpl.LOCUS_NAMES):
            result = self.generate_signal(sequences[idx])
            saving_path = os.path.join(output_path, i+".txt")
            np.savetxt(saving_path, result, fmt='%f')
            gen_signals.append(result)

        self.save_all_images(gen_signals, tmpl.LOCUS_INFO_SUBDIR, output_path)

    def save_all_images(self, data, name, upper_path):
        """
        Store all Squiggler generated images
        """
        out_path = os.path.join(upper_path, name+".png")
        total_images = len(data)
        fig, axs = plt.subplots(total_images, figsize=(
            16, 3*total_images), sharey=True)
        for idx, i in enumerate(data):
            axs[idx].plot(i, '-o', label='expected signal')
            axs[idx].set_title(tmpl.LOCUS_NAMES[idx])
            axs[idx].set_ylabel('Normalized current')
            axs[idx].legend()
        plt.savefig(out_path, bbox_inches='tight')
        plt.close()

    def save_all(self, output_path, sequences):
        """
        Store all Squiggler used sequences
        """
        fasta_path = os.path.join(output_path, tmpl.LOCUS_FLANKS)
        with open(fasta_path, "w") as file:
            file.write("type,sequence\n")
            for idx, i in enumerate(tmpl.LOCUS_NAMES):
                file.write(i+","+sequences[idx]+"\n")
