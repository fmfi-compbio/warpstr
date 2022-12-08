import h5py
import numpy as np
from scipy.signal import medfilt

from src.input_handler.input import caller_config


class Fast5:
    """
    Class for fast5 files and methods for processing them
    """
    block_size: int
    strand_start: int
    fasta: str
    moves: np.ndarray

    def __init__(self, fast5path: str, get_data_only: bool = False):
        self.handle = h5py.File(fast5path, 'r')
        self.data = None
        self.norm = None
        self.tag = None
        if not (get_data_only):
            self.use_basecall()

    def use_basecall(self):
        for j in [i for i in self.handle['Analyses'].keys() if i.startswith('Basecall')]:
            if 'BaseCalled_template' in self.handle['Analyses'][j].keys():
                self.tag = j[-3:]
        self.bc_tag = 'Basecall_1D_'+self.tag
        self.seg_tag = 'Segmentation_'+self.tag

    def get_tr_extract_reqs(self):
        """
        Getting required information of fast5 for the extraction  of TR signal
        """
        self.fasta = self.load_fasta()
        self.moves = np.asarray(
            self.handle['Analyses'][self.bc_tag]['BaseCalled_template']['Move'])
        self.block_size = int(self.handle['Analyses'][self.bc_tag]['Summary']
                              ['basecall_1d_template'].attrs['block_stride'])
        self.strand_start = int(self.handle['Analyses'][self.seg_tag]['Summary']
                                ['segmentation'].attrs['first_sample_template'])

    def get_data_processed(self, position=None):
        """
        Gets the data from fast5 to use for analysis
        """
        if self.data is None:
            rname = list(self.handle['Raw']['Reads'].keys())[0]
            self.data = np.asarray(
                self.handle['Raw']['Reads'][rname]['Signal'])
            self.remove_spikes()
            self.norm = normalize_signal_mad(self.data)
        if position is not None:
            return self.norm[position[0]:position[1]+1]
        return self.norm

    def get_basecalled_seq(self, position=None):
        """
        Gets the basecalled sequence or subsequence from fast5
        """
        self.fasta = self.load_fasta()
        if position is not None:
            return self.fasta[position[0]:position[1]]
        return self.fasta

    def remove_spikes(self):
        """
        Removes outliers in the fast5 data
        """
        if caller_config.spike_removal == 'median3':
            self.data = medfilt(self.data, 3)
        elif caller_config.spike_removal == 'median5':
            self.data = medfilt(self.data, 5)
        elif caller_config.spike_removal == 'Brute':
            self.data = brute_remove(self.data)

    def load_fasta(self):
        """
        Loads fasta sequence stored in the fast5 data
        """
        if self.fasta is None:
            try:
                fasta = self.handle['Analyses'][self.bc_tag]['BaseCalled_template']['Fastq'][()].decode('ascii').split('\n')[
                    1]
            except KeyError as e:
                print(f'Could not acquire FASTA sequence from .fast5 file due to error={e}')
            else:
                return fasta
        return self.fasta


def brute_remove(data):
    """
    Removes outliers in raw signal data using constants
    :input data: numpy array with raw signal data containing outliers
    :output out: raw signal data with outliers replaced by median of neighbours
    """
    out = data.copy()
    thresh = np.where((data > 1000) | (data < 250))[0]
    for i in thresh:
        if i > 2:
            out[i] = np.median(out[i-2:i+3])
    return out


def normalize_signal_mad(data):
    """
    Normalizes input raw data using MAD
    :input data: numpy array with raw signal data to be normalized
    :output norm_data: normalized numpy array
    """
    robust_quants = (46.5, 53.5)
    shift = np.mean(np.percentile(data, robust_quants))
    scale = np.median(np.abs(data - shift))
    norm_data = (data - shift) / scale
    return np.asarray(norm_data)
