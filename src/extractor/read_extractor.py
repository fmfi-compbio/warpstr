import os
import sys
from dataclasses import dataclass
from multiprocessing import Pool
from typing import List

import h5py
import pysam

import src.templates as tmpl
from src.config import main_config


@dataclass
class Read:
    name: str
    is_reverse: bool
    sam_position: int


class Extractor:
    """
    This class is for extractors of fast5files mapped to specific genomic coordinates
    """
    fast5s: List[Read]

    def __init__(self, root_path: str, bamfile, coord: str, sample: str, full_sample: str):
        """
        :input root_path: lower level path having .bam files and .fast5 files
        :input bamfile: path to bam file
        :input coord: genomic coordinates of desired locus
        :input sample: string id of sample
        :input full_sample: string id of sample and run
        """

        self.path = root_path
        self.bamfile = bamfile
        self.sample = sample
        self.full_sample = full_sample
        self.coord = coord
        self.fast5s = self.get_fast5_list()
        # find fast5s in this bam, mapped to genomic coordinates

        self.summary = []

    def get_fast5_list(self):
        """
        Returns the list of fast5 readnames mapped to the input genomic coordinates
        """
        # fetch reads mapped to the locus using pysam
        fast5s = []
        samfile = pysam.AlignmentFile(self.bamfile, 'r')
        chrom, start, end = process_coord(self.coord)
        for x in samfile.fetch(chrom, start, end):
            if (x.flag >> 11 & 1) == 0 and (x.flag >> 8 & 1) == 0:
                sam_position = start - x.reference_start
                if x.is_reverse:
                    sam_position = len(x.query_sequence)-sam_position
                fast5s.append(Read(
                    name=x.query_name,
                    is_reverse=x.is_reverse,
                    sam_position=sam_position)
                )
        return fast5s

    def extract_from_multifast5s(self, out_path: str):
        """
        Extract single fast5 files using the list of found fast5 files
        :param out_path: path where to extract single fast5 files
        """

        # append prefix to each found read name in the list
        fast5names_lst = ['read_'+f.name for f in self.fast5s]
        dbg_msg = 'There are {num} in {path}'.format(
            num=len(fast5names_lst), path=self.path)
        handle_msg_dbg(dbg_msg)
        # find all multi fast5 files and extract only single fast5s
        for dirpath, dirnames, filenames in os.walk(self.path):
            for filename in [f for f in filenames if f.endswith('.fast5')]:
                fast5path = os.path.join(dirpath, filename)
                fast5file = h5py.File(fast5path, 'r')

                # check if read exists in this multi fast5 file
                found_lst = []
                for readname in fast5names_lst:
                    if readname in fast5file:
                        new_path = os.path.join(
                            out_path, '{}.fast5'.format(readname[5:]))
                        try:
                            copy_fast5(fast5file[readname], new_path)
                        except Exception as e:
                            err_msg = '{err} when copying {read} from {path}'.format(
                                err=e, read=readname, path=fast5file)
                            handle_msg_err(err_msg)

                        found_lst.append(readname)
                fast5names_lst = [name for name in fast5names_lst if name not in found_lst]

            if len(fast5names_lst) == 0:
                dbg_msg = 'All raw fast5 files have been successfully extracted'
                handle_msg_dbg(dbg_msg)
                break

        if len(fast5names_lst) != 0:
            for i in fast5names_lst:
                rname = i[5:]
                to_delete = None
                for j in self.fast5s:
                    if j.name == rname:
                        to_delete = j
                if to_delete is not None:
                    self.fast5s.remove(to_delete)
            return False

    def get_summary(self):
        """
        Prepare summary about reads extracted from .bam file
        """

        summary = []
        basic_info = ','+self.sample+','+self.full_sample+','
        for x in self.fast5s:
            summary.append(x.name+basic_info+str(x.is_reverse)+','+str(x.sam_position))

        self.summary = summary


def get_bam_list(path: str, coord: str, samples: List[str], outpath: str):
    """
    Creates list for later processing, which contains:
    :param path: path with data for a specific sample
    :param bamfile: path to the .bam file
    :param samples: list of samples for which to find reads
    """

    bamfiles = []

    # for each sample find all .bam files
    for sample in samples:
        for full_sample_name in [i for i in os.listdir(path) if i.startswith(sample)]:
            sample_path = os.path.join(path, full_sample_name)
            for dirpath, _, filenames in os.walk(sample_path):
                for filename in [f for f in filenames if f.endswith('.bam')]:
                    bamfile = os.path.join(dirpath, filename)
                    bamfiles.append(
                        (sample_path, bamfile, coord, sample, full_sample_name, outpath))

    dbg_msg = 'found these bamfiles {bamfiles} in {path}'.format(
        bamfiles=bamfiles, path=path)
    handle_msg_dbg(dbg_msg)
    return bamfiles


def process_bam_file(bam_file):
    """
    Wrapper for extracting single fast5s and getting summary
    :param bam_file: list acquired by get_bam_list() passed to extractor init
    :param out_path: root path where extracted single fast5s will be stored
    """
    ex = Extractor(bam_file[0], bam_file[1],
                   bam_file[2], bam_file[3], bam_file[4])
    out_path = bam_file[5]

    # prepare directory in out path for sample associated with this bam file
    sampleout_path = os.path.join(out_path, bam_file[4])
    if os.path.isdir(sampleout_path) is False:
        os.makedirs(sampleout_path)

    # call single fast5 extraction and get summary
    try:
        ex.extract_from_multifast5s(sampleout_path)
    except Exception as e:
        err_msg = '{err} when extracting single reads from {path}'.format(
            err=e, path=sampleout_path)
        handle_msg_err(err_msg)

    ex.get_summary()
    return ex


def process_coord(coord: str):
    """
    Break down genomic coordinates represented as string into separate coordinates
    :param coord: genomic coordinates
    """
    chromosome = coord.split(':')[0]
    start = int(coord.split(':')[1].split('-')[0].replace(',', ''))
    end = int(coord.split('-')[1].replace(',', ''))
    return chromosome, start, end


def copy_fast5(old_read, new_path):
    """
    Copy fast5 read to new path
    :param old_read: fast5 filehandle to the original read
    :param new_path: path where the original read will be copied
    """

    with h5py.File(new_path, 'w') as new_single_read:
        new_single_read.attrs['file_version'] = 2.0

        for group in old_read:
            if group == 'Raw':
                read_number = old_read['Raw'].attrs['read_number']
                new_single_read.copy(
                    old_read[group], 'Raw/Reads/Read_{}'.format(read_number))
            elif group in ('channel_id', 'context_tags', 'tracking_id'):
                if 'UniqueGlobalKey' not in new_single_read:
                    new_single_read.create_group('UniqueGlobalKey')
                new_single_read.copy(
                    old_read[group], 'UniqueGlobalKey/{}'.format(group))
            else:
                new_single_read.copy(old_read[group], group)


def extract_reads(path: str, coord: str, samples: List[str], out_path: str):
    """
    Wrapper for finding and processing bam files in the path
    :param path: high level path containg .bam files and .fast5 files
    :param coord: genomic coordinates as string, for which reads will be extracted
    :param samples: list of samples for which try to extract single reads
    :param out_path: high level path where to copy single fast5s
    :param threads: the number of threads to use
    """

    # find lower level filepaths containing .bam files
    fast5_out_path = os.path.join(out_path, tmpl.FAST5_SUBDIR)
    x = get_bam_list(path, coord, samples, fast5_out_path)

    # process these filepaths in threads
    dbg_msg = 'Running with {thr}'.format(thr=main_config.threads)
    handle_msg_dbg(dbg_msg)

    bam_extracts = []
    if main_config.threads > 1:
        with Pool(main_config.threads) as p:
            bam_extracts = p.map(process_bam_file, x)
    elif main_config.threads == 1:
        bam_extracts = [process_bam_file(i) for i in x]

    # prepare file for writing summary
    overview_file = os.path.join(out_path, tmpl.OVERVIEW_NAME)
    with open(overview_file, 'a') as f:
        for bam_ex in bam_extracts:
            for read in bam_ex.summary:
                f.write(read+'\n')


def handle_msg_err(err_msg):
    print(tmpl.ERR_MSG.format(lvl='2_EXT', msg=err_msg), file=sys.stderr)


def handle_msg_dbg(dbg_msg):
    print(tmpl.DBG_MSG.format(lvl='2_EXT', msg=dbg_msg), file=sys.stderr)
