import os
import sys
from dataclasses import dataclass
from multiprocessing import Pool
from typing import List, Optional, Tuple

import numpy as np
import pandas as pd
from Bio import pairwise2

import src.templates as tmpl
from src.config import alignment_config, main_config
from src.schemas import Fast5, Locus
from src.squiggler.dna_sequence import get_reverse_strand


@dataclass
class ReadForExtraction:
    name: str
    run_id: str
    reverse: bool
    approx_location: int


@dataclass
class Mapping:
    ref: str = '-'
    mapping: str = '-'
    query: str = '-'


@dataclass
class Position:
    start: int
    end: int

    @property
    def valid(self) -> bool:
        if self.start > self.end:
            return False
        elif self.start == -1 or self.end == -1:
            return False
        else:
            return True


@dataclass
class Alignment:
    score: int = -1
    identity: float = -1.0
    position: Position = Position(-1, -1)
    mapping: Mapping = Mapping()

    def __post_init__(self):
        if self.score <= alignment_config.accuracy_factor and self.identity <= alignment_config.identity_factor:
            self.position = Position(-1, -1)

    @property
    def score_string(self) -> str:
        return f'Score: {self.score} Identity: {self.identity}'

    @property
    def found(self) -> bool:
        return self.position.valid


@dataclass
class Flank:
    left: str
    right: str


@dataclass
class Flanks:
    template: Flank
    reverse: Flank

    def get_by_reverse(self, reverse: bool) -> Flank:
        if reverse:
            return self.reverse
        return self.template


@dataclass
class FlankInRead:
    read_id: str
    lflank_raw: Position
    rflank_raw: Position
    l_alignment: Alignment
    r_alignment: Alignment
    sequence: Optional[str]

    @property
    def alignment_mappings(self):
        return [
            self.l_alignment.mapping,
            self.r_alignment.mapping
        ]

    @property
    def valid(self):
        if self.l_alignment.found and self.r_alignment.found and self.lflank_raw.valid and self.rflank_raw.valid:
            return 1
        else:
            return 0


def load_flanks(locus_path: str) -> Flanks:
    """
    Function loading template and reverse flanks
    :param locus_path: Path to the extracted data for selected locus
    :returns (lft,rft): left and right flanking template sequences
    :returns (lfr,rfr): left and right flanking reverse sequences
    """
    path = os.path.join(locus_path, tmpl.LOCUS_INFO_SUBDIR, tmpl.LOCUS_FLANKS)
    if not os.path.exists(path):
        raise FileNotFoundError(f'File with flanks not found in path={path}')

    with open(path, 'r') as f:
        colnames = f.readline().split(',')
        seqid = 0
        for idx, i in enumerate(colnames):
            if i.rstrip() == 'sequence':
                seqid = idx

        lft = f.readline().split(',')[seqid].rstrip()
        rft = f.readline().split(',')[seqid].rstrip()
        lfr = f.readline().split(',')[seqid].rstrip()
        rfr = f.readline().split(',')[seqid].rstrip()

    return Flanks(
        template=Flank(
            left=lft,
            right=rft
        ),
        reverse=Flank(
            left=lfr,
            right=rfr
        )
    )


def handle_msg_dbg(dbg_msg):
    print(tmpl.DBG_MSG.format(lvl='3_TREX', msg=dbg_msg), file=sys.stderr)


def transform_moves(moves: np.ndarray) -> np.ndarray:
    """
    Preparation of Move table from guppy for easier indexing.
    Guppy Move table contains info whether current window has the same context or not
    This information is stored as zeros and ones for easy compression
    However for easier indexing we transform it into table of contexts
    :param moves: Move table
    :return moves_r: transformed Moves table
    """

    moves_r = np.zeros(len(moves), dtype=np.int32)
    moves_r[0] = 0
    for idx, i in enumerate(moves[1:], start=1):
        if i:
            moves_r[idx] = 1+moves_r[idx-1]
        else:
            moves_r[idx] = moves_r[idx-1]
    return moves_r


def extract_from_moves(moves_r: np.ndarray, pos: Position, strand_start: int, block_stride: int) -> Position:
    """
    Computing raw signal positions corresponding approximately to the input basecalled nucleotides
    Input positions represent start and end of flanks in basecalled sequence
    This method finds which raw signal values correspond to those positions
    :param moves_r: transformed Moves table for easy indexing
    :param pos1: position of the first context
    :param pos2: position of the second context
    :param strand_start: position in raw signal from which basecalling started
    :param block_stride: guppy parameter saying how many values overlap in block
    :return start: start position where the first context occured
    :return end: end position where the last context occured
    """
    decoded_pos1 = np.where(moves_r == pos.start)[0]
    decoded_pos2 = np.where(moves_r == pos.end)[0]

    if len(decoded_pos1) > 0:
        start = strand_start + decoded_pos1[0]*block_stride
    else:
        start = -1

    if len(decoded_pos2) > 0:
        end = strand_start + decoded_pos2[-1]*block_stride
    else:
        end = -1

    return Position(start, end)


def find_sequence(seq1: str, seq2: str, origin_offset: int = 0) -> Alignment:
    """
    Extraction of tandem repeat region using flank sequences.
    :param seq1: text sequence
    :param seq2: pattern sequence
    :return score: score of the alignment
    :return real_start: start of the alignment in the basecalled sequence
    :return end: calculated alignments of flanks with basecalled read
    :return (score_string,ref,mapping,query): calculated alignments of flanks with basecalled read
    """
    al1: str
    al2: str
    score: int
    start: int
    end: int

    # get the best alignment using local alignment from pairwise lib
    al1, al2, score, start, end = pairwise2.align.localms(
        seq1, seq2, alignment_config.match_score, alignment_config.mismatch_score, alignment_config.gap_open_score,
        alignment_config.gap_extend_score, one_alignment_only=True)[0]

    # decode the start and end positions of the alignment so alignment can be written to file later
    # this also handles cases when alignment starts or ends earlier
    # and ensures that the alignment is full length with the input pattern
    nums_gaps: int = al1[start:end].count('-')
    nums_gaps2: int = al2[start:end].count('-')
    real_start: int = al2[:start].count('-')
    end = real_start+len(seq2)+nums_gaps2-nums_gaps

    # prepare alignment as list of strings
    ref: str = al1[real_start:end]
    query: str = al2[real_start:end]
    mapping = ''
    identity = 0
    for i, j in zip(ref, query):
        if i == j:
            mapping += tmpl.ALIGNMENT_MATCH_CHAR
            identity += 1
        else:
            mapping += tmpl.ALIGNMENT_MISMATCH_CHAR

    # Save also score into the alignment
    diff = (len(seq2)-(len(query)-nums_gaps2))*alignment_config.gap_extend_score
    score = int(score + diff)
    identity = identity/(len(ref))
    return Alignment(
        score=score,
        identity=identity,
        position=Position(real_start+origin_offset, end+origin_offset),
        mapping=Mapping(
            ref=ref,
            mapping=mapping,
            query=query
        )
    )


def align_seq(read: str, flank: Flank) -> Tuple[Alignment, Alignment]:
    """Align both flanks to the read sequence

    Args:
        read (str): read sequence where to find flanks
        flanks (Flank): left and right flank sequence

    Returns:
        Tuple[Alignment, Alignment]: Found alignments of flanks
    """

    left_alignment = find_sequence(read, flank.left)

    # find right alignment only in the part of the sequence after the left alignment
    if (left_alignment.position.end+len(flank.right)) > len(read):
        right_alignment = Alignment()
    else:
        right_alignment = find_sequence(read[left_alignment.position.end:], flank.right, left_alignment.position.end)

    return left_alignment, right_alignment


def extract_tr(overview_row: Tuple[ReadForExtraction, Locus]) -> FlankInRead:
    readname: str = overview_row[0].name
    runid: str = overview_row[0].run_id
    reverse: bool = overview_row[0].reverse
    approx_location: int = overview_row[0].approx_location
    locus: Locus = overview_row[1]

    fast5path = os.path.join(locus.path, tmpl.FAST5_SUBDIR, str(
        runid), tmpl.ANNOT_SUBDIR, readname+'.fast5')

    # decode flanks according to strand
    flanks = load_flanks(locus.path)
    flank = flanks.get_by_reverse(reverse)
    # load fast5 file and guppy info
    try:
        fast5 = Fast5(fast5path)
        fast5.get_tr_extract_reqs()
    except Exception as e:
        raise ValueError('Error={err} when loading fast5 from {path}'.format(err=e, path=fast5path))

    if approx_location > len(fast5.fasta):
        approx_location = len(fast5.fasta)
    elif approx_location < 0:
        approx_location = 0

    start = max(int(approx_location-0.05*len(fast5.fasta)-5000), 0)
    end = min(int(approx_location+0.05*len(fast5.fasta)+5000), len(fast5.fasta)-1)

    ref_text = fast5.fasta[start:end]
    left_alignment, right_alignment = align_seq(ref_text, flank)

    if not left_alignment.found or not right_alignment.found:
        left_raw = Position(-1, -1)
        right_raw = Position(-1, -1)
        seq = None
    else:
        left_alignment.position.start += start
        right_alignment.position.start += start
        left_alignment.position.end += start
        right_alignment.position.end += start

        # prepare Move table from guppy for easier indexing
        # index Move table using aligned flank sequences
        moves_r = transform_moves(fast5.moves)
        left_raw = extract_from_moves(
            moves_r, left_alignment.position, fast5.strand_start, fast5.block_size)
        right_raw = extract_from_moves(
            moves_r, right_alignment.position, fast5.strand_start, fast5.block_size)

        seq = fast5.fasta[left_alignment.position.end:right_alignment.position.start]
        if seq and reverse:
            seq = get_reverse_strand(seq)

    return FlankInRead(
        read_id=readname,
        lflank_raw=left_raw,
        rflank_raw=right_raw,
        l_alignment=left_alignment,
        r_alignment=right_alignment,
        sequence=seq
    )


def load_overview_df(overview_path: str):
    df = pd.read_csv(overview_path)
    if df is None:
        raise RuntimeError(f'Loaded empty overview file from path={overview_path}')
    df.set_index('read_name', inplace=True)
    df.columns = df.columns.map(str)
    return df


def extract_tr_all(locus: Locus):
    """
    Extraction of tandem repeat regions using flanks of the locus.
    In the overview file of the locus following is stored:
    -Starting and ending positions of both flanks in the basecalling sequence
    -Starting and ending positions of both flanks in the raw fast5
    Alignments of flanks with the basecalled reads are also stored
    :param locus_path: Path to the extracted data for selected locus
    """

    # get path to overview file, where positions will be outputted
    overview_path = os.path.join(locus.path, tmpl.OVERVIEW_NAME)

    # get path where alignments will be stored
    align_path = os.path.join(locus.path,  tmpl.ALIGN_SUBDIR)

    # load existing overview file where is basic info about reads
    df = load_overview_df(overview_path)

    # prepare reads into list for parallelization
    read_list = []
    for row in df.itertuples():
        region = ReadForExtraction(
            name=row.Index,
            run_id=row.run_id,
            reverse=row.reverse,
            approx_location=row.sam_dist
        )
        read_list.append((region, locus))

    # run TR extraction
    dbg_msg = 'Extracting TRs - running with {thr} threads'.format(thr=main_config.threads)
    handle_msg_dbg(dbg_msg)
    if main_config.threads > 1:
        with Pool(main_config.threads) as p:
            results = p.map(extract_tr, read_list)
    else:
        results = [extract_tr(row) for row in read_list]

    new_col = []
    for row in df.itertuples():
        for res in results:
            if res.read_id == row.Index:

                # parse result for future storage in overview file
                new_col.append(process_row(res))

                # store alignment result as well
                save_alignment(row, res.alignment_mappings, align_path)

    # store positions of extracted tandem repeat regions in the overview file
    new_col = np.array(new_col).reshape((-1, len(new_col[0]))).T
    df = df.assign(saved=new_col[0],
                   l_start_raw=new_col[1],
                   l_end_raw=new_col[2],
                   r_start_raw=new_col[3],
                   r_end_raw=new_col[4],
                   l_seq_start=new_col[5],
                   l_seq_end=new_col[6],
                   r_seq_start=new_col[7],
                   r_seq_end=new_col[8])

    df.to_csv(os.path.join(overview_path))

    filename = os.path.join(locus.path, tmpl.PREDICTIONS_SUBDIR, 'basecalls', 'all.fasta')
    store_fasta(filename, results)


def store_fasta(filename: str, fasta: List[FlankInRead]):
    with open(filename, 'w') as file:
        for read in fasta:
            if read.sequence:
                file.write(f'>{read.read_id}\n')
                file.write(f'{read.sequence}\n\n')


def process_row(res: FlankInRead):
    return (
        res.valid,
        res.lflank_raw.start,
        res.lflank_raw.end,
        res.rflank_raw.start,
        res.rflank_raw.end,
        res.l_alignment.position.start,
        res.l_alignment.position.end,
        res.r_alignment.position.start,
        res.r_alignment.position.end
    )


def save_alignment(row, alignments: List[Mapping], align_path: str):
    """
    Write alignment into file using sample for discrimination
    :param row: Path to the extracted data for selected locus
    :param alignments: alignments to be written, list containing usually formatted alignment and string with score
    :param align_path: path where alignment will be appended
    """
    m_path = os.path.join(align_path, row.sample+'mappings.txt')
    if os.path.exists(m_path) is False:
        with open(m_path, 'w') as f:
            f.write(tmpl.ALIGNMENT_SEPARATOR)
            f.write(tmpl.FASTA_READNAME_ID.format(
                read_id=row.Index, reverse=str(row.reverse)))
            for a in alignments:
                f.write(tmpl.ALIGNMENT.format(
                    text=a.ref, equal=a.mapping, query=a.query))
            f.write(tmpl.ALIGNMENT_SEPARATOR)
    else:
        with open(m_path, 'a') as f:
            f.write(tmpl.FASTA_READNAME_ID.format(
                read_id=row.Index, reverse=str(row.reverse)))
            for a in alignments:
                f.write(tmpl.ALIGNMENT.format(
                    text=a.ref, equal=a.mapping, query=a.query))
            f.write(tmpl.ALIGNMENT_SEPARATOR)
