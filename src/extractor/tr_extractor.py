import os
import sys
from multiprocessing import Pool
from typing import Tuple

import numpy as np
import pandas as pd
from Bio import pairwise2
from mappy import revcomp

import src.templates as tmpl
from src.input_handler.fast5 import Fast5


def save_alignment(row, alignments, align_path, config):
    """
    Write alignment into file using sample for discrimination
    :param row: Path to the extracted data for selected locus
    :param alignments: alignments to be written, list containing usually formatted alignment and string with score
    :param align_path: path where alignment will be appended
    """
    m_path = os.path.join(align_path, row.sample+tmpl.ALIGNMENT_SUFFIX_FILE)
    if os.path.exists(m_path) is False:
        with open(m_path, 'w') as f:
            f.write(tmpl.ALIGNMENT_SEPARATOR)
            f.write(tmpl.FASTA_HEAD.format(match=config['match_score'], mismatch=config['mismatch_score'],
                                           gap_open=config['gap_open_score'], gap_extend=config['gap_extend_score'],
                                           acc_factor=config['accuracy_factor']))
            f.write(tmpl.ALIGNMENT_SEPARATOR)
            f.write(tmpl.FASTA_READNAME_ID.format(
                read_id=row.Index, reverse=str(row.reverse)))
            for a in alignments:
                f.write(tmpl.ALIGNMENT.format(
                    score=a[0], text=a[1], equal=a[2], query=a[3]))
            f.write(tmpl.ALIGNMENT_SEPARATOR)
    else:
        with open(m_path, 'a') as f:
            f.write(tmpl.FASTA_READNAME_ID.format(
                read_id=row.Index, reverse=str(row.reverse)))
            for a in alignments:
                f.write(tmpl.ALIGNMENT.format(
                    score=a[0], text=a[1], equal=a[2], query=a[3]))
            f.write(tmpl.ALIGNMENT_SEPARATOR)


def load_flanks(locus_path: str):
    """
    Function loading template and reverse flanks
    :param locus_path: Path to the extracted data for selected locus
    :returns (lft,rft): left and right flanking template sequences
    :returns (lfr,rfr): left and right flanking reverse sequences
    """
    path = os.path.join(locus_path, tmpl.LOCUS_INFO_SUBDIR, tmpl.LOCUS_FLANKS)
    lft, rft, lfr, rfr = '', '', '', ''

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

    return ((lft, rft), (lfr, rfr))


def handle_msg_dbg(dbg_msg):
    print(tmpl.DBG_MSG.format(lvl='3_TREX', msg=dbg_msg), file=sys.stderr)


def transform_moves(moves):
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


def extract_from_moves(moves_r, pos1, pos2, strand_start, block_stride):
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
    decoded_pos1 = np.where(moves_r == pos1)[0]
    decoded_pos2 = np.where(moves_r == pos2)[0]

    if len(decoded_pos1) > 0:
        start = strand_start + decoded_pos1[0]*block_stride
    else:
        start = -1

    if len(decoded_pos2) > 0:
        end = strand_start + decoded_pos2[-1]*block_stride
    else:
        end = -1

    return start, end


def find_sequence(seq1, seq2, match_score, mis_score, gap_open, gap_ext):
    """
    Extraction of tandem repeat region using flank sequences.
    :param seq1: text sequence
    :param seq2: pattern sequence
    :param match_score: score if two nucleotides matched
    :param mis_score: penalty if two nucleotides do not match
    :param gap_open: penalty for opening gap
    :param gap_ext: penalty for extending gap
    :return score: score of the alignment
    :return real_start: start of the alignment in the basecalled sequence
    :return end: calculated alignments of flanks with basecalled read
    :return (score_string,ref,mapping,query): calculated alignments of flanks with basecalled read
    """

    # get the best alignment using local alignment from pairwise lib
    al1, al2, score, start, end = pairwise2.align.localms(
        seq1, seq2, match_score, mis_score, gap_open, gap_ext, one_alignment_only=True)[0]

    # decode the start and end positions of the alignment so alignment can be written to file later
    # this also handles cases when alignment starts or ends earlier
    # and ensures that the alignment is full length with the input pattern
    nums_gaps = al1[start:end].count('-')
    nums_gaps2 = al2[start:end].count('-')
    real_start = al2[:start].count('-')
    end = real_start+len(seq2)+nums_gaps2-nums_gaps

    # prepare alignment as list of strings
    ref = al1[real_start:end]
    query = al2[real_start:end]
    mapping = ''
    identity = 0
    for i, j in zip(ref, query):
        if i == j:
            mapping += tmpl.ALIGNMENT_MATCH_CHAR
            identity += 1
        else:
            mapping += tmpl.ALIGNMENT_MISMATCH_CHAR
    # Save also score into the alignment
    diff = (len(seq2)-(len(query)-nums_gaps2))*gap_ext
    score = score + diff
    score_string = 'Score: '+str(score)+' Identity: '+str(identity)+'/'+str(len(ref))
    identity = identity/(len(ref))
    return (score, identity), real_start, end, (score_string, ref, mapping, query)


def align_seq(read, flanks, config, return_alignment=False):
    """
    Extraction of tandem repeat region using flank sequences.
    :param read: basecalled read sequence
    :param flanks: list of flank sequences
    :param config: dictionary of scoring strategy and accuracy factor determining min. quality of alignments
    :return (start_l, end_l): start and end position of left flank sequence
    :return (start_r, end_r): start and end position of right flank sequence
    :return (alignment_l,alignment_r): calculated alignments of flanks with basecalled read
    """
    match = int(config['match_score'])
    mis = int(config['mismatch_score'])
    gap_o = int(config['gap_open_score'])
    gap_e = int(config['gap_extend_score'])

    # align basecalled sequence with flanks
    score_l, start_l, end_l, alignment_l = find_sequence(
        read, flanks[0], match, mis, gap_o, gap_e)

    if (end_l+len(flanks[1])) > len(read):
        score_r, start_r, end_r, alignment_r = (-1, -1), -1, -1, ('not found', '-', '-', '-')
    else:
        score_r, start_r, end_r, alignment_r = find_sequence(
            read[end_l:], flanks[1], match, mis, gap_o, gap_e)
    # check if alignment was good enough using accuracy factor
    if score_l[0] <= len(flanks[0])*float(config['accuracy_factor']) and score_l[1] <= float(config['identity_factor']):
        start_l, end_l = -1, -1
    if score_r[0] <= len(flanks[1])*float(config['accuracy_factor']) and score_r[1] <= float(config['identity_factor']):
        start_r, end_r = -1, -1
    else:
        start_r += end_l
        end_r += end_l

    return (start_l, end_l), (start_r, end_r), (alignment_l, alignment_r)


def extract_tr(overview_row):
    readname = overview_row[0]
    runid = overview_row[1]
    reverse = overview_row[2]
    approx_location = overview_row[3]

    fast5path = os.path.join(curr_locus_path, tmpl.FAST5_SUBDIR, str(
        runid), tmpl.ANNOT_SUBDIR, readname+'.fast5')

    dbg_msg = 'processing {read}'.format(read=readname)
    handle_msg_dbg(dbg_msg)

    # decode flanks according to strand
    flanks = flanks_rev if reverse else flanks_tem

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

    dbg_msg = 'aligning flank with seq[{start}:{end}] of total {ln} in {read}'.format(
        start=start, end=end, ln=len(fast5.fasta), read=readname)
    handle_msg_dbg(dbg_msg)

    try:
        ref_text = fast5.fasta[start:end]
    except Exception as e:
        print('ERROR', str(e))

    # align flanks to the read
    (l_begin_seq, l_end_seq), (r_begin_seq, r_end_seq), alignments = align_seq(
        ref_text, flanks, curr_config_align)
    if l_begin_seq != -1:
        l_begin_seq += start
    if l_end_seq != -1:
        l_end_seq += start
    if r_begin_seq != -1:
        r_begin_seq += start
    if r_end_seq != -1:
        r_end_seq += start

    if l_begin_seq == -1 or l_end_seq == -1 or r_begin_seq == -1 or r_end_seq == -1:
        l_start, l_end = -1, -1
        r_start, r_end = -1, -1
        seq = ''
    else:
        # prepare Move table from guppy for easier indexing
        # index Move table using aligned flank sequences
        dbg_msg = 'delineating TR [{b1},{e1}],[{b2},{e2}] for {read}'.format(
            b1=l_begin_seq, e1=l_end_seq, b2=r_begin_seq, e2=r_end_seq, read=readname)
        handle_msg_dbg(dbg_msg)

        moves_r = transform_moves(fast5.moves)

        l_start, l_end = extract_from_moves(
            moves_r, l_begin_seq, l_end_seq, fast5.strand_start, fast5.block_size)
        r_start, r_end = extract_from_moves(
            moves_r, r_begin_seq, r_end_seq, fast5.strand_start, fast5.block_size)

        seq = fast5.fasta[l_end_seq:r_begin_seq]
        if reverse:
            seq = revcomp(seq)

    dbg_msg = 'finished TR extraction for {read}'.format(read=readname)
    handle_msg_dbg(dbg_msg)

    return {'read_id': readname,
            'lflank': (l_start, l_end),
            'rflank': (r_start, r_end),
            'align_res': ((l_begin_seq, l_end_seq), (r_begin_seq, r_end_seq), alignments),
            'sequence': seq,
            'reverse': reverse}


def load_overview_df(overview_path: str):
    df = pd.read_csv(overview_path)
    if not df:
        raise RuntimeError(f'Loaded empty overview file from path={overview_path}')
    df.set_index('read_name', inplace=True)
    df.columns = df.columns.map(str)
    return df


def extract_tr_all(locus_path: str, config_align, threads: int):
    """
    Extraction of tandem repeat regions using flanks of the locus.
    In the overview file of the locus following is stored:
    -Starting and ending positions of both flanks in the basecalling sequence
    -Starting and ending positions of both flanks in the raw fast5
    Alignments of flanks with the basecalled reads are also stored
    :param locus_path: Path to the extracted data for selected locus
    """
    global flanks_tem
    global flanks_rev
    global curr_config_align
    global curr_locus_path

    curr_locus_path = locus_path
    curr_config_align = config_align

    flanks = load_flanks(locus_path)
    if flanks:
        flanks_tem, flanks_rev = flanks[0], flanks[1]
    else:
        return False

    # get path to overview file, where positions will be outputted
    overview_path = os.path.join(locus_path, tmpl.OVERVIEW_NAME)

    # get path where alignments will be stored
    align_path = os.path.join(locus_path,  tmpl.ALIGN_SUBDIR)

    # load existing overview file where is basic info about reads
    df = load_overview_df(overview_path)

    # prepare reads into list for parallelization
    read_list = []
    for row in df.itertuples():
        read_list.append((row.Index, row.run_id, row.reverse, row.sam_dist))

    # run TR extraction
    dbg_msg = 'Extracting TRs - running with {thr} threads'.format(thr=threads)
    handle_msg_dbg(dbg_msg)
    if threads > 1:
        with Pool(threads) as p:
            results = p.map(extract_tr, read_list)
    else:
        results = [extract_tr(row) for row in read_list]

    new_col = []
    for row in df.itertuples():
        for res in results:
            if res['read_id'] == row.Index:

                # parse result for future storage in overview file
                new_col.append(process_row(res))

                # store alignment result as well
                save_alignment(row, res['align_res'][2], align_path, config_align)

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

    filename = os.path.join(locus_path, tmpl.PREDICTIONS_SUBDIR, 'basecalls', 'all.fasta')
    store_fasta(filename, results)

    filename = os.path.join(locus_path, tmpl.PREDICTIONS_SUBDIR, 'basecalls', 'basecalls_reverse.fasta')
    results_rev = [res for res in results if res['reverse']]
    store_fasta(filename, results_rev)

    filename = os.path.join(locus_path, tmpl.PREDICTIONS_SUBDIR, 'basecalls', 'basecalls_template.fasta')
    results_temp = [res for res in results if res['reverse'] is False]
    store_fasta(filename, results_temp)

    return True


def store_fasta(filename: str, fasta):
    with open(filename, 'w') as file:
        for i in fasta:
            if len(str(i['sequence'])) > 0:
                file.write('>'+i['read_id']+'\n')
                file.write(str(i['sequence'])+'\n\n')


def valid_region(start: int, end: int) -> bool:
    """
    Checks if input region is valid
    :param start: start of the region
    :param end: end of the region
    """
    if start > end:
        return False
    elif start == -1 or end == -1:
        return False
    else:
        return True


def is_valid(left: Tuple[int, int], right: Tuple[int, int], l_start: int, l_end: int, r_start: int, r_end: int):
    """
    Checks if TR region extracted is valid
    :param left: tuple of start and end coordinate of left flank seq
    :param right: tuple of start and end coordinate of right flank seq
    :param l_start: start of left flank in raw signal
    :param l_end: end of left flank in raw signal
    :param r_start: start of right flank in raw signal
    :param r_end: end of right flank in raw signal
    """
    if valid_region(left[0], left[1]) and valid_region(right[0], right[1])\
            and valid_region(l_start, l_end) and valid_region(r_start, r_end)\
            and valid_region(l_start, r_end):
        saved = 1
    else:
        saved = 0
    return saved


def process_row(res):
    """
    Processes alignment result and saves interim results
    :param res: dict result of flanks alignment with read
    """

    left, right, alignments = res['align_res']
    (l_start, l_end) = res['lflank']
    (r_start, r_end) = res['rflank']

    saved = is_valid(left, right, l_start, l_end, r_start, r_end)

    # save positions in the list for future use
    return (saved, l_start, l_end, r_start, r_end, left[0], left[1], right[0], right[1])
