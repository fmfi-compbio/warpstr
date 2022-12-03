import pysam
from Bio.Seq import Seq


def get_sequences(coord: str, ref_path: str, flank_length: int):
    """
    Gets all reference sequences
    """
    lf_t, rf_t = get_flanks(coord, ref_path, flank_length, reverse=False)
    lf_r, rf_r = get_flanks(coord, ref_path, flank_length, reverse=True)
    tmpp, revp = get_ref_pattern(coord, ref_path)
    return (lf_t, rf_t, lf_r, rf_r, tmpp, revp)


def process_coord(coord: str):
    """
    Break down genomic coordinates represented as string into separate coordinates
    :param coord: genomic coordinates
    """
    chromosome = coord.split(':')[0]
    start = int(coord.split(':')[1].split('-')[0].replace(',', ''))
    end = int(coord.split('-')[1].replace(',', ''))
    return chromosome, start, end


def get_reverse_strand(seq: str) -> str:
    """
    Reverses sequence
    """
    return str(Seq(seq).reverse_complement())


def get_ref_pattern(coord: str, ref: str):
    """
    Gets reference pattern
    """
    temp_ref_pattern = pysam.faidx(ref, coord)
    temp_ref_pattern = ''.join(temp_ref_pattern.splitlines()[1:]).upper()
    rev_ref_pattern = get_reverse_strand(temp_ref_pattern)
    return temp_ref_pattern, rev_ref_pattern


def get_flanks(coord: str, ref: str, flank_length: int, reverse: bool, motifpart: int = 0):
    """
    Gets flanking sequences
    """
    chrom, start, end = process_coord(coord)

    lflank_coord = chrom+':'+str(start-flank_length)+'-'+str(start+motifpart-1)
    lflank = pysam.faidx(ref, lflank_coord)
    rflank_coord = chrom+':'+str(end+1-motifpart)+'-'+str(end+flank_length)
    rflank = pysam.faidx(ref, rflank_coord)

    lflank = ''.join(lflank.splitlines()[1:]).upper()
    rflank = ''.join(rflank.splitlines()[1:]).upper()

    if reverse:
        lflank, rflank = get_reverse_strand(rflank), get_reverse_strand(lflank)
    return lflank, rflank
