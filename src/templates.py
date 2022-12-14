TIME_PRINT = '''--- Initial time of WarpSTR:{start: %Y-%m-%d %H:%M:%S} ---'''
LOG_MSG_TIME = '''{start:%Y-%m-%d %H:%M:%S} {log}'''
DURATION_PRINT = '''  Duration for {process}: {hours:02}h {minutes:02}m {seconds:02}s'''

ERR_MSG = '''ERROR_{lvl}: {msg}'''
DBG_MSG = '''DEBUG_{lvl}: {msg}'''

# subdir names in output data structure
FAST5_SUBDIR = 'fast5'  # subdirectory where single fast5s are extracted
ANNOT_SUBDIR = 'annot'  # subdirectory where guppy outputs annotated fast5s
OVERVIEW_NAME = 'overview.csv'  # where info about reads in the locus path is stored
ALIGN_SUBDIR = 'alignments'  # subdirectory where alignments are stored
PREDICTIONS_SUBDIR = 'predictions'
WARPS = 'DTW_alignments'
LOCUS_INFO_SUBDIR = 'expected_signals'
SUMMARY_SUBDIR = 'summaries'
COMPLEX_SUBDIR = 'complexSTR_analysis'

# squiggler part
LOCUS_FLANKS = 'sequences.csv'
LOCUS_NAMES = ['left_flank_template', 'right_flank_template', 'left_flank_reverse', 'right_flank_reverse',
               'temp_ref_pattern', 'rev_ref_pattern']

# alignment part
FASTA_READNAME_ID = '''>{read_id}_{reverse}\n'''
ALIGNMENT_MATCH_CHAR = '|'
ALIGNMENT_MISMATCH_CHAR = '-'
ALIGNMENT = '''{text}\n{equal}\n{query}\n'''
ALIGNMENT_SEPARATOR = '*'*100+'\n'

# DNA alphabet conversion
DNA_DICT = {
    'R': ['A', 'G'],
    'Y': ['C', 'T'],
    'S': ['G', 'C'],
    'W': ['A', 'T'],
    'K': ['G', 'T'],
    'M': ['A', 'C'],
    'B': ['C', 'G', 'T'],
    'D': ['A', 'G', 'T'],
    'H': ['A', 'C', 'T'],
    'V': ['A', 'C', 'G'],
    'N': ['A', 'C', 'G', 'T']
}
ENCODING_DICT = {
    '(': ')',
    ')': '(',
    '{': '}',
    '}': '{',
    'A': 'T',
    'T': 'A',
    'G': 'C',
    'C': 'G',
    'M': 'K',
    'K': 'M',
    'N': 'N',
    'W': 'W',
    'S': 'S',
    'R': 'Y',
    'Y': 'R',
    'B': 'V',
    'D': 'H',
    'H': 'D',
    'V': 'B'
}
