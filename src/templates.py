TIME_PRINT = '''--- Initial time of WarpSTR:{start: %Y-%m-%d %H:%M:%S} ---'''
LOG_MSG_TIME = '''{start:%Y-%m-%d %H:%M:%S} {log}'''
DURATION_PRINT = '''  Duration for {process}: {hours:02}h {minutes:02}m {seconds:02}s'''
ERR_CREATE_DIR = '''Error when creating directory: {directory}'''
ERR_CREATE_FILE = '''Error when creating file: {file}'''
ERR_MSG = '''ERROR_{lvl}: {msg}'''
DBG_MSG = '''DEBUG_{lvl}: {msg}'''

#subdir names in output data structure
FAST5_SUBDIR = "fast5" #subdirectory where single fast5s are extracted
ANNOT_SUBDIR = "annot" #subdirectory where guppy outputs annotated fast5s
OVERVIEW_NAME = "overview.csv" #where info about reads in the locus path is stored
ALIGN_SUBDIR = "alignments" #subdirectory where alignments are stored
PREDICTIONS_SUBDIR = "predictions"
WARPS = "DTW_alignments"
LOCUS_INFO_SUBDIR = "expected_signals"
REPORTS_SUBDIR = "reports"
SUMMARY_SUBDIR = "summaries"
COMPLEX_SUBDIR = "complexSTR_analysis"

#reports
REPORT_FILENAME='''{locus}_report.html'''
REPORT_HEADER='''Main report for {locus} locus'''
REPORT_SUCCESS='''  Report saved as {path}'''
REPORT_FAILURE='''Error when saving report as {path}'''

#reports params
MSA_SCRIPT='<script src="https://s3.eu-central-1.amazonaws.com/cdn.bio.sh/msa/latest/msa.min.gz.js"></script>'

#squiggler part
LOCUS_FLANKS = "sequences.csv"
LOCUS_NAMES = ["left_flank_template", "right_flank_template", "left_flank_reverse","right_flank_reverse",\
               "temp_ref_pattern", "rev_ref_pattern"]

#alignment part
FASTA_HEAD = '''***Settings for alignment:\n*Match score: {match}\n*Mismatch score: {mismatch}\n*Gap open penalty: {gap_open}\n*Gap extend penalty: {gap_extend}\n*Accuracy factor: {acc_factor}\n'''
FASTA_READNAME_ID = '''>{read_id}_{reverse}\n'''
ALIGNMENT_MATCH_CHAR = "|"
ALIGNMENT_MISMATCH_CHAR = "-"
ALIGNMENT_SUFFIX_FILE = "_mappings.txt"
ALIGNMENT = '''{score}\n{text}\n{equal}\n{query}\n'''
ALIGNMENT_SEPARATOR = "*"*100+"\n"

#graph for state automata
NODE_SIZE = 400
NODE_ALPHA = 0.5
EDGE_WIDTH = 1
EDGE_ALPHA = 0.5
ARROW_STYLE = "->"
ARROW_SIZE = 25
COLOR_MAP = {0:"yellow",1:"red",2:"blue",3:"green",4:"cyan",5:"k",6:"indigo",7:"gold",8: "lavender",9:"maroon"} 
COLOR_EDGES = {"straight": "blue","repeat":"indigo","forward":"green"}
CONN_STYLES = {"straight":None,"repeat":"arc3,rad=-0.25","forward":"arc3,rad=-0.25"}
LINE_STYLES = {"straight":"solid","repeat":"dashed","forward":"dashdot"}


#DNA alphabet conversion
DNA_DICT = {
    "R":["A","G"],
    "Y":["C","T"],
    "S":["G","C"],
    "W":["A","T"],
    "K":["G","T"],
    "M":["A","C"],
    "B":["C","G","T"],
    "D":["A","G","T"],
    "H":["A","C","T"],
    "V":["A","C","G"],
    "N":["A","C","G","T"]
}
ENCODING_DICT = {
    "(" : ")",
    ")" : "(",
    "{" : "}",
    "}" : "{",
    "A" : "T",
    "T" : "A",
    "G" : "C",
    "C" : "G",
    "M" : "K",
    "K" : "M",
    "N" : "N",
    "W" : "W",
    "S" : "S",
    "R" : "Y",
    "Y" : "R",
    "B" : "V",
    "D" : "H",
    "H" : "D",
    "V" : "B"
}