import os
from datetime import datetime

import src.helpers as aux
import src.templates as tmpl
from src.dtw_automata import main_wrapper
from src.extractor import extract_reads, extract_tr_all
from src.guppy import guppy_annotate
from src.input_handler.input import inputs, main_config, warpstr_config
from src.input_handler.locus import Locus
from src.report import run_genotyping_complex, run_genotyping_overview
from src.squiggler import Squiggler

# Main TRcall algorithm code starts below

start_time = aux.get_start_time()
script_path, _ = os.path.split(os.path.abspath(__file__))

# prepare output path
if os.path.isdir(main_config.output) is False:
    os.makedirs(main_config.output)
if main_config.exp_signal_generation:
    squiggler = Squiggler(main_config.reference_path)
    print('vytvaram squiggler')
if main_config.genotyping:
    muscle_path = os.path.join(script_path, 'example', 'muscle3.8.31_i86linux64')

# run main extraction of fast5 files per each locus
total_loci = len(warpstr_config['loci'])
for idx, obj in enumerate(warpstr_config['loci'], start=1):
    log = f"Processing {idx} of {total_loci} Locus name: {obj['name']}"
    locus = Locus(**obj)
    start_locus_time = datetime.now()
    last_time = start_locus_time
    print(tmpl.LOG_MSG_TIME.format(start=start_locus_time, log=log))

    # prepare all possible output subdirs
    # aux.prepare_subdirs(locus_path, config)

    # 1st step: finding reads mapped to the desired locus and extracting them as single fast5s
    if main_config.single_read_extraction:
        for inp in inputs:
            samples = inp.runs.split(',')
            extract_reads(inp.path, locus.coord, samples, locus.path)
        print('Finished extraction of single fast5 reads mapped to the locus', locus.name)
        # last_time = aux.print_time_duration(start_locus_time, 'read extraction', last_time)

    # 2nd step: call guppy basecalling with annotation
    if main_config.guppy_annotation:
        guppy_annotate(script_path, locus.path, main_config.threads)
        # last_time = aux.print_time_duration(start_locus_time, 'guppy annotation', last_time)

    # generate expected signals
    if main_config.exp_signal_generation:
        print('pustam squiggler')
        squiggler.process_locus(locus)
        # last_time = aux.print_time_duration(start_locus_time, 'expected signal generation', last_time)

    # extract tandem repeat regions
    if main_config.tr_region_extraction:
        extract_tr_all(locus)
        # last_time = aux.print_time_duration(start_locus_time, 'tr extraction', last_time)

    # run tandem repeat length calling
    overview_df, collapsed_df = None, None
    if main_config.tr_region_calling:
        overview_df, collapsed_df = main_wrapper(locus)
        # last_time = aux.print_time_duration(start_locus_time, 'tr calling', last_time)

    # run genotyping
    if main_config.genotyping:
        run_genotyping_overview(overview_df, locus.path, muscle_path)
        run_genotyping_complex(locus.path, locus, collapsed_df)
        last_time = aux.print_time_duration(start_locus_time, 'genotyping', last_time)

    aux.print_time_duration(start_locus_time, f'locus={locus.name} processing', None)


end_time = datetime.now()
t_seconds = (end_time-start_time).total_seconds()
hours = t_seconds//3600
minutes = (t_seconds-(hours*3600))//60
seconds = (t_seconds - minutes * 60)
print(tmpl.DURATION_PRINT.format(process='whole', hours=int(hours), minutes=int(minutes), seconds=int(seconds)))
