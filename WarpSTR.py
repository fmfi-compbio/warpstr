import os
from datetime import datetime
from subprocess import call

from src.squiggler import Squiggler
from src.input_handler import load_args
from src.extractor import extract_reads, extract_tr_all
from src.dtw_automata import main_wrapper
from src.report import create_report, run_genotyping_overview, run_genotyping_complex
import src.templates as tmpl
import src.helpers as aux
from src.prepare_sequence import prepare_sequence

# Main TRcall algorithm code starts below

start_time = aux.get_start_time()
script_path,_ = os.path.split(os.path.abspath(__file__))

#load args from config
config, config_path = load_args(script_path)
threads = config['threads']
verbose = config['verbose']

#prepare output path
main_out_path = config['output']
if os.path.isdir(main_out_path) is False:
    os.makedirs(main_out_path)
    if os.path.isdir(main_out_path) is False:
        raise FileNotFoundError('Error when creating directory: {directory}')

reports_path = os.path.join(main_out_path,tmpl.REPORTS_SUBDIR)
if os.path.exists(reports_path) is False:
    os.mkdir(reports_path)

src_report_path = os.path.join(script_path,"src","report")

#run main extraction of fast5 files per each locus
loci_counter = 0
for locus in config['loci']:
    loci_counter += 1
    log = "Processing {count} of {total} Locus name: {locus_name}".format(
        count = loci_counter,
        total = len(config['loci']),
        locus_name = locus['name'])
    start_locus_time = datetime.now()
    last_time = start_locus_time
    print(tmpl.LOG_MSG_TIME.format(start=start_locus_time, log = log))

    if 'flank_length' not in locus:
        locus['flank_length'] = config['flank_length']
        print(f'Flank length not set for locus - using default value of {config["flank_length"]}')

    if 'sequence' not in locus:
        locus = prepare_sequence(locus, config['reference_path'])
        print(f'Sequence was not set for {locus["name"]}. Automatic configuration defined sequence as:\n{locus["sequence"]}\nderived from reference sequence {locus["noting"]}')
    
    #prepare all possible output subdirs
    locus_path = os.path.join(main_out_path,locus['name'])
    fast5_out_path = os.path.join(locus_path,tmpl.FAST5_SUBDIR)
    aux.prepare_subdirs(locus_path, config)

    #1st step: finding reads mapped to the desired locus and extracting them as single fast5s
    if config['single_read_extraction']:
        for data_path in config['inputs']:
            samples = data_path['samples'].split(',')
            results_dict = extract_reads(data_path['path'],locus['coord'],samples, locus_path,threads=threads) 
            print("Finished extraction of single fast5 reads mapped to the locus",locus['name'])
            if verbose==1:
                for result in results_dict:
                    for key, value in result:
                        print(key, ' : ', value)
    
            last_time = aux.print_time_duration(start_locus_time,"read extraction",last_time,verbose)
    
    #2nd step: call guppy basecalling with annotation
    if config['guppy_annotation']:
        guppy_wrapper = os.path.join(script_path,"src","guppy_annotate","wrapper.sh")
        guppy = config['guppy_config']
        if os.path.exists(guppy['path']) is False:
            raise ValueError(f'Could not load Guppy from {guppy["path"]}')

        try:
            ret_call = call([guppy_wrapper,"-g",guppy['path'],"-r",fast5_out_path,"-f",guppy['flowcell'],"-k",guppy['kit'] \
                           ,"-a",tmpl.ANNOT_SUBDIR,"-t",str(threads)])
            if ret_call < 0:
                print("Execution of guppy failed in locus",locus_path)
                raise ValueError("Running Guppy was terminated by signal", -ret_call)
        except OSError as e:
            print(str(e),"Execution of guppy failed in locus",locus_path)
            raise

        last_time = aux.print_time_duration(start_time,"guppy annotation",last_time,verbose)

    # generate expected signals
    if config['exp_signal_generation']:
        squiggler = Squiggler(config['pore_model_path'], config['reference_path'])
        squiggler.process_locus(locus_path, locus['flank_length'], locus['coord'])
        last_time = aux.print_time_duration(start_locus_time,"expected signal generation",last_time,verbose)

    # extract tandem repeat regions
    if config['tr_region_extraction']:
        extract_tr_all(locus_path, config['alignment'], threads)
        last_time = aux.print_time_duration(start_locus_time,"tr extraction",last_time,verbose)

    # run tandem repeat length calling
    overview_df, collapsed_df = None, None
    if config['tr_region_calling']:
        tr_results, overview_df, collapsed_df = main_wrapper(locus_path, config['pore_model_path'],\
                                        locus,config['tr_calling_config'],threads)
        create_report(main_out_path,tr_results,locus,locus_path,src_report_path)
        last_time = aux.print_time_duration(start_locus_time,"tr calling",last_time,verbose)

    # run genotyping
    if config['genotyping']:
        muscle_path = os.path.join(script_path,"example","muscle3.8.31_i86linux64")    
        run_genotyping_overview(overview_df,locus_path, config['genotyping_config'], muscle_path)
        run_genotyping_complex(locus_path, config['genotyping_config'], locus, collapsed_df)
        last_time = aux.print_time_duration(start_locus_time,"genotyping",last_time,verbose)

    aux.print_time_duration(start_locus_time,"locus processing",None, 5)


end_time = datetime.now()
t_seconds = (end_time-start_time).total_seconds()
hours = t_seconds//3600
minutes = (t_seconds-(hours*3600))//60
seconds = (t_seconds - minutes * 60)
print(tmpl.DURATION_PRINT.format(process="whole",hours=int(hours),minutes=int(minutes),seconds=int(seconds)))
