import os
import shutil
from datetime import datetime

from src.config import main_config

from . import templates as tmpl


def get_start_time():
    start_time = datetime.now()
    print(tmpl.TIME_PRINT.format(start=start_time))
    return start_time


def print_time_duration(start_time: datetime, process: str, second_time: datetime):
    if second_time is None or second_time == start_time:
        use_time = start_time
    else:
        use_time = second_time

    end_time = datetime.now()
    t_seconds = (end_time-use_time).total_seconds()
    hours = t_seconds//3600
    min = (t_seconds-(hours*3600))//60
    sec = (t_seconds - min * 60)
    if main_config.verbose > 0:
        print(tmpl.DURATION_PRINT.format(process=process, hours=int(hours), minutes=int(min), seconds=int(sec)))
    return end_time


def handle_subdir(subdir_path: str, overwrite_flag: bool):
    if os.path.exists(subdir_path) is False:
        os.mkdir(subdir_path)
    elif overwrite_flag:
        shutil.rmtree(subdir_path)
        os.mkdir(subdir_path)


def prepare_subdirs(locus_path: str):
    predictions_subdir = os.path.join(locus_path, tmpl.PREDICTIONS_SUBDIR)
    align_subdir = os.path.join(locus_path,  tmpl.ALIGN_SUBDIR)
    locus_info_subdir = os.path.join(locus_path, tmpl.LOCUS_INFO_SUBDIR)
    fast5_subdir = os.path.join(locus_path, tmpl.FAST5_SUBDIR)
    summaries_subdir = os.path.join(locus_path, tmpl.SUMMARY_SUBDIR)

    if os.path.exists(locus_path) is False:
        os.mkdir(locus_path)
    if os.path.exists(summaries_subdir) is False:
        os.mkdir(summaries_subdir)
    if main_config.single_read_extraction:
        handle_subdir(fast5_subdir, main_config.force_overwrite)
        overview_file = os.path.join(locus_path, tmpl.OVERVIEW_NAME)
        if os.path.exists(overview_file):
            os.remove(overview_file)
        with open(overview_file, 'w') as f:
            f.write('read_name,sample,run_id,reverse,sam_dist\n')
    if main_config.tr_region_extraction:
        handle_subdir(align_subdir, main_config.force_overwrite)
    if main_config.exp_signal_generation:
        handle_subdir(locus_info_subdir, main_config.force_overwrite)
    if main_config.tr_region_calling:
        handle_subdir(predictions_subdir, main_config.force_overwrite)
        DTW_alignments = os.path.join(predictions_subdir, tmpl.WARPS)
        handle_subdir(DTW_alignments, main_config.force_overwrite)
        basecalls = os.path.join(predictions_subdir, 'basecalls')
        handle_subdir(basecalls, main_config.force_overwrite)
        seqs = os.path.join(predictions_subdir, 'sequences')
        handle_subdir(seqs, main_config.force_overwrite)
