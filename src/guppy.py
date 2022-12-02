import os
import src.templates as tmpl
from typing import Dict
from subprocess import call

def guppy_annotate(script_path:str, locus_path: str, threads: int, guppy:Dict[str,str]):
    """Call script running Guppy annotate for all fast5 in locus path

    Args:
        script_path (str): path to the WarpSTR source directory
        locus_path (str): path to the locus output directory
        threads (int): number of threads to use
        guppy (Dict[str,str]): guppy config

    Raises:
        FileNotFoundError: Path to Guppy was incorrectly set
        ValueError: Calling Guppy exited with failure status using subprocess
        OSError: Calling subprocess failed with OSError due to I/O problems
    """
    guppy_wrapper = os.path.join(script_path, "src","guppy_annotate","wrapper.sh")
    fast5_out_path = os.path.join(locus_path, tmpl.FAST5_SUBDIR)
    if os.path.exists(guppy['path']) is False:
        raise FileNotFoundError(f'Could not load Guppy from {guppy["path"]}')

    command = [
        guppy_wrapper,
        "-g",guppy['path'],
        "-r",fast5_out_path,
        "-f",guppy['flowcell'],
        "-k",guppy['kit'],
        "-a",tmpl.ANNOT_SUBDIR,
        "-t",str(threads)
        ]
    try:
        ret_call = call(command)
        if ret_call < 0:
            raise ValueError(f"Execution of guppy failed in locus={locus_path} by return signal {-ret_call}")
    except OSError as e:
        raise OSError(f"Execution of guppy failed in locus={locus_path} with err={e}")
