import argparse
import os
import pathlib
from dataclasses import dataclass
from typing import Optional

import yaml


def load_args(script_path: str):
    """
    Parses input args and returns parsed config
    """

    parser = argparse.ArgumentParser(description='WarpSTR')
    parser.add_argument('config', metavar='Config', type=is_valid_file, help='config file in .yaml')
    args = parser.parse_args()
    config = load_yaml(args.config)
    if config is None:
        raise ValueError(f'Error when loading config file from {args.config}')

    defaults_path = os.path.join(script_path, 'src', 'default.yaml')
    defaults = load_yaml(os.path.join(script_path, 'src', 'default.yaml'))
    if defaults is None:
        raise ValueError(f'Error when loading config file from {defaults_path}')

    add_defaults(config, defaults)
    return config


def add_defaults(config, default):
    """
    Look through the parameters and adds the default ones recursively.
    :param config: dict - dictionary with parameters
    :param default: dict - dictionary with default parameters
    :return: None
    """

    for key in default:
        if isinstance(default[key], dict):
            # go deeper
            if key not in config:
                config[key] = {}
            add_defaults(config[key], default[key])
        elif isinstance(default[key], list):
            # add it to all the configs
            if key in config:
                for item in config[key]:
                    add_defaults(item, default[key][0])
        elif isinstance(default[key], str) and default[key] == 'required':
            # required - check if it is available
            if key not in config:
                raise KeyError('Required argument not available - %s' % key)
            if config[key] is None or config[key].strip() == '':
                raise KeyError('Required argument is empty - %s' % key)
        else:
            # add it
            if key not in config:
                config[key] = default[key]


def is_valid_file(arg: str) -> str:
    """
    Checks if input config file exists
    :return: arg - path to config file
    """
    if not os.path.exists(arg):
        raise argparse.ArgumentTypeError(f'The config file {arg} does not exist!')
    return arg


def load_yaml(yaml_path: str):
    """
    Reads config info from the given yaml config file
    :param yaml_path: str - path to yaml parameters file
    :return: dict - dictionary with read params
    """
    loaded_yaml = None

    try:
        with open(yaml_path, 'r') as stream:
            loaded_yaml = yaml.load(stream, Loader=yaml.FullLoader)
    except yaml.YAMLError as exc:
        raise ValueError(f'Incorrect YAML format in config path {yaml_path},err={exc}')
    except IOError as err:
        raise IOError(f'IOError when processing config path {yaml_path}, err={err}')

    return loaded_yaml


@dataclass
class RescalerConfig:
    reps_as_one: bool = False
    threshold: float = 0.5
    max_std: float = 0.5
    method: str = 'mean'

    def __post_init__(self):
        assert self.threshold > 0
        assert self.max_std > 0
        assert self.method == 'mean' or self.method == 'median'


@dataclass
class CallerConfig:
    spike_removal: str = 'Brute'
    min_values_per_state: int = 4
    states_in_segment: int = 6
    min_state_similarity: float = 0.75
    visualize_alignment: bool = True
    visualize_phase: bool = True
    visualize_strand: bool = True
    visualize_cost: bool = True

    def __post_init__(self):
        assert self.min_values_per_state > 1
        assert self.states_in_segment > 1
        assert self.min_state_similarity > 0
        assert self.spike_removal in ['None', 'median3', 'median5', 'Brute']


@dataclass
class GenotypingConfig:
    min_weight: float = 0.2
    std_filter: float = 2
    visualize: bool = True
    msa: bool = False

    def __post_init__(self):
        assert self.min_weight > 0 and self.min_weight < 1
        assert self.std_filter > 1


@dataclass
class AlignmentConfig:
    accuracy_factor: float = 1.15
    identity_factor: float = 0.85
    match_score: int = 2
    mismatch_score: int = -3
    gap_open_score: int = -3
    gap_extend_score: int = -3


@dataclass
class GuppyConfig:
    path: Optional[str] = None
    flowcell: Optional[str] = None
    kit: Optional[str] = None


@dataclass
class Config:
    output: str
    reference_path: str
    pore_model_path: str
    single_read_extraction: bool
    guppy_annotation: bool
    exp_signal_generation: bool
    tr_region_extraction: bool
    tr_region_calling: bool
    genotyping: bool
    flank_length: int
    threads: int
    verbose: int = 0
    force_overwrite: bool = False


@dataclass
class InputConfig:
    path: str
    runs: str


ROOT_DIR_PATH = pathlib.Path().resolve()
warpstr_config = load_args(ROOT_DIR_PATH)

if 'loci' not in warpstr_config:
    raise KeyError('No loci defined in the config')

main_config = Config(
    output=warpstr_config['output'],
    pore_model_path=warpstr_config['pore_model_path'],
    exp_signal_generation=warpstr_config['exp_signal_generation'],
    reference_path=warpstr_config['reference_path'],
    single_read_extraction=warpstr_config['single_read_extraction'],
    guppy_annotation=warpstr_config['guppy_annotation'],
    tr_region_calling=warpstr_config['tr_region_calling'],
    tr_region_extraction=warpstr_config['tr_region_extraction'],
    genotyping=warpstr_config['genotyping'],
    threads=warpstr_config['threads'],
    verbose=warpstr_config['verbose'],
    force_overwrite=warpstr_config['force_overwrite'],
    flank_length=warpstr_config['flank_length']
)

inputs = [InputConfig(**i) for i in warpstr_config.get('inputs', [])]


print('NACITAVAM NASTAVENIE')
rescaler_config = RescalerConfig(**warpstr_config['rescaling'])
caller_config = CallerConfig(**warpstr_config['tr_calling_config'])
alignment_config = AlignmentConfig(**warpstr_config['alignment'])
genotyping_config = GenotypingConfig(**warpstr_config['genotyping_config'])

if 'guppy_config' in warpstr_config:
    guppy_config = GuppyConfig(**warpstr_config['guppy_config'])
else:
    guppy_config = GuppyConfig()
