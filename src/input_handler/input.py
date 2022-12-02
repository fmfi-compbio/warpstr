import argparse
import os
import sys
import yaml


def load_args(script_path: str):
    """
    Parses input args and returns parsed config
    """

    try:
        parser = argparse.ArgumentParser(description="WarpSTR")
        parser.add_argument('config', metavar='Config', type=is_valid_file, help='config file in .yaml')
        args = parser.parse_args()
    except argparse.ArgumentTypeError:
        print("Error when parsing input")
        sys.exit(-1)

    config = load_yaml(args.config)
    if config is None:
        raise ValueError(f'Error when loading config file from {args.config}')
    
    defaults_path = os.path.join(script_path,'src','default.yaml')
    defaults = load_yaml(os.path.join(script_path,'src','default.yaml'))
    if defaults is None:
        raise ValueError(f'Error when loading config file from {defaults_path}')

    add_defaults(config, defaults)
    validate_config(config)

    return config, args.config


def validate_config(config):
    """
    Validates input config
    :return: arg - path to config file
    """

    tr_config = config['tr_calling_config']
    if config['tr_region_calling']:
        if 'normalization' in tr_config and tr_config['normalization'] != 'MAD':
            raise ValueError(
                'Currently only "MAD" is supported for normalization')
        if 'min_values_per_state' in tr_config and not (isinstance(tr_config['min_values_per_state'], int) and tr_config['min_values_per_state'] > 0):
            raise ValueError(
                'Value of "min_values_per_state" must be integer and larger than zero')
        if 'states_in_segment' in tr_config and not (isinstance(tr_config['states_in_segment'], int) and tr_config['states_in_segment'] > 0):
            raise ValueError(
                'Value of "states_in_segment" must be integer and larger than zero')
        if 'min_state_similarity' in tr_config and not (isinstance(tr_config['min_state_similarity'], (float, int)) and tr_config['min_state_similarity'] > 0):
            raise ValueError(
                'Value of "min_state_similarity" must be float and larger than zero')
        if 'spike_removal' in tr_config and tr_config['spike_removal'] not in ['None', 'median3', 'median5', 'Brute']:
            raise ValueError(
                'Supported values for "spike_removal" are "None", "median3", "median5", "Brute"')
        if 'normalization' in tr_config and tr_config['normalization'] not in ['MAD', 'Without']:
            raise ValueError(
                'Supported values for "normalization" are "MAD", "Without"')
        if 'rescaling' not in tr_config:
            raise ValueError('Missing config for rescaling')
        rc_conf = tr_config['rescaling']
        if 'threshold' in rc_conf and not (isinstance(rc_conf['threshold'], (float, int)) and rc_conf['threshold'] > 0):
            raise ValueError(
                'Value of "threshold" in rescaling must be float and larger than zero')
        if 'max_std' in rc_conf and not (isinstance(rc_conf['max_std'], (float, int)) and rc_conf['max_std'] > 0):
            raise ValueError(
                'Value of "max_std" in rescaling must be float and larger than zero')
        if 'method' in rc_conf and rc_conf['method'] not in ['mean', 'median']:
            raise ValueError(
                'Supported values for "method" in rescaling are "mean", "median"')

    gt_config = config['genotyping_config']
    if config['genotyping']:
        if 'min_weight' in gt_config and not (isinstance(gt_config['min_weight'], (float, int)) and gt_config['min_weight'] >= 0 and gt_config['min_weight'] <= 1):
            raise ValueError(
                'Value of "min_weight" in rescaling must be float and larger than zero')
        if 'std_filter' in gt_config and not (isinstance(gt_config['std_filter'], (float, int)) and gt_config['std_filter'] > 1):
            raise ValueError(
                'Value of "std_filter" in rescaling must be float and larger or equal than 1')


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
            if config[key] is None or config[key].strip() == "":
                raise KeyError('Required argument is empty - %s' % key)
        else:
            # add it
            if key not in config:
                config[key] = default[key]

def is_valid_file(arg):
    """
    Checks if input config file exists
    :return: arg - path to config file
    """

    if not os.path.exists(arg):
        err_msg = (f"The file {arg} does not exist!")
        raise argparse.ArgumentTypeError(err_msg)
    return arg


def load_yaml(yaml_path):
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
        print(f'Incorrect YAML format in config path {yaml_path}')
        raise
    except IOError as err:
        print(f'IOError when processing config path {yaml_path}')
        raise

    return loaded_yaml
