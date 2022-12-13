import argparse
import csv
import os

import h5py
import yaml


def copy_fast5(old_read, new_path):
    with h5py.File(new_path, 'w') as new_single_read:
        new_single_read.attrs['file_version'] = 2.0
        for group in old_read:
            if group == 'Raw':
                read_number = old_read['Raw'].attrs['read_number']
                new_single_read.copy(
                    old_read[group], 'Raw/Reads/Read_{}'.format(read_number))
            elif group in ('channel_id', 'context_tags', 'tracking_id'):
                if 'UniqueGlobalKey' not in new_single_read:
                    new_single_read.create_group('UniqueGlobalKey')
                new_single_read.copy(
                    old_read[group], 'UniqueGlobalKey/{}'.format(group))
            else:
                new_single_read.copy(old_read[group], group)


parser = argparse.ArgumentParser(description='WarpSTR')
parser.add_argument('--config', help='input config', required=True)
parser.add_argument('--file', help='csv file', required=True)
args = parser.parse_args()

if not os.path.exists(args.config):
    raise FileNotFoundError(f'Config File={args.config} does not exists')

if not os.path.exists(args.file):
    raise FileNotFoundError(f'CSV File={args.file} does not exists')

with open(args.config, 'r') as stream:
    config = yaml.load(stream, Loader=yaml.FullLoader)

header = []
rows = {}
add_run_id = False
with open(args.file, 'r') as file:
    csvreader = csv.reader(file)
    header = next(csvreader)

    required_cols = ['fast5_path', 'locus', 'read_name', 'reverse', 'l_start_raw', 'r_end_raw']
    if (any(col not in required_cols for col in header)):
        raise ValueError(f'Not all required columns present in input CSV file. Required fields are: {required_cols}')

    locus_idx = header.index('locus')
    fast5_idx = header.index('fast5_path')
    name_idx = header.index('read_name')

    if 'run_id' not in header:
        header.append('run_id')
        add_run_id = True
    else:
        run_idx = header.index('run_id')

    for row in csvreader:
        locus_name = row[locus_idx]

        if add_run_id:
            run_id = 'run_0'
        else:
            run_id = row[run_idx]
        dest = os.path.join(
            config['output'],
            locus_name,
            'fast5',
            run_id,
            'annot'
        )
        if not os.path.exists(dest):
            os.makedirs(dest)

        readname = row[name_idx]
        fast5file = h5py.File(row[fast5_idx], 'r')
        readfullname = f'read_{readname}'

        if readfullname in fast5file:
            new_path = os.path.join(dest, f'{readfullname[5:]}.fast5'.format())
            copy_fast5(fast5file[readfullname], new_path)
        else:
            raise ValueError(f'Read {readname} not found in fast5 file {row[fast5_idx]}')

        if locus_name not in rows:
            rows[locus_name] = []
        rows[locus_name].append(row)

for locus in rows:
    overview_file = os.path.join(
        config['output'],
        locus,
        'overview.csv'
    )
    header.append('saved')
    with open(overview_file, 'w') as file:
        csvwriter = csv.writer(file, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
        csvwriter.writerow(header)
        for row in rows[locus]:
            if add_run_id:
                row.append('run_0')
            row.append('1')
            csvwriter.writerow(row)

print(f'Finished preparing WarpSTR document structure for the input config {args.config} and csv file {args.file}')
print('Continue with calling:')
print(f'>>>python WarpSTR.py {args.config}')
print('where steps "single_read_extraction", "guppy_annotation" and "tr_region_extraction" are set to "False"')
