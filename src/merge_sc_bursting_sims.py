#!/usr/bin/env python3
import h5py

import pathlib
from datetime import datetime, timedelta

if __name__ == '__main__':
    i = 0
    with h5py.File(f'../output/modeling_paper/modeling_sc_bursting_combined.h5','w') as outfile:
        for f in pathlib.Path('../output/modeling_paper').glob(f'fig_sc_bursting*.h5'):
            outgroup = outfile.create_group(f.stem)
            with h5py.File(f, 'r') as infile:
                for group in infile.keys():
                    if group.startswith('tangles_sc_event_run'):
                        infile.copy(infile[group], outgroup, f'tangles_sc_event_run.{i:08}')
                        i += 1
            print(f'Done with {f.stem}')
