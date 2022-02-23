#!/usr/bin/env python3
import h5py

import pathlib
from datetime import datetime, timedelta

if __name__ == '__main__':
    i = 0
    with h5py.File(f'../output/modeling_paper/modeling_toggles_combined.h5','w') as outfile:
        for f in pathlib.Path('../output/modeling_paper').glob(f'fig_toggles_sims*.h5'):
            outgroup = outfile.create_group(f.stem)
            with h5py.File(f, 'r') as infile:
                for group in infile.keys():
                    if group.startswith('tangles_discrete_run'):
                        if i % 5 == 0:
                            infile.copy(infile[group], outgroup, f'tangles_discrete_run.{i:08}')
                        i += 1
            print(f'Done with {f.stem}')
    with h5py.File(f'../output/modeling_paper/modeling_toggles_topo_combined.h5','w') as outfile:
        for f in pathlib.Path('../output/modeling_paper').glob(f'fig_toggles_topo_sims*.h5'):
            outgroup = outfile.create_group(f.stem)
            with h5py.File(f, 'r') as infile:
                for group in infile.keys():
                    if group.startswith('tangles_discrete_run'):
                        if i % 5 == 0:
                            infile.copy(infile[group], outgroup, f'tangles_discrete_run.{i:08}')
                        i += 1
            print(f'Done with {f.stem}')
