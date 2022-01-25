#!/usr/bin/env python3
import h5py

import pathlib
from datetime import datetime, timedelta

if __name__ == '__main__':
    i = 0
    for fig in [2,4]:
        with h5py.File(f'../output/modeling_paper/modeling_fig{fig}_combined.h5','w') as outfile:
            for f in pathlib.Path('../output/modeling_paper').glob(f'fig{fig}*.h5'):
                # Skip HDF files written within the last day
                if datetime.now() - datetime.fromtimestamp(f.stat().st_mtime) < timedelta(days=1):
                    continue
                outgroup = outfile.create_group(f.stem)
                with h5py.File(f, 'r') as infile:
                    for group in infile.keys():
                        if group.startswith('tangles_discrete_run'):
                            infile.copy(infile[group], outgroup, f'tangles_discrete_run.{i:08}')
                            i += 1
                print(f'Done with {f.stem}')
