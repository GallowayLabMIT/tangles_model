#!/usr/bin/env python3
import h5py

import pathlib
from datetime import datetime, timedelta

if __name__ == '__main__':
    i = 0
    with h5py.File('../output/modeling_paper/energy_torque_perturbations.h5','w') as outfile:
        for f in pathlib.Path('../output/modeling_paper').glob('energy_torque_perturbations-node*.h5'):
            outgroup = outfile.create_group(f.stem)
            with h5py.File(f, 'r') as infile:
                for group in infile.keys():
                    if group.startswith('tangles_summarized_run.'):
                        infile.copy(infile[group], outgroup, f'tangles_summarized_run.{i:08}')
                        i += 1
                    if group.startswith('tangles_sc_event_run.'):
                        infile.copy(infile[group], outgroup, f'tangles_sc_event_run.{i:08}')
                        i += 1
            print(f'Done with {f.stem}')
