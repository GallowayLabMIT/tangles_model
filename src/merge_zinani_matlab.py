#!/usr/bin/env python3
import pathlib
import re

import numpy as np
import pandas as pd

if __name__ == '__main__':
    dfs = []
    for f in pathlib.Path('../output/modeling_paper/zinani').glob(f'Gene*_Run*.csv'):
        match = re.match(r'Gene(?P<condition>\w+)_Run(?P<idx>\d+)\.csv', f.name)
        df = pd.DataFrame(np.genfromtxt(f, delimiter=',', skip_header=1), columns=['time_min', 'her1', 'her7'])
        df['condition'] = match['condition']
        df['idx'] = match['idx']
        dfs.append(df)
    pd.concat(dfs, ignore_index=True).to_csv('../output/modeling_paper/zinani_matlab.csv')
