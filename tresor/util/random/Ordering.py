__version__ = "v1.0"
__copyright__ = "Copyright 2024"
__license__ = "MIT"
__author__ = "Adam Cribbs lab"

import time
import numpy as np
import pandas as pd
from functools import wraps
from tresor.util.Console import Console


class Ordering:

    def __init__(
            self,
            method,
    ):
        self.method = method
        self.console = Console()

    def __call__(self, deal):
        if self.method == 'shuffle':
            order = self.shuffle
        elif self.method == 'permute':
            order = self.permute
        else:
            order = self.shuffle
        @wraps(deal)
        def switch(dself, *args, **kwargs):
            # print(args)
            # print(kwargs)
            if 'verbose' in kwargs['params'].keys():
                self.console.verbose = kwargs['params']['verbose']
            else:
                self.console.verbose = True
            self.console.print('======>ordering...')
            res2p = deal(dself, **kwargs)
            # print(res2p['data'].shape[0])
            res2p['data'] = order(
                data=res2p['data'],
                # use_seed=kwargs['data'],
                # seed=kwargs['data'],
            )
            # print(res2p['data'])
            return res2p
        return switch

    def shuffle(self, data, use_seed=True, seed=1):
        num_samples = len(data)
        self.console.print('=========>The number of objects to be shuffling: {}'.format(num_samples))
        if isinstance(data, pd.DataFrame):
            if use_seed:
                state = np.random.RandomState(seed)
                data = data.iloc[state.shuffle(num_samples)].reset_index(drop=True)
                self.console.print(data)
            else:
                data = data.iloc[np.random.shuffle(num_samples)].reset_index(drop=True)
        elif type(data) is np.ndarray:
            if use_seed:
                self.console.print('=========>start shuffling...')
                stime = time.time()
                state = np.random.RandomState(seed)
                state.shuffle(data)
                self.console.print('=========>shuffling time: {}'.format(time.time() - stime))
            else:
                np.random.shuffle(data)
        else:
            if use_seed:
                state = np.random.RandomState(seed)
                ids = state.shuffle(data)
                data = [data[i - 1] for i in ids]
            else:
                ids = np.random.shuffle(data)
                data = [data[i - 1] for i in ids]
        return data

    def permute(self, ):
        pass