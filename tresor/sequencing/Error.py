__version__ = "0.0.1"
__copyright__ = "Copyright 2025"
__license__ = "GPL-3.0"
__developer__ = "Jianfeng Sun"
__maintainer__ = "Jianfeng Sun"
__email__ = "jianfeng.sunmt@gmail.com"


import time
import numpy as np
import pandas as pd
from tresor.util.random.Number import Number as rannum
from tresor.read.error.Library import Library as errlib
from tresor.util.sequence.symbol.Single import Single as dnasgl
from tresor.util.Console import Console


class Error:

    def __init__(
            self,
            method,
    ):
        self.method = method
        self.errlib = errlib()

        self.console = Console()

    def __call__(self, deal):
        from functools import wraps
        # if self.method == 'default':
        #     func = self.postable
        if self.method == 'default':
            func = self.table2D
        @wraps(deal)
        def indexing(ph, *args, **kwargs):
            res2p = deal(ph, **kwargs)
            # print(func(res2p))
            # print(kwargs['params'].keys())

            if 'verbose' in kwargs['params'].keys():
                self.console.verbose = kwargs['params']['verbose']
            else:
                self.console.verbose = True
            self.console.print('======>Generate sequencing errors...')

            return func(res2p)
        return indexing

    def postable(self, res2p, kind='index_by_lambda'):
        # print(res2p.keys())
        seq_stime = time.time()
        data_seq = pd.DataFrame(res2p['data_spl'], columns=['read', 'sam_id', 'source'])
        del res2p['data']
        del res2p['data_spl']
        self.console.print('=========>There are {} reads to be sequenced'.format(data_seq.shape[0]))
        self.console.print('=========>The position table construction starts')
        pcr_postable_stime = time.time()
        # print(data_seq['read'][0])
        if kind == 'index_by_same_len':
            seq_pos_ids, seq_ids = self.postableIndexBySameLen(
                seq_len=len(data_seq['read'][0]),
                num_seq=data_seq.shape[0],
            )
        elif kind == 'index_by_lambda':
            seq_ids = []
            seq_pos_ids = []
            data_seq.apply(lambda x: self.postableLambda(x, seq_ids, seq_pos_ids), axis=1)
        else:
            seq_pos_ids, seq_ids = self.postableIndexBySameLen(
                seq_len=len(data_seq['read'][0]),
                num_seq=data_seq.shape[0],
            )
        pos_table = {'seq_ids': seq_ids, 'seq_pos_ids': seq_pos_ids}
        print(pos_table)
        # print(len(seq_ids))
        pcr_postable_etime = time.time()
        self.console.print('=========>Time for constructing the position table: {time:.3f}s'.format(time=pcr_postable_etime - pcr_postable_stime))
        seq_nt_num = len(seq_ids)
        self.console.print('=========>There are {} nucleotides to be sequenced'.format(seq_nt_num))
        self.console.print('=========>Determination of number of sequencing errors start')
        seq_err_num_simu_stime = time.time()
        if res2p['err_num_met'] == 'binomial':
            seq_err_num = rannum().binomial(n=seq_nt_num, p=res2p['seq_error'], use_seed=True, seed=1)
        elif res2p['err_num_met'] == 'nbinomial':
            seq_err_num = rannum().nbinomial(
                n=seq_nt_num * (1 - res2p['seq_error']),
                p=1 - res2p['seq_error'],
                use_seed=True,
                seed=1,
            )
            print(seq_err_num)
        else:
            seq_err_num = rannum().binomial(n=seq_nt_num, p=res2p['seq_error'], use_seed=True, seed=1)

        self.console.print('============>There are {} nucleotide errors during sequencing'.format(seq_err_num))
        err_lin_ids = rannum().uniform(low=0, high=seq_nt_num, num=seq_err_num, use_seed=True, seed=1)
        # print(err_lin_ids)
        arr_err_pos = []# [[row1, col1], [row2, col2], ...]
        for i in err_lin_ids:
            arr_err_pos.append([pos_table['seq_ids'][i], pos_table['seq_pos_ids'][i]])
        pseudo_nums = rannum().uniform(low=0, high=3, num=seq_err_num, use_seed=False)
        # print(pseudo_nums)
        seq_err_num_simu_etime = time.time()
        self.console.print('=========>Time for determining sequencing errors: {time:.3f}s'.format(time=seq_err_num_simu_etime - seq_err_num_simu_stime))
        self.console.print('=========>Sequencing error Assignment starts')
        seq_err_assign_stime = time.time()
        data_seq['read'] = data_seq.apply(lambda x: list(x['read']), axis=1)
        for pos_err, pseudo_num in zip(arr_err_pos, pseudo_nums):
            # print(pos_err[0])
            # print(pos_err[1])
            pcr_err_base = data_seq.loc[pos_err[0], 'read'][pos_err[1]]
            dna_map = dnasgl().todict(
                nucleotides=dnasgl().getEleTrimmed(
                    ele_loo=pcr_err_base,
                    universal=True,
                ),
                reverse=True,
            )
            # print('before', data_seq.loc[pos_err[0], 'read'][pos_err[1]])
            data_seq.loc[pos_err[0], 'read'][pos_err[1]] = dna_map[pseudo_num]
            # print('after', data_seq.loc[pos_err[0], 'read'][pos_err[1]])
        del arr_err_pos
        del pseudo_nums
        seq_err_assign_etime = time.time()
        self.console.print('=========>Time for assigning sequencing errors: {time:.2f}s'.format(time=seq_err_assign_etime - seq_err_assign_stime))
        data_seq['read'] = data_seq.apply(lambda x: ''.join(x['read']), axis=1)
        if res2p['seq_deletion']:
            data_seq['read'] = data_seq['read'].apply(lambda x: self.deletion(read=x, del_rate=res2p['seq_del_rate']))
        if res2p['seq_insertion']:
            data_seq['read'] = data_seq['read'].apply(lambda x: self.insertion(read=x, ins_rate=res2p['seq_ins_rate']))
        res2p['data'] = data_seq.values
        del data_seq
        seq_etime = time.time()
        self.console.print('=========>Sequencing time: {time:.2f}s'.format(time=seq_etime - seq_stime))
        return res2p

    def table2D(self, res2p):
        seq_stime = time.time()
        data_seq = pd.DataFrame(res2p['data_spl'], columns=['read', 'sam_id', 'source'])
        del res2p['data']
        del res2p['data_spl']
        self.console.print('=========>There are {} reads to be sequenced'.format(data_seq.shape[0]))

        data_seq['read_len'] = data_seq['read'].apply(lambda x: len(x))
        if res2p['err_num_met'] == 'binomial':
            data_seq['num_err_per_read'] = data_seq['read_len'].apply(lambda x: rannum().binomial(
                n=x, p=res2p['seq_error'], use_seed=False, seed=1
            ))
            res2p['num_err_per_read_dict']['sequencing'] = data_seq['num_err_per_read']
            if res2p['vis_num_err_per_read']:
                from tresor.pcr.Binomial import validate_binomial_and_plot
                results = validate_binomial_and_plot(
                    data_seq['num_err_per_read'],
                    n_trials=None, p_true=None,
                    title="sequencing per-read error counts (Binomial)"
                )
                for k, v in results.items():
                    print(
                        f"{v.model}: n={v.n}, params={v.params}, loglik={v.loglik:.2f}, "
                        f"AIC={v.aic:.2f}, BIC={v.bic:.2f}, chi2={v.chi2:.2f}, df={v.df}, p={v.p_value:.3g}, bins={v.n_bins_used}. {v.notes}"
                    )
        elif res2p['err_num_met'] == 'nbinomial':
            data_seq['num_err_per_read'] = data_seq['read_len'].apply(lambda x: rannum().nbinomial(
                n=x*(1-res2p['seq_error']), p=1-res2p['seq_error'], use_seed=False, seed=1
            ))
            res2p['num_err_per_read_dict']['sequencing'] = data_seq['num_err_per_read']
            if res2p['vis_num_err_per_read']:
                from tresor.pcr.NBinomial import validate_nbinom_and_plot
                results = validate_nbinom_and_plot(
                    data_seq['num_err_per_read'],
                    fit_poisson=False,
                    title="sequencing per-read error counts (NBinomial)"
                )
                for k, v in results.items():
                    print(
                        f"{v.model}: {v.params}, loglik={v.loglik:.2f}, AIC={v.aic:.2f}, BIC={v.bic:.2f}, "
                        f"chi2={v.chi2:.2f}, df={v.df}, p={v.p_value:.3g}, bins={v.n_bins_used}. {v.notes}"
                    )

        data_seq['pos_err_per_read'] = data_seq.apply(lambda x: rannum().uniform(
            low=0, high=x['read_len'], num=x['num_err_per_read'], use_seed=False, seed=1
        ), axis=1)
        data_seq['base_roll_per_read'] = data_seq['num_err_per_read'].apply(lambda x: rannum().uniform(
            low=0, high=3, num=x, use_seed=False
        ))
        data_seq['read'] = data_seq.apply(lambda x: self.errlib.change(
            read=x['read'],
            pos_list=x['pos_err_per_read'],
            base_list=x['base_roll_per_read'],
        ), axis=1)
        # print(data_seq)
        # print(data_seq[['read', 'read_pcr', 'pos_err_per_read']])

        data_seq['read'] = data_seq.apply(lambda x: ''.join(x['read']), axis=1)
        if res2p['seq_deletion']:
            data_seq['read'] = data_seq['read'].apply(lambda x: self.deletion(read=x, del_rate=res2p['seq_del_rate']))
        if res2p['seq_insertion']:
            data_seq['read'] = data_seq['read'].apply(lambda x: self.insertion(read=x, ins_rate=res2p['seq_ins_rate']))

        data_seq = data_seq[['read', 'sam_id', 'source']]
        res2p['data'] = data_seq.values
        del data_seq
        seq_etime = time.time()
        self.console.print('=========>Sequencing time: {time:.2f}s'.format(time=seq_etime - seq_stime))
        return res2p

    def deletion(self, read, del_rate):
        num_err_per_read = rannum().binomial(
            n=len(read), p=del_rate, use_seed=False, seed=False
        )
        pos_list = rannum().choice(
            high=len(read), num=num_err_per_read, use_seed=False, seed=False, replace=False,
        )
        for _, pos in enumerate(pos_list):
            read = read[:pos] + read[pos + 1:]
        return read

    def insertion(self, read, ins_rate):
        num_err_per_read = rannum().binomial(
            n=len(read), p=ins_rate, use_seed=False, seed=False
        )
        pos_list = rannum().choice(
            high=len(read), num=num_err_per_read, use_seed=False, seed=False, replace=False,
        )
        base_list = rannum().uniform(
            low=0, high=4, num=num_err_per_read, use_seed=False
        )
        for i, pos in enumerate(pos_list):
            dna_map = dnasgl().todict(
                nucleotides=dnasgl().get(
                    universal=True,
                ),
                reverse=True,
            )
            read = read[:pos] + dna_map[base_list[i]] + read[pos:]
            ### read
            ### pos, base_list[i], dna_map[base_list[i]], dna_map
            ### read
            # 5 3 G {0: 'A', 1: 'T', 2: 'C', 3: 'G'}
            # TTTTTTTTTGGGCCCGGGAAAAAACCCAAAGGGGGG
            # TTTTTGTTTTGGGCCCGGGAAAAAACCCAAAGGGGGG
            # 9 0 A {0: 'A', 1: 'T', 2: 'C', 3: 'G'}
            # CCCTTTCCCTTTGGGTTTGGGTTTCCCGGGAAACCC
            # CCCTTTCCCATTTGGGTTTGGGTTTCCCGGGAAACCC
            # 3 0 A {0: 'A', 1: 'T', 2: 'C', 3: 'G'}
            # AAATTTTTTAAACCCAAAAAAAAAAAATTTTTTCCC
            # AAAATTTTTTAAACCCAAAAAAAAAAAATTTTTTCCC
        return read

    def postableIndexBySameLen(self, seq_len, num_seq):
        nt_ids = [i for i in range(seq_len)]
        seq_pos_ids = nt_ids * num_seq
        seq_ids = np.array([[i] * seq_len for i in range(num_seq)]).ravel().tolist()
        return seq_pos_ids, seq_ids

    def postableLambda(self, x, ids, pos_ids):
        l = len(x['read'])
        for i in range(l):
            pos_ids.append(i)
        for i in [x.name] * l:
            ids.append(i)
        return ids, pos_ids

    def todict(self, bases, reverse=False):
        aa_dict = {}
        for k, v in enumerate(bases):
            aa_dict[v] = k
        if reverse:
            aa_dict = {v: k for k, v in aa_dict.items()}
        return aa_dict