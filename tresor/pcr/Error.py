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
from tresor.util.sequence.symbol.Single import Single as dnasgl
from tresor.read.error.Library import Library as errlib
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
        @wraps(deal)
        def indexing(ph, *args, **kwargs):
            res2p = deal(ph, **kwargs)
            # print(func(res2p))

            if 'verbose' in kwargs['params'].keys():
                self.console.verbose = kwargs['params']['verbose']
            else:
                self.console.verbose = True

            self.console.print('======>Generate PCR errors...')

            if self.method == 'err1d':
                self.console.print('=========>The error 2D method is being used')
                func = self.table1D
            elif self.method == 'err2d':
                self.console.print('=========>The error 2D method is being used')
                func = self.table2D
            elif self.method == 'bftree':
                self.console.print('=========>The boolean flag tree method is being used')
                func = self.table_tree
            elif self.method == 'sptree':
                self.console.print('=========>The ShortPath method is being used')
                func = self.table_tree
            elif self.method == 'mutation_table_minimum':
                self.console.print('=========>The mutation_table_minimum method is being used')
                func = self.table_mutation_minimum
            elif self.method == 'mutation_table_complete':
                self.console.print('=========>The mutation_table_complete method is being used')
                func = self.table_mutation_complete
            else: # default
                self.console.print('=========>The ResimPy method is being used')
                func = self.table_tree
            return func(res2p)
        return indexing

    def table1D(
            self,
            res2p,
            kind='index_by_same_len',
    ):
        """
        ..  test:
            -----
            # dna = ['A', 'T', 'G', 'C']
            # dna.remove(pcr_err_base)
            # dna_map = self.todict(dna, reverse=True)
            # data_pcr.loc[pos_err[0], 'read'][pos_err[1]] = dna_map[pseudo_num]
            # pcr_err_base = data_pcr[pos_err[0]][0][pos_err[1]]

            # d = [
            #     [i, j] for i in ampl_read_len_df.index for j in np.arange(ampl_read_len_df[i])
            # ]
            # print(len(d))
            # ampl_read_len_list = np.arange(len(d))[:, np.newaxis]
            # map_tab = np.concatenate((ampl_read_len_list, d), axis=1)
            # print(self.tactic1(map_tab))

        ..  Method:
            -------
            =>string concat
                for pos_err, pseudo_num in zip(arr_err_pos, pseudo_nums):
                    tar_read = data_pcr.loc[pos_err[0], 'read']
                    pcr_err_base = tar_read[pos_err[1]]
                    dna_map = dnasgl().todict(
                        nucleotides=dnasgl().getEleTrimmed(
                            ele_loo=pcr_err_base,
                            universal=True,
                        ),
                        reverse=True,
                    )
                    # print('before', data_pcr.loc[pos_err[0], 'read'][pos_err[1]])
                    data_pcr.loc[pos_err[0], 'read'] = tar_read[:pos_err[1]] + dna_map[pseudo_num] + tar_read[pos_err[1]+1:]
        :param res2p:
        :return:
        """
        # print(res2p.keys())
        pcr_stime = time.time()
        data_pcr = pd.DataFrame(res2p['data_spl'], columns=['read', 'sam_id', 'source'])
        del res2p['data_spl']
        # print(len(data_pcr['read'][0]))
        # print(data_pcr.shape[0])
        res2p['recorder_pcr_read_num'].append(data_pcr.shape[0])
        self.console.print('======>{} reads to be amplified'.format(data_pcr.shape[0]))
        self.console.print('======>constructing the position table starts...')
        # print(data_pcr)
        pcr_postable_stime = time.time()
        if kind == 'index_by_same_len':
            seq_pos_ids, seq_ids = self.postableIndexBySameLen(seq_len=len(data_pcr['read'][0]), num_seq=data_pcr.shape[0])
        elif kind == 'index_by_lambda':
            seq_ids = []
            seq_pos_ids = []
            data_pcr.apply(lambda x: self.postableLambda(x, seq_ids, seq_pos_ids), axis=1)
        else:
            seq_pos_ids, seq_ids = self.postableIndexBySameLen(seq_len=len(data_pcr['read'][0]), num_seq=data_pcr.shape[0])
        # print(seq_ids)
        pos_table = {'seq_ids': seq_ids, 'seq_pos_ids': seq_pos_ids}
        self.console.print('======>time for constructing the position table  {time:.3f}s'.format(time=time.time() - pcr_postable_stime))
        ampl_nt_num = len(seq_ids)
        self.console.print('======>{} nucleotides to be amplified'.format(ampl_nt_num))
        res2p['recorder_nucleotide_num'].append(ampl_nt_num)
        self.console.print('======>determining PCR error numbers starts...')
        pcr_err_num_simu_stime = time.time()
        if res2p['err_num_met'] == 'binomial':
            pcr_err_num = rannum().binomial(n=ampl_nt_num, p=res2p['pcr_error'], use_seed=True, seed=res2p['ipcr'] + 1)
        elif res2p['err_num_met'] == 'nbinomial':
            pcr_err_num = rannum().nbinomial(
                n=ampl_nt_num*(1-res2p['pcr_error']),
                p=1-res2p['pcr_error'],
                use_seed=True,
                seed=res2p['ipcr'] + 1
            )
            # print(pcr_err_num)
        else:
            pcr_err_num = rannum().binomial(n=ampl_nt_num, p=res2p['pcr_error'], use_seed=True, seed=res2p['ipcr'] + 1)
        self.console.print('======>{} nucleotides to be erroneous at this PCR'.format(pcr_err_num))
        res2p['recorder_pcr_err_num'].append(pcr_err_num)
        spl_nt_ids = rannum().uniform(low=0, high=ampl_nt_num, num=pcr_err_num, use_seed=True, seed=res2p['ipcr'] + 1)
        arr_err_pos = []
        for i in spl_nt_ids:
            arr_err_pos.append([pos_table['seq_ids'][i], pos_table['seq_pos_ids'][i]])
        pseudo_nums = rannum().uniform(low=0, high=3, num=pcr_err_num, use_seed=False)
        # print(pseudo_nums)
        self.console.print('======>time for determining PCR error numbers  {time:.3f}s'.format(time=time.time() - pcr_err_num_simu_stime))
        self.console.print('======>assigning PCR errors starts...')
        pcr_err_assign_stime = time.time()
        data_pcr['read'] = data_pcr.apply(lambda x: list(x['read']), axis=1)
        # caser = self.todict5(arr_err_pos)
        # data_pcr['read'] = data_pcr.apply(lambda x: self.worktable(x, caser), axis=1)
        for pos_err, pseudo_num in zip(arr_err_pos, pseudo_nums):
            pcr_err_base = data_pcr.loc[pos_err[0], 'read'][pos_err[1]]
            dna_map = dnasgl().todict(
                nucleotides=dnasgl().getEleTrimmed(
                    ele_loo=pcr_err_base,
                    universal=True,
                ),
                reverse=True,
            )
            # print('before', data_pcr.loc[pos_err[0], 'read'][pos_err[1]])
            data_pcr.loc[pos_err[0], 'read'][pos_err[1]] = dna_map[pseudo_num]
            # print('after', data_pcr.loc[pos_err[0], 'read'][pos_err[1]])
        self.console.print('======>time for assigning PCR errors {time:.2f}s'.format(time=time.time() - pcr_err_assign_stime))
        self.console.print('======>merging the PCR duplicates and all previous sequences starts...')
        del arr_err_pos
        del spl_nt_ids
        pcr_merge_stime = time.time()
        data_pcr['read'] = data_pcr.apply(lambda x: ''.join(x['read']), axis=1)
        data_pcr['source'] = 'pcr-' + str(res2p['ipcr'] + 1)
        # print(data_pcr)
        # print(res2p)
        data_pcr['read'] = data_pcr['read'].apply(lambda x: self.errlib.insertion(read=x, ins_rate=res2p['pcr_ins_rate']))
        data_pcr['read'] = data_pcr['read'].apply(lambda x: self.errlib.deletion(read=x, del_rate=res2p['pcr_del_rate']))

        data_pcr = np.array(data_pcr)
        # print(data_pcr)
        # print(data_pcr.shape)

        res2p['data'] = np.concatenate((res2p['data'], data_pcr), axis=0)
        del data_pcr
        # print(res2p['recorder_pcr_err_num'])
        self.console.print('======>time for merging sequences {time:.2f}s'.format(time=time.time() - pcr_merge_stime))
        self.console.print('======>Summary report:')
        self.console.print('=========>PCR time: {time:.2f}s'.format(time=time.time() - pcr_stime))
        self.console.print('=========>the dimensions of the data: number of reads: {}'.format(res2p['data'].shape))
        self.console.print('=========>the number of reads at this PCR: {}, '.format(res2p['recorder_pcr_read_num']))
        self.console.print('=========>the number of nucleotides at this PCR: {}, '.format(res2p['recorder_nucleotide_num']))
        self.console.print('=========>the number of errors at this PCR: {}, '.format(res2p['recorder_pcr_err_num']))
        return res2p

    def table2D(self, res2p):
        pcr_stime = time.time()
        data_pcr = pd.DataFrame(res2p['data_spl'], columns=['read', 'sam_id', 'source'])
        data_pcr['read_len'] = data_pcr['read'].apply(lambda x: len(x))
        if res2p['err_num_met'] == 'binomial':
            data_pcr['num_err_per_read'] = data_pcr['read_len'].apply(lambda x: rannum().binomial(
                n=x, p=res2p['pcr_error'], use_seed=False, seed=res2p['ipcr'] + 1
            ))
            if 'num_err_per_read_dict' in res2p.keys():
                res2p['num_err_per_read_dict'][res2p['ipcr'] + 1] = data_pcr['num_err_per_read']
            if 'vis_num_err_per_read' in res2p.keys():
                if res2p['vis_num_err_per_read']:
                    from tresor.vis.Binomial import validate_binomial_and_plot
                    results = validate_binomial_and_plot(
                        data_pcr['num_err_per_read'],
                        n_trials=None, p_true=None,
                        title="PCR per-read error counts (Binomial)"
                    )
                    for k, v in results.items():
                        print(
                            f"{v.model}: n={v.n}, params={v.params}, loglik={v.loglik:.2f}, "
                            f"AIC={v.aic:.2f}, BIC={v.bic:.2f}, chi2={v.chi2:.2f}, df={v.df}, p={v.p_value:.3g}, bins={v.n_bins_used}. {v.notes}"
                        )
        elif res2p['err_num_met'] == 'nbinomial':
            data_pcr['num_err_per_read'] = data_pcr['read_len'].apply(lambda x: rannum().nbinomial(
                n=x*(1-res2p['pcr_error']), p=1-res2p['pcr_error'], use_seed=False, seed=res2p['ipcr'] + 1
            ))
            if 'num_err_per_read_dict' in res2p.keys():
                res2p['num_err_per_read_dict'][res2p['ipcr'] + 1] = data_pcr['num_err_per_read']
            if 'vis_num_err_per_read' in res2p.keys():
                if res2p['vis_num_err_per_read']:
                    from tresor.vis.NBinomial import validate_nbinom_and_plot
                    results = validate_nbinom_and_plot(
                        data_pcr['num_err_per_read'],
                        fit_poisson=False,
                        title="PCR per-read error counts (NBinomial)"
                    )
                    for k, v in results.items():
                        print(
                            f"{v.model}: {v.params}, loglik={v.loglik:.2f}, AIC={v.aic:.2f}, BIC={v.bic:.2f}, "
                            f"chi2={v.chi2:.2f}, df={v.df}, p={v.p_value:.3g}, bins={v.n_bins_used}. {v.notes}"
                        )
        elif res2p['err_num_met'] == 'betabinomial':
            data_pcr['num_err_per_read'] = data_pcr['read_len'].apply(lambda x: rannum().beta_binomial(
                n=x, p=res2p['pcr_error'], use_seed=False, seed=1
            ))
            if 'num_err_per_read_dict' in res2p.keys():
                res2p['num_err_per_read_dict'][res2p['ipcr'] + 1] = data_pcr['num_err_per_read']
            if 'vis_num_err_per_read' in res2p.keys():
                if res2p['vis_num_err_per_read']:
                    from tresor.vis.BetaBinomial import validate_nbinom_and_plot
                    results = validate_nbinom_and_plot(
                        data_pcr['num_err_per_read'],
                        fit_poisson=False,
                        title="PCR per-read error counts (BetaBinomial)"
                    )
                    for k, v in results.items():
                        print(
                            f"{v.model}: {v.params}, loglik={v.loglik:.2f}, AIC={v.aic:.2f}, BIC={v.bic:.2f}, "
                            f"chi2={v.chi2:.2f}, df={v.df}, p={v.p_value:.3g}, bins={v.n_bins_used}. {v.notes}"
                        )
        elif res2p['err_num_met'] == 'poisson':
            data_pcr['num_err_per_read'] = data_pcr['read_len'].apply(lambda x: rannum().poisson(
                n=x, p=res2p['pcr_error'], use_seed=False, seed=1
            ))
            if 'num_err_per_read_dict' in res2p.keys():
                res2p['num_err_per_read_dict'][res2p['ipcr'] + 1] = data_pcr['num_err_per_read']
            if 'vis_num_err_per_read' in res2p.keys():
                if res2p['vis_num_err_per_read']:
                    from tresor.vis.Poisson import validate_nbinom_and_plot
                    results = validate_nbinom_and_plot(
                        data_pcr['num_err_per_read'],
                        fit_poisson=False,
                        title="PCR per-read error counts (Poisson)"
                    )
                    for k, v in results.items():
                        print(
                            f"{v.model}: {v.params}, loglik={v.loglik:.2f}, AIC={v.aic:.2f}, BIC={v.bic:.2f}, "
                            f"chi2={v.chi2:.2f}, df={v.df}, p={v.p_value:.3g}, bins={v.n_bins_used}. {v.notes}"
                        )

        data_pcr['pos_err_per_read'] = data_pcr.apply(lambda x: rannum().uniform(
            low=0, high=x['read_len'], num=x['num_err_per_read'], use_seed=False, seed=res2p['ipcr'] + 1
        ), axis=1)
        data_pcr['base_roll_per_read'] = data_pcr['num_err_per_read'].apply(lambda x: rannum().uniform(
            low=0, high=3, num=x, use_seed=False
        ))
        data_pcr['read_pcr'] = data_pcr.apply(lambda x: self.errlib.change(
            read=x['read'],
            pos_list=x['pos_err_per_read'],
            base_list=x['base_roll_per_read'],
        ), axis=1)
        # print(data_pcr[['read', 'read_pcr', 'pos_err_per_read']])
        del res2p['data_spl']
        pcr_merge_stime = time.time()
        data_pcr['sam_id'] = data_pcr['sam_id'].apply(lambda x: x + '_' + str(res2p['ipcr'] + 1))
        data_pcr['source'] = 'pcr-' + str(res2p['ipcr'] + 1)
        ### data_pcr.columns
        # ['read', 'sam_id', 'source', 'read_len', 'num_err_per_read',
        #         'pos_err_per_read', 'base_roll_per_read', 'read_pcr']
        ### data_pcr
        # read  ...                                           read_pcr
        # 0     CGCGTTAGTAATTCATTTTTTTCCCTTTAAAGGGTTTAAATTTGGG...  ...  CGCGTTAGTAATTCATTTTTTTCCCTTTAAAGGGTTTAAATTTGGG...
        # 1     CGCGTTAGTAATTCATCCCTTTAAAAAATTTGGGTTTGGGAAAAAA...  ...  CGCGTTAGTAATTCATCCCTTTAAAAAATTTGGGTTTGGGAAAAAA...
        # 2     CGCGTTAGTAATTCATGGGCCCTTTTTTAAATTTAAACCCCCCAAA...  ...  CGCGTTAGTAATTCATGGGCCCTTTTTTAAATTTAAACCCCCCAAA...
        # ...                                                 ...  ...                                                ...
        # 8692  GAAATCATGTAGTTCGCCCCCCAAATTTAAACCCAAATTTCCCAAA...  ...  GAAATCATGTAGTTCGCCCCCCAAATTTAAACCCAAATTTCCCAAA...
        # 8693  CGCGTTAGTAATTCATTTTGGGGGGAAACCCAAACCCCCCCCCGGG...  ...  CGCGTTAGTAATTCATTTTGGGGGGAAACCCAAACCCCCCCCCGGG...
        # [8694 rows x 8 columns]
        # print(data_pcr)
        # print(res2p)
        data_pcr['read'] = data_pcr['read'].apply(
            lambda x: self.errlib.insertion(read=x, ins_rate=res2p['pcr_ins_rate']))
        data_pcr['read'] = data_pcr['read'].apply(
            lambda x: self.errlib.deletion(read=x, del_rate=res2p['pcr_del_rate']))

        data_pcr = np.array(data_pcr[['read_pcr', 'sam_id', 'source']])
        res2p['data'] = np.concatenate((res2p['data'], data_pcr), axis=0)
        del data_pcr
        # print(res2p)
        # print(res2p['recorder_pcr_read_num'])
        self.console.print('======>time for merging sequences {time:.2f}s'.format(time=time.time() - pcr_merge_stime))
        self.console.print('======>Summary report:')
        self.console.print('=========>PCR time: {time:.2f}s'.format(time=time.time() - pcr_stime))
        self.console.print('=========>the dimensions of the data: number of reads: {}'.format(res2p['data'].shape))
        self.console.print('=========>the number of reads at this PCR: {}, '.format(res2p['recorder_pcr_read_num']))
        self.console.print('=========>the number of nucleotides at this PCR: {}, '.format(res2p['recorder_nucleotide_num']))
        self.console.print('=========>the number of errors at this PCR: {}, '.format(res2p['recorder_pcr_err_num']))
        return res2p

    def table_mutation_minimum(self, res2p):
        """

        Notes
        -----

        Parameters
        ----------
        res2p

        Returns
        -------

        """
        pcr_stime = time.time()
        df_mut_info = pd.DataFrame()
        data_pcr = pd.DataFrame(res2p['data_spl'], columns=['read_len', 'sam_id', 'source'])
        data_pcr['num_err_per_read'] = data_pcr['read_len'].apply(lambda x: rannum().binomial(
            n=int(x), p=res2p['pcr_error'], use_seed=False, seed=res2p['ipcr'] + 1
        ))
        # data_pcr['num_err_per_read']
        # index num_err_per_read
        # 0     0
        # 1     0
        # 2     1
        # ...
        # 41    0
        # 42    0

        df_mut_info['pos_err_per_read'] = data_pcr.apply(lambda x: rannum().uniform(
            low=0, high=x['read_len'], num=x['num_err_per_read'], use_seed=False, seed=res2p['ipcr'] + 1
        ), axis=1)
        df_mut_info['base_roll_per_read'] = data_pcr['num_err_per_read'].apply(lambda x: rannum().uniform(
            low=0, high=3, num=x, use_seed=False
        ))

        # print(data_pcr[['read', 'read_pcr', 'pos_err_per_read']])
        del res2p['data_spl']
        pcr_merge_stime = time.time()
        data_pcr['sam_id'] = data_pcr['sam_id'].apply(lambda x: x + '_' + str(res2p['ipcr'] + 1))
        data_pcr['source'] = 'pcr-' + str(res2p['ipcr'] + 1)

        # print(df_mut_info)
        df_mut_info['sam_id'] = data_pcr['sam_id'].copy()
        # print(data_pcr)
        df_mut_info['mark'] = df_mut_info['pos_err_per_read'].apply(lambda x: 1 if x.size == 0 else 0)
        # print(df_mut_info)
        df_mut_info = df_mut_info.drop(df_mut_info.loc[df_mut_info['mark'] == 1].index)
        df_mut_info = df_mut_info.drop(['mark'], axis=1)
        ### df_mut_info
        # pos_err_per_read base_roll_per_read                          sam_id
        # 201                 [9]                [1]           154*c*1*g*4*_3_4_9_10
        # 264                [47]                [1]        53*c*1*g*4*_1_2_6_7_9_10
        # 870                [28]                [0]           7*c*1*g*3*_3_4_8_9_10
        # 946                [13]                [2]        13*c*1*g*0*_1_3_4_8_9_10
        # 1312               [40]                [1]      68*c*1*g*4*_3_5_6_7_8_9_10
        # ...                 ...                ...                             ...
        # 152258             [48]                [2]             115*c*1*g*4*_4_5_10
        # 152548              [6]                [1]         164*c*1*g*4*_3_4_6_7_10
        # 152656             [30]                [2]         173*c*1*g*4*_1_5_7_9_10
        # 152914             [33]                [0]  46*c*1*g*4*_1_2_3_4_5_6_7_9_10
        # 153185              [5]                [2]         177*c*1*g*4*_2_5_6_8_10
        #
        # [792 rows x 3 columns]
        data_pcr = np.array(data_pcr[['read_len', 'sam_id', 'source']])
        res2p['data'] = np.concatenate((res2p['data'], data_pcr), axis=0)
        # print(np.concatenate((res2p['mut_info'], np.array(df_mut_info)), axis=0))
        res2p['mut_info'] = np.concatenate((res2p['mut_info'], np.array(df_mut_info)), axis=0)
        del data_pcr
        self.console.print('======>Time for merging sequences {time:.2f}s'.format(time=time.time() - pcr_merge_stime))
        self.console.print('======>Summary report:')
        self.console.print('=========>PCR time: {time:.2f}s'.format(time=time.time() - pcr_stime))
        self.console.print('=========>Dimension of the data: number of reads: {}'.format(res2p['data'].shape))
        self.console.print('=========>Number of reads at the PCR: {}, '.format(res2p['recorder_pcr_read_num']))
        self.console.print('=========>Number of nucleotides at the PCR: {}, '.format(res2p['recorder_nucleotide_num']))
        self.console.print('=========>Number of errors at the PCR: {}, '.format(res2p['recorder_pcr_err_num']))
        return res2p

    def table_mutation_complete(self, res2p):
        """

        Notes
        -----

        Parameters
        ----------
        res2p

        Returns
        -------

        """
        pcr_stime = time.time()
        df_mut_info = pd.DataFrame()
        # print(res2p)
        # print(res2p['data_spl'])

        data_pcr = pd.DataFrame(res2p['data_spl'], columns=[
            'read_len',
            'sam_id',
            'source',
            'bead_mut',
            'bead_del',
            'bead_ins',
            'pcr_err_mark',
        ])
        # print(data_pcr)
        data_pcr['num_err_per_read'] = data_pcr['read_len'].apply(lambda x: rannum().binomial(
            n=int(x), p=res2p['pcr_error'], use_seed=False, seed=res2p['ipcr'] + 1
        ))
        data_pcr['pcr_err_mark'] = data_pcr['num_err_per_read'].apply(lambda x: False if x == 0 else True)
        # print(data_pcr['pcr_err_mark'])
        # @@ data_pcr['num_err_per_read']
        # index num_err_per_read
        # 0     0
        # 1     0
        # 2     1
        # ...
        # 41    0
        # 42    0

        df_mut_info['pos_err_per_read'] = data_pcr.apply(lambda x: rannum().uniform(
            low=0, high=x['read_len'], num=x['num_err_per_read'], use_seed=False, seed=res2p['ipcr'] + 1
        ), axis=1)
        df_mut_info['base_roll_per_read'] = data_pcr['num_err_per_read'].apply(lambda x: rannum().uniform(
            low=0, high=3, num=x, use_seed=False
        ))

        # print(data_pcr[['read', 'read_pcr', 'pos_err_per_read']])
        del res2p['data_spl']
        pcr_merge_stime = time.time()
        data_pcr['sam_id'] = data_pcr['sam_id'].apply(lambda x: x + '_' + str(res2p['ipcr'] + 1))
        data_pcr['source'] = 'pcr-' + str(res2p['ipcr'] + 1)

        ### df_mut_info
        # pos_err_per_read base_roll_per_read
        # 0                    []                 []
        # 1                    []                 []
        # 2                    []                 []
        # 3                    []                 []
        # 4                    []                 []
        # ...                 ...                ...
        # 153595               []                 []
        # 153596               []                 []
        # 153597               []                 []
        # 153598               []                 []
        # 153599               []                 []
        df_mut_info['sam_id'] = data_pcr['sam_id'].copy()
        df_mut_info['bead_mut'] = data_pcr['bead_mut'].copy()
        df_mut_info['bead_del'] = data_pcr['bead_del'].copy()
        df_mut_info['bead_ins'] = data_pcr['bead_ins'].copy()
        df_mut_info['pcr_err_mark'] = data_pcr['pcr_err_mark'].copy()

        data_pcr = np.array(data_pcr[[
            'read_len',
            'sam_id',
            'source',
            'bead_mut',
            'bead_del',
            'bead_ins',
            'pcr_err_mark',
        ]])
        # print(data_pcr)
        # print(res2p['data'])

        # print(res2p['mut_info'])
        # print(df_mut_info)
        res2p['data'] = np.concatenate((res2p['data'], data_pcr), axis=0)
        # print(np.concatenate((res2p['mut_info'], np.array(df_mut_info)), axis=0))
        res2p['mut_info'] = np.concatenate((res2p['mut_info'], np.array(df_mut_info)), axis=0)
        self.console.print('======>Time for merging sequences {time:.2f}s'.format(time=time.time() - pcr_merge_stime))
        self.console.print('======>Summary report:')
        self.console.print('=========>PCR time: {time:.2f}s'.format(time=time.time() - pcr_stime))
        self.console.print('=========>Dimension of the data: number of reads: {}'.format(res2p['data'].shape))
        self.console.print('=========>Number of reads at the PCR: {}, '.format(res2p['recorder_pcr_read_num']))
        self.console.print('=========>Number of nucleotides at the PCR: {}, '.format(res2p['recorder_nucleotide_num']))
        self.console.print('=========>Number of errors at the PCR: {}, '.format(res2p['recorder_pcr_err_num']))
        del data_pcr
        return res2p

    def table_tree(self, res2p):
        pcr_stime = time.time()
        data_pcr = pd.DataFrame(
            res2p['data_spl'],
            columns=[
                'sam_id',
                'source',
            ]
        )
        del res2p['data_spl']
        data_pcr['sam_id'] = data_pcr['sam_id'].apply(lambda x: x + '_' + str(res2p['ipcr'] + 1))
        data_pcr['source'] = 'pcr-' + str(res2p['ipcr'] + 1)
        # print(data_pcr)
        data_pcr = np.array(data_pcr)
        res2p['data'] = np.concatenate((res2p['data'], data_pcr), axis=0)
        # print(res2p['data'])
        # print(res2p['data'].shape)
        del data_pcr
        self.console.print('======>Summary report:')
        self.console.print('=========>Number of reads at this PCR: {}'.format(res2p['data'].shape))
        self.console.print('=========>Time for PCR tree construction: {time:.2f}s'.format(time=time.time() - pcr_stime))
        return res2p

    def postableIndexBySameLen(self, seq_len, num_seq):
        """
        If each read has 10 positions, this function returns
        seq_ids
        1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, ...
        seq_pos_ids
        0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, ...

        Parameters
        ----------
        seq_len
        num_seq

        Returns
        -------
        seq_ids
            1d list, the id of the ith read
        seq_pos_ids
            1d list, the base positions (starting from 0 every time for each read) of reads

        """
        nt_ids = [i for i in range(seq_len)]
        seq_pos_ids = nt_ids * num_seq
        seq_ids = np.array([[i] * seq_len for i in range(num_seq)]).ravel().tolist()
        return seq_pos_ids, seq_ids

    def postableLambda(self, x, ids, pos_ids):
        for i in range(x['read_len']):
            pos_ids.append(i)
        ids = ids + [x.name] * x['read_len']
        return ids, pos_ids

    def worktable(self, x, ids):
        if x.name not in ids:
            return x['read']
        else:
            tt = x['read']
            for i in ids[x.name]:
                pcr_err_base = tt[i]
                dna_map = dnasgl().todict(
                    nucleotides=dnasgl().getEleTrimmed(
                        ele_loo=pcr_err_base,
                        universal=True,
                    ),
                    reverse=True,
                )
                tt = tt[:i] + dna_map[2] + tt[i+1:]
            return tt

    def postable_deprecate(self, x, ids, pos_ids):
        for i in range(x['read_len']):
            pos_ids.append(i)
        for i in [x.name] * x['read_len']:
            ids.append(i)
        return ids, pos_ids

    def epos_deprecate(self, spl_nt_ids, ampl_read_len_list):
        """
        ..  @call:
            ------
            arr_err_pos = self.epos(spl_nt_ids=spl_nt_ids, ampl_read_len_list=ampl_read_len_list)

        :param spl_nt_ids: is a list of ids of nucleotides to be modified.
        :param ampl_read_len_list: a list of lengths of reads to be amplified.
        :return:
        """
        pos_err = []
        for i in spl_nt_ids:
            cc = 0
            for id, j in enumerate(ampl_read_len_list):
                cc += j
                if cc >= i:
                    # print(cc, i)
                    pos_err.append([id, (j - 1) - (cc % i)])
                    break
        return pos_err

    def epos2_deprecate(self, ampl_read_len_list):
        ids = []
        pos = []
        for k, e in enumerate(ampl_read_len_list):
            l = [k] * e
            ids = ids + l
            p = [i for i in range(e)]
            pos = pos + p
        return np.array([pos, ids])

    def tactic1(self, arr_2d):
        result = {}
        len_arr = len(arr_2d[0])
        if len_arr == 2:
            for item in arr_2d:
                result[item[0]] = item[1]
        else:
            for item in arr_2d:
                result[item[0]] = item[1:]
        return result

    def todict(self, nucleotides, reverse=False):
        aa_dict = {}
        for k, v in enumerate(nucleotides):
            aa_dict[v] = k
        if reverse:
            aa_dict = {v: k for k, v in aa_dict.items()}
        return aa_dict

    def todict5(self, arr_2d):
        """
        ..  @description:
            -------------
            2d arr to 1d dict, each key in the 1d dict being an list

        :param arr_2d:
        :return:
        """
        result = {}
        for item in arr_2d:
            result[item[0]] = []
        for item in arr_2d:
            if item[0] in result.keys():
                result[item[0]].append(item[1])
                # print(result[item[0]])
        return result