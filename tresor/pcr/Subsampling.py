__version__ = "v1.0"
__copyright__ = "Copyright 2024"
__license__ = "MIT"
__developer__ = "Jianfeng Sun"
__maintainer__ = "Jianfeng Sun"
__email__ = "jianfeng.sunmt@gmail.com"
__lab__ = "Cribbslab"

import numpy as np
from tresor.pcr.Amplify import Amplify as pcr
from tresor.sequencing.Calling import Calling as seq
from tresor.read.error.Library import Library as errlib
from tresor.util.sequence.fastq.Write import write as wfastq
from tresor.util.file.Reader import Reader as pfreader
from tresor.util.random.Number import Number as rannum
from tresor.util.random.Sampling import Sampling as ranspl
from tresor.util.Console import Console
from tresor.util.Kit import tactic6


class Subsampling:

    def __init__(
            self,
            verbose=False
    ):
        self.pcr = pcr
        self.seq = seq
        self.wfastq = wfastq
        self.rannum = rannum
        self.ranspl = ranspl
        self.errlib = errlib()
        self.console = Console()
        self.console.verbose = verbose

    def mutation_table_minimum(
            self,
            pcr_dict,
            replace=False,
    ):
        """

        Parameters
        ----------
        pcr_dict

        Returns
        -------

        """
        self.console.print('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
        self.console.print('======>Substitutions of nucleotides by PCR errors using mutation_table_minimum')
        self.console.print('======>Read PCR amplified reads')
        df_seq_lib = pfreader().generic(pcr_dict['read_lib_fpn'])
        # umi_map = df_seq_lib[0].to_dict()
        import pandas as pd
        umi_map = pd.Series(df_seq_lib[0].values, index=df_seq_lib[1].astype(str)).to_dict()
        # print(umi_map)

        self.console.print('======>Sampling reads to be sequenced')
        num_all_pcr_ampl_reads = pcr_dict['mut_info'].shape[0]
        self.console.print(
            '=========>There are a total number of {} PCR amplified reads'.format(num_all_pcr_ampl_reads))
        # print(num_all_pcr_ampl_reads)
        # print(pcr_dict.keys())
        if pcr_dict['seq_sub_spl_number'] is not None:
            num_reads_for_sequencing = pcr_dict['seq_sub_spl_number']
        else:
            num_reads_for_sequencing = pcr_dict['seq_sub_spl_rate'] * num_all_pcr_ampl_reads
        # spl_ids = rannum().choice(
        #     # low=0,
        #     high=num_all_pcr_ampl_reads,
        #     num=num_reads_for_sequencing,
        #     use_seed=pcr_dict['use_seed'],
        #     seed=pcr_dict['seed'],
        #     replace=replace,
        # )
        spl_ids = rannum().uniform(
            low=0, high=num_all_pcr_ampl_reads, num=num_reads_for_sequencing, use_seed=False, seed=1
        )
        print(spl_ids)
        # print(pcr_dict['data'])
        spl_id_map = tactic6(pcr_dict['data'][:, [1, 2]])
        # print(pcr_dict['data'])
        # print(spl_id_map)
        # print(len(spl_id_map))
        spl_mut_info = pcr_dict['mut_info'][spl_ids]
        # print(pcr_dict['mut_info'].shape)
        # print(spl_mut_info)
        # print(len(spl_mut_info))
        keys = spl_mut_info[:, 2]
        # keys = pcr_dict['data'][spl_ids][:, 1]
        # print(keys)
        # print(np.sort(keys) == np.sort(spl_mut_info[:, 2]))
        # print(len(keys))
        pos_dict = tactic6(pcr_dict['mut_info'][:, [2, 0]])
        base_dict = tactic6(pcr_dict['mut_info'][:, [2, 1]])
        # print(pos_dict)
        # print(base_dict)
        # print(tactic6(base_np))
        res_data = []
        for key in keys:
            mol_id = key.split('_')[0]
            k = key.split('_')[1:]
            # print('kkk', key, k)
            read = umi_map[mol_id]
            # print(read)
            for i in range(len(k)):
                # print('id', i)
                sub_k = mol_id + '_' + '_'.join(k[: i + 1]) if k != [] else mol_id
                # print(sub_k)
                # print(pos_dict[sub_k], base_dict[sub_k])
                if sub_k in pos_dict.keys():
                    read = self.errlib.change(read, pos_list=pos_dict[sub_k], base_list=base_dict[sub_k])
            #     print(read)
            # print(read)
            res_data.append([
                read,  # read
                str(mol_id) + '_' + '_'.join(k) if k != [] else str(mol_id),  # sam id
                spl_id_map[str(mol_id) + '_' + '_'.join(k)] if k != [] else 'init',  # source
            ])
            # print(read)
        # print(np.array(res_data).shape)
        return np.array(res_data)

    def mutation_table_complete(
            self,
            pcr_dict,
            replace=False,
    ):
        # print(pcr_dict)
        print(pcr_dict['mut_info'])
        self.console.print('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
        self.console.print('======>Substitutions of nucleotides by PCR errors using mutation_table_complete')
        self.console.print('======>Read PCR amplified reads')
        df_seq_lib = pfreader().generic(pcr_dict['read_lib_fpn'])
        # umi_map = df_seq_lib[0].to_dict()
        import pandas as pd
        umi_map = pd.Series(df_seq_lib[0].values, index=df_seq_lib[1].astype(str)).to_dict()
        # print(umi_map)
        # num_all_pcr_ampl_reads = pcr_dict['data'].shape[0]
        # self.console.print('=========>There are a total number of {} PCR amplified reads'.format(num_all_pcr_ampl_reads))
        # spl_ids = rannum().uniform(
        #     low=0, high=num_all_pcr_ampl_reads, num=9000, use_seed=False, seed=1
        # )
        self.console.print('======>Sampling reads to be sequenced')
        num_all_pcr_ampl_reads = pcr_dict['data'].shape[0]
        self.console.print('=========>There are a total number of {} PCR amplified reads'.format(num_all_pcr_ampl_reads))
        # print(num_all_pcr_ampl_reads)
        # print(pcr_dict.keys())
        if pcr_dict['seq_sub_spl_number'] is not None:
            num_reads_for_sequencing = pcr_dict['seq_sub_spl_number']
        else:
            num_reads_for_sequencing = pcr_dict['seq_sub_spl_rate'] * num_all_pcr_ampl_reads
        spl_ids = rannum().choice(
            # low=0,
            high=num_all_pcr_ampl_reads,
            num=num_reads_for_sequencing,
            use_seed=pcr_dict['use_seed'],
            seed=pcr_dict['seed'],
            replace=replace,
        )
        # print(spl_ids)
        asd1 = pd.DataFrame(pcr_dict['data'][:, 3:6]).astype(bool)
        print(asd1)
        print(asd1[asd1[1]].shape[0] / asd1.shape[0])
        # print(pcr_dict['data'].shape)
        spl_id_map = tactic6(pcr_dict['data'][:, [1, 2]])
        # print(spl_id_map)
        # print(len(spl_id_map))
        spl_mut_info = pcr_dict['mut_info'][spl_ids]
        # print(pcr_dict['mut_info'].shape)
        asd = pd.DataFrame(spl_mut_info)[[3, 4, 5]].astype(bool)
        print(asd[asd[4]])
        print(asd[asd[4]].shape[0] / asd.shape[0])

        # print(pd.DataFrame(spl_mut_info)[[3, 4, 5]].dtypes)
        # print(len(spl_mut_info))
        keys = spl_mut_info[:, 2]
        # print(len(keys))
        pos_dict = tactic6(pcr_dict['mut_info'][:, [2, 0]])
        base_dict = tactic6(pcr_dict['mut_info'][:, [2, 1]])
        # print(pos_dict)
        # print(base_dict)
        # print(tactic6(base_np))
        res_data = []
        for key in keys:
            mol_id = key.split('_')[0]
            k = key.split('_')[1:]
            # print('kkk', key, k)
            read = umi_map[mol_id]
            # print(read)
            for i in range(len(k)):
                # print('id', i)
                sub_k = mol_id + '_' + '_'.join(k[: i+1]) if k != [] else mol_id
                # print(sub_k)
                # print(pos_dict[sub_k], base_dict[sub_k])
                read = self.errlib.change(read, pos_list=pos_dict[sub_k], base_list=base_dict[sub_k])
            #     print(read)
            # print(read)
            res_data.append([
                read,  # read
                str(mol_id) + '_' + '_'.join(k) if k != [] else str(mol_id),  # sam id
                spl_id_map[str(mol_id) + '_' + '_'.join(k)] if k != [] else 'init',  # source
            ])
            # print(read)
        # print(np.array(res_data).shape)
        return np.array(res_data)

    def sptree(
            self,
            pcr_dict,
            replace=False,
    ):
        """
        short path tree method
        Parameters
        ----------
        pcr_dict

        Returns
        -------

        """
        self.console.print('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
        self.console.print('======>Substitutions of nucleotides by PCR errors using the PCR tree')
        self.console.print('======>Read PCR amplified reads')
        df_seq_lib = pfreader().generic(pcr_dict['read_lib_fpn'])
        # umi_map = df_seq_lib[0].to_dict()
        import pandas as pd
        umi_map = pd.Series(df_seq_lib[0].values, index=df_seq_lib[1].astype(str)).to_dict()
        # @@ single locus
        # {0: 'GGGAAATTTAAACCCTTTAAAGGGAAAAAAGGGCCC', 1: 'GGGTTTAAACCCCCCCCCGGGAAATTTTTTGGGTTT',
        # 2: 'AAATTTGGGCCCGGGAAAGGGCCCAAAAAAGGGAAA', ...,
        # 48: 'CCCAAAAAAGGGCCCAAAGGGCCCTTTGGGGGGCCC', 49: 'TTTTTTAAATTTAAAAAAGGGAAAGGGGGGGGGCCC'}
        # @@ single cell
        # {'0*c*0*g*1*': 'GAAATCATGTAGTTCGGGGGGGCCCTTTTTTTTTTTTAAAAAACCCAAAGGG',
        # '0*c*0*g*4*': 'GAAATCATGTAGTTCGCCCCCCAAATTTAAACCCAAATTTCCCAAAAAAAAA',
        # '1*c*0*g*4*': 'GAAATCATGTAGTTCGAAACCCGGGGGGCCCTTTCCCTTTGGGAAATTTCCC',
        # ...
        # '221*c*1*g*4*': 'CGCGTTAGTAATTCATAAAGGGGGGCCCAAACCCGGGGGGGGGGGGTTTCCC',
        # '222*c*1*g*4*': 'CGCGTTAGTAATTCATAAATTTCCCGGGCCCCCCGGGAAAGGGCCCCCCTTT'}
        # @@ gene
        # {'0*s*0*g*3*': 'GGGGGGCCCTTTTTTTTTTTTAAAAAACCCAAAGGG',
        # '1*s*0*g*3*': 'GGGTTTGGGAAAAAAAAAAAATTTTTTGGGCCCGGG',
        # '2*s*0*g*3*': 'CCCTTTCCCGGGGGGTTTGGGCCCCCCTTTAAAAAA', ...,
        # '421*s*0*g*5*': 'TTTGGGGGGGGGGGGGGGGGGAAAAAAAAAAAATTT',
        # '422*s*0*g*5*': 'AAATTTTTTCCCCCCCCCGGGAAAAAAGGGAAATTT'}

        self.console.print('======>Sampling reads to be sequenced')
        num_all_pcr_ampl_reads = pcr_dict['data'].shape[0]
        self.console.print('=========>There are a total number of {} PCR amplified reads'.format(num_all_pcr_ampl_reads))
        # print(num_all_pcr_ampl_reads)
        # print(pcr_dict.keys())
        if pcr_dict['seq_sub_spl_number'] is not None:
            num_reads_for_sequencing = pcr_dict['seq_sub_spl_number']
        else:
            num_reads_for_sequencing = int(pcr_dict['seq_sub_spl_rate'] * num_all_pcr_ampl_reads)
        spl_ids = rannum().choice(
            # low=0,
            high=num_all_pcr_ampl_reads,
            num=num_reads_for_sequencing,
            use_seed=pcr_dict['use_seed'],
            seed=pcr_dict['seed'],
            replace=replace,
        )
        ### spl_ids
        #@@ single locus
        # [ 5221 12317 12284 ...  3697  7994  7549]
        #@@ single cell
        # [  5587  59571  16821 ... 5900  75148]
        # @@ gene
        # [170263  25800 275989 ... 134570 133424]
        ### len(spl_ids)
        # num_reads_for_sequencing
        spl_data = pcr_dict['data'][spl_ids]
        self.console.print('=========>{} reads selected for sequencing'.format(num_reads_for_sequencing))
        ### spl_data
        #@@ single locus
        # [['23_1_2_5_6_8' 'pcr-8']
        #  ['6_3_4_5_8' 'pcr-8']
        #  ['13_1_2_5_7' 'pcr-7']
        #  ...
        #  ['13_3_4_5_7_9' 'pcr-9']
        #  ['4_3_7_9' 'pcr-9']
        #  ['47_2_5_6_7_9' 'pcr-9']]
        #@@ single cell
        # [['77*c*1*g*4*_2_4_5_6_7' 'pcr-7']
        #  ['32*c*1*g*4*_1_2_4_7' 'pcr-7']
        #  ['196*c*1*g*4*_1' 'pcr-1']
        # ...
        #  ['159*c*1*g*4*_1_2_4_5_9' 'pcr-9']
        #  ['49*c*1*g*4*_2_3_4_5_6_7' 'pcr-7']]
        # @@ gene
        # [['261*s*0*g*4*_1_3_5_7_8_9' 'pcr-9']
        #  ['175*s*0*g*4*_1_4_7' 'pcr-7']
        #  ['193*s*0*g*5*_1_2_6_7_10' 'pcr-10']
        # ...
        #  ['182*s*0*g*5*_1_4_8_9' 'pcr-9']
        #  ['226*s*0*g*4*_3_5_6_9' 'pcr-9']]
        ### len(spl_data)
        # num_reads_for_sequencing
        spl_id_map = tactic6(spl_data)
        ### spl_id_map
        #@@ single locus
        # {'23_1_2_5_6_8': 'pcr-8', '6_3_4_5_8': 'pcr-8', '13_1_2_5_7': 'pcr-7', ..., '34_1_2_4': 'pcr-4', '43_1_8': 'pcr-8'}
        #@@ single cell
        # {'77*c*1*g*4*_2_4_5_6_7': 'pcr-7', '32*c*1*g*4*_1_2_4_7': 'pcr-7', '196*c*1*g*4*_1': 'pcr-1', ..., '159*c*1*g*4*_1_2_4_5_9': 'pcr-9', '49*c*1*g*4*_2_3_4_5_6_7': 'pcr-7'}
        #@@ gene
        # {'261*s*0*g*4*_1_3_5_7_8_9': 'pcr-9', '175*s*0*g*4*_1_4_7': 'pcr-7', '193*s*0*g*5*_1_2_6_7_10': 'pcr-10', ..., '182*s*0*g*5*_1_4_8_9': 'pcr-9', '226*s*0*g*4*_3_5_6_9': 'pcr-9'}
        ### len(spl_id_map)
        # num_reads_for_sequencing
        trees = spl_data[:, 0].ravel().tolist()
        ### trees
        #@@ single locus
        # ['23_1_2_5_6_8', '6_3_4_5_8', '13_1_2_5_7', ..., '34_1_2_4', '43_1_8']
        # @@ single cell
        # ['77*c*1*g*4*_2_4_5_6_7', '32*c*1*g*4*_1_2_4_7', '196*c*1*g*4*_1', ..., '159*c*1*g*4*_1_2_4_5_9', '49*c*1*g*4*_2_3_4_5_6_7']
        # @@ gene
        # ['261*s*0*g*4*_1_3_5_7_8_9', '175*s*0*g*4*_1_4_7', '193*s*0*g*5*_1_2_6_7_10', ..., '182*s*0*g*5*_1_4_8_9', '226*s*0*g*4*_3_5_6_9']
        ### len(trees)
        # num_reads_for_sequencing
        self.console.print('=========>Sampled reads contain {} unrepeated PCR amplified molecules'.format(np.unique(trees).size))

        mol_map_to_all_its_pcr_trees = {}
        for tree in trees:
            mol_map_to_all_its_pcr_trees[tree.split('_')[0]] = []
            # print(tree, tree.split('_')[0])
        for tree in trees:
            mol_map_to_all_its_pcr_trees[tree.split('_')[0]].append(tree.split('_')[1:])
        ### mol_map_to_all_its_pcr_trees
        #@@ single locus
        # {'23': [['1', '2', '5', '6', '8'], ['1', '3', '4', '8', '9', '10'],
        # ['1', '4', '7', '9'], ['2', '3', '4', '5', '7', '9'], ['1', '3', '6'],
        # ['1', '2', '4', '7', '8']], '6': [['3', '4', '5', '8'], ['2', '3', '4', '5', '6', '10'],
        # ['5', '6', '7'], ['1', '2', '4', '6', '7', '9'], ['1', '2', '9']],
        # '13': [['1', '2', '5', '7'],
        # ['1', '4', '7', '8', '10'], ['1', '4', '8', '9'], ['6', '7', '9']],
        # ...,
        # '29': [['2', '5', '9', '10'], ['1', '2', '3', '5', '6', '8', '10'],
        # ['3', '4', '7', '9', '10']],
        # '39': [['6', '7', '9']]}
        #@@ single cell
        # {'77*c*1*g*4*': [['2', '4', '5', '6', '7'], ['1', '3', '4', '7', '8', '10'],
        # ['3', '5', '9', '10'], ['1', '2', '3', '5', '6', '7', '10']], '32*c*1*g*4*': [['1', '2', '4', '7'],
        # ['7', '8', '9', '10']], '196*c*1*g*4*': [['1']], ..., '55*c*1*g*4*': [['2', '3', '8', '9']],
        # '49*c*1*g*4*': [['2', '3', '4', '5', '6', '7']]}
        #@@ gene
        # {'261*s*0*g*4*': [['1', '3', '5', '7', '8', '9']], '175*s*0*g*4*':
        # [['1', '4', '7'], ['1', '4', '5', '8', '9']], '193*s*0*g*5*':
        # [['1', '2', '6', '7', '10']], ..., '171*s*0*g*4*': [['2', '3', '4', '6', '8', '9']],
        # '226*s*0*g*4*': [['3', '5', '6', '9']]}
        ### len(mol_map_to_all_its_pcr_trees)
        # number of unique molecules

        mol_sub_tree_map = {}
        for k, v in mol_map_to_all_its_pcr_trees.items():
            mol_sub_tree_map[k] = []
            for j in v:
                mol_sub_tree_map[k].append('_'.join(j))
        ### mol_sub_tree_map
        #@@ single locus
        # {'23': ['1_2_5_6_8', '1_3_4_8_9_10', '1_4_7_9', '2_3_4_5_7_9',
        # '1_3_6', '1_2_4_7_8'], '6': ['3_4_5_8', '2_3_4_5_6_10', '5_6_7',
        # '1_2_4_6_7_9', '1_2_9'], '13': ['1_2_5_7', '1_4_7_8_10', '1_4_8_9',
        # '6_7_9'], ..., '29': ['2_5_9_10', '1_2_3_5_6_8_10', '3_4_7_9_10'], '39': ['6_7_9']}
        #@@ single cell
        # {'77*c*1*g*4*': ['2_4_5_6_7', '1_3_4_7_8_10', '3_5_9_10', '1_2_3_5_6_7_10'],
        #  '32*c*1*g*4*': ['1_2_4_7', '7_8_9_10'], '196*c*1*g*4*': ['1'], ...,
        #  '55*c*1*g*4*': ['2_3_8_9'], '49*c*1*g*4*': ['2_3_4_5_6_7']}
        #@@ gene
        # {'261*s*0*g*4*': ['1_3_5_7_8_9'], '175*s*0*g*4*': ['1_4_7', '1_4_5_8_9'],
        # '193*s*0*g*5*': ['1_2_6_7_10'], ..., '171*s*0*g*4*': ['2_3_4_6_8_9'],
        # '226*s*0*g*4*': ['3_5_6_9']}

        res_data = []
        for k, sub_trees in mol_sub_tree_map.items():
            ### k, sub_trees
            #@@ single locus
            # 23 ['1_2_5_6_8', '1_3_4_8_9_10', '1_4_7_9', '2_3_4_5_7_9', '1_3_6', '1_2_4_7_8']
            # 6 ['3_4_5_8', '2_3_4_5_6_10', '5_6_7', '1_2_4_6_7_9', '1_2_9']
            # 13 ['1_2_5_7', '1_4_7_8_10', '1_4_8_9', '6_7_9']
            # ...
            # 29 ['2_5_9_10', '1_2_3_5_6_8_10', '3_4_7_9_10']
            # 39 ['6_7_9']
            #@@ single cell
            # 77*c*1*g*4* ['2_4_5_6_7', '1_3_4_7_8_10', '3_5_9_10', '1_2_3_5_6_7_10']
            # 32*c*1*g*4* ['1_2_4_7', '7_8_9_10']
            # 196*c*1*g*4* ['1']
            # ...
            # 55*c*1*g*4* ['2_3_8_9']
            # 49*c*1*g*4* ['2_3_4_5_6_7']
            #@@ gene
            # 261*s*0*g*4* ['1_3_5_7_8_9']
            # 175*s*0*g*4* ['1_4_7', '1_4_5_8_9']
            # 193*s*0*g*5* ['1_2_6_7_10']
            # ...
            # 171*s*0*g*4* ['2_3_4_6_8_9']
            # 226*s*0*g*4* ['3_5_6_9']

            read = umi_map[k]
            # read = umi_map[int(k)]
            read_cache = {}
            for sub_tree in sub_trees:
                ### k, sub_tree
                #@@ single locus
                # 23 1_2_5_6_8
                # 23 1_3_4_8_9_10
                # 23 1_4_7_9
                # 23 2_3_4_5_7_9
                # 23 1_3_6
                # 23 1_2_4_7_8
                # 6 3_4_5_8
                # 6 2_3_4_5_6_10
                # 6 5_6_7
                # 6 1_2_4_6_7_9
                # 6 1_2_9
                # 13 1_2_5_7
                # 13 1_4_7_8_10
                # 13 1_4_8_9
                # 13 6_7_9
                # ...
                # 29 2_5_9_10
                # 29 1_2_3_5_6_8_10
                # 29 3_4_7_9_10
                # 39 6_7_9

                #@@ single cell
                # 77*c*1*g*4* ['2_4_5_6_7', '1_3_4_7_8_10', '3_5_9_10', '1_2_3_5_6_7_10']
                # 77*c*1*g*4* 2_4_5_6_7
                # 77*c*1*g*4* 1_3_4_7_8_10
                # ...
                # 49*c*1*g*4* ['2_3_4_5_6_7']
                # 49*c*1*g*4* 2_3_4_5_6_7

                #@@ gene
                # 261*s*0*g*4* 1_3_5_7_8_9
                # 175*s*0*g*4* 1_4_7
                # 175*s*0*g*4* 1_4_5_8_9
                # ...
                # 171*s*0*g*4* 2_3_4_6_8_9
                # 226*s*0*g*4* 3_5_6_9

                if len(sub_tree) > 0:
                    read_id = k + '_' + sub_tree
                else: # sub_tree can be ''. len('') == 0
                    read_id = k
                ### read_id
                #@@ single locus
                # 23_1_2_5_6_8
                # 23_1_3_4_8_9_10
                # 23_1_4_7_9
                # 23_2_3_4_5_7_9
                # 23_1_3_6
                # 23_1_2_4_7_8
                # 6_3_4_5_8
                # 6_2_3_4_5_6_10
                # 6_5_6_7
                # 6_1_2_4_6_7_9
                # 6_1_2_9
                # 13_1_2_5_7
                # 13_1_4_7_8_10
                # 13_1_4_8_9
                # 13_6_7_9
                # ...
                # 29_2_5_9_10
                # 29_1_2_3_5_6_8_10
                # 29_3_4_7_9_10
                # 39_6_7_9

                #@@ single cell
                # 77*c*1*g*4*_2_4_5_6_7
                # 77*c*1*g*4*_1_3_4_7_8_10
                # 77*c*1*g*4*_3_5_9_10
                # ...
                # 55*c*1*g*4*_2_3_8_9
                # 49*c*1*g*4*_2_3_4_5_6_7

                #@@ gene
                # 261*s*0*g*4*_1_3_5_7_8_9
                # 175*s*0*g*4*_1_4_7
                # 175*s*0*g*4*_1_4_5_8_9
                # ...
                # 171*s*0*g*4*_2_3_4_6_8_9
                # 226*s*0*g*4*_3_5_6_9

                k_ = k
                sub_tree_pcr_arr = sub_tree.split('_')
                # print(sub_tree_pcr_arr)
                for pcr in sub_tree_pcr_arr:
                    if pcr != '':
                        k_ = k_ + '_' + pcr
                    ### k, '<======>', k_, 'read id:', read_id
                    ### k, k_
                    #@@ single locus
                    # 23 <======> 23_1 read id: 23_1_2_5_6_8
                    # 23 <======> 23_1_2 read id: 23_1_2_5_6_8
                    # 23 <======> 23_1_2_5 read id: 23_1_2_5_6_8
                    # 23 <======> 23_1_2_5_6 read id: 23_1_2_5_6_8
                    # 23 <======> 23_1_2_5_6_8 read id: 23_1_2_5_6_8
                    # 23 <======> 23_1 read id: 23_1_3_4_8_9_10
                    # 23 <======> 23_1_3 read id: 23_1_3_4_8_9_10
                    # 23 <======> 23_1_3_4 read id: 23_1_3_4_8_9_10
                    # 23 <======> 23_1_3_4_8 read id: 23_1_3_4_8_9_10
                    # 23 <======> 23_1_3_4_8_9 read id: 23_1_3_4_8_9_10
                    # 23 <======> 23_1_3_4_8_9_10 read id: 23_1_3_4_8_9_10
                    # 23 <======> 23_1 read id: 23_1_4_7_9
                    # 23 <======> 23_1_4 read id: 23_1_4_7_9
                    # 23 <======> 23_1_4_7 read id: 23_1_4_7_9
                    # 23 <======> 23_1_4_7_9 read id: 23_1_4_7_9
                    # 23 <======> 23_2 read id: 23_2_3_4_5_7_9
                    # ...

                    #@@ single cell
                    # 55*c*1*g*4* <======> 55*c*1*g*4*_2 read id: 55*c*1*g*4*_2_3_8_9
                    # 55*c*1*g*4* <======> 55*c*1*g*4*_2_3 read id: 55*c*1*g*4*_2_3_8_9
                    # 55*c*1*g*4* <======> 55*c*1*g*4*_2_3_8 read id: 55*c*1*g*4*_2_3_8_9
                    # 55*c*1*g*4* <======> 55*c*1*g*4*_2_3_8_9 read id: 55*c*1*g*4*_2_3_8_9
                    # 49*c*1*g*4* <======> 49*c*1*g*4*_2 read id: 49*c*1*g*4*_2_3_4_5_6_7
                    # 49*c*1*g*4* <======> 49*c*1*g*4*_2_3 read id: 49*c*1*g*4*_2_3_4_5_6_7
                    # 49*c*1*g*4* <======> 49*c*1*g*4*_2_3_4 read id: 49*c*1*g*4*_2_3_4_5_6_7
                    # 49*c*1*g*4* <======> 49*c*1*g*4*_2_3_4_5 read id: 49*c*1*g*4*_2_3_4_5_6_7
                    # 49*c*1*g*4* <======> 49*c*1*g*4*_2_3_4_5_6 read id: 49*c*1*g*4*_2_3_4_5_6_7
                    # 49*c*1*g*4* <======> 49*c*1*g*4*_2_3_4_5_6_7 read id: 49*c*1*g*4*_2_3_4_5_6_7
                    # ...

                    #@@ gene
                    # 261*s*0*g*4* <======> 261*s*0*g*4*_1 read id: 261*s*0*g*4*_1_3_5_7_8_9
                    # 261*s*0*g*4* <======> 261*s*0*g*4*_1_3 read id: 261*s*0*g*4*_1_3_5_7_8_9
                    # 261*s*0*g*4* <======> 261*s*0*g*4*_1_3_5 read id: 261*s*0*g*4*_1_3_5_7_8_9
                    # ...

                    if k_ in read_cache.keys():
                        read = read_cache[k_]
                    else:
                        read = self.errlib.mutated(
                            read=read,
                            pcr_error=pcr_dict['pcr_error'],
                        )
                    if pcr_dict['pcr_deletion']:
                        read = self.errlib.deletion(read=read, del_rate=pcr_dict['pcr_del_rate'])
                    if pcr_dict['pcr_insertion']:
                        read = self.errlib.insertion(read=read, ins_rate=pcr_dict['pcr_ins_rate'])
                    read_cache[k_] = read
                # print(read_cache)
                res_data.append([
                    read,
                    read_id,
                    spl_id_map[read_id]
                ])
        # print(res_data)
        # print(len(res_data))
        return np.array(res_data)

    def bftree(
            self,
            pcr_dict,
            replace=False,
    ):
        """
        boolean flag tree method

        Notes
        -----
            bool_flags_ = [[] for _ in range(len(uniq_mol_map_new))]
            realvalued_flags_ = [[] for _ in range(len(uniq_mol_map_new))]
            uniq_mol_map_new_ = [[] for _ in range(len(uniq_mol_map_new))]
            for ii, u in enumerate(uniq_mol_map_new):
                for jj, v in enumerate(u):
                    if v != '0':
                        bool_flags_[ii].append(bool_flag_table[ii][jj])
                        realvalued_flags_[ii].append(realvalued_flag_table[ii][jj])
                        uniq_mol_map_new_[ii].append(uniq_mol_map_new[ii][jj])
            print(bool_flags_)
            print(realvalued_flags_)
            print(uniq_mol_map_new_)

        Other Parameters
        ----------------
        trees (i.e., PCR tree)
            [['2', '5', '6', '12', '13', '15'],
             ['3', '4', '7', '9', '12'],
             ['2', '6', '9', '10', '11', '12', '13', '15'],
             ['2', '4', '6', '8', '11'],
             ['3', '4', '7', '9', '13', '14']]
        trees_np
            [['2' '5' '6' '12' '13' '15' '0' '0']
             ['3' '4' '7' '9' '12' '0' '0' '0']
             ['2' '6' '9' '10' '11' '12' '13' '15']
             ['2' '4' '6' '8' '11' '0' '0' '0']
             ['3' '4' '7' '9' '13' '14' '0' '0']]
        bool_flag_table
            [[ True False False False False False False False]
             [ True  True  True  True False False False False]
             [ True False False False False False False False]
             [ True  True False False False False False False]
             [ True  True  True  True False False False False]]
        realvalued_flag_table
            [['2' '-1' '-1' '-1' '-1' '-1' '-1' '-1']
             ['3' '4' '7' '9' '-1' '-1' '-1' '-1']
             ['2' '-1' '-1' '-1' '-1' '-1' '-1' '-1']
             ['2' '4' '-1' '-1' '-1' '-1' '-1' '-1']
             ['3' '4' '7' '9' '-1' '-1' '-1' '-1']]

        Returns
        -------

        """
        self.console.print('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
        self.console.print('======>Substitutions of nucleotides by PCR errors using the PCR tree')
        self.console.print('======>Read PCR amplified reads')
        df_seq_lib = pfreader().generic(pcr_dict['read_lib_fpn'])
        # umi_map = df_seq_lib[0].to_dict()
        import pandas as pd
        umi_map = pd.Series(df_seq_lib[0].values, index=df_seq_lib[1].astype(str)).to_dict()
        # print(umi_map)

        self.console.print('======>Sampling reads to be sequenced')
        num_all_pcr_ampl_reads = pcr_dict['data'].shape[0]
        self.console.print('=========>There are a total number of {} PCR amplified reads'.format(num_all_pcr_ampl_reads))
        # print(num_all_pcr_ampl_reads)
        # print(pcr_dict.keys())
        if pcr_dict['seq_sub_spl_number'] is not None:
            num_reads_for_sequencing = pcr_dict['seq_sub_spl_number']
        else:
            num_reads_for_sequencing = pcr_dict['seq_sub_spl_rate'] * num_all_pcr_ampl_reads
        spl_ids = rannum().choice(
            # low=0,
            high=num_all_pcr_ampl_reads,
            num=num_reads_for_sequencing,
            use_seed=pcr_dict['use_seed'],
            seed=pcr_dict['seed'],
            replace=replace,
        )
        ### spl_ids
        # [ 5221 12317 12284 ...  3697  7994  7549]
        ### len(spl_ids)
        # num_reads_for_sequencing

        spl_data = pcr_dict['data'][spl_ids]
        self.console.print('=========>{} reads selected for sequencing'.format(num_reads_for_sequencing))
        ### spl_data
        # [['23_1_2_5_6_8' 'pcr-8']
        #  ['6_3_4_5_8' 'pcr-8']
        #  ['13_1_2_5_7' 'pcr-7']
        #  ...
        #  ['13_3_4_5_7_9' 'pcr-9']
        #  ['4_3_7_9' 'pcr-9']
        #  ['47_2_5_6_7_9' 'pcr-9']]
        ### len(spl_data)
        # num_reads_for_sequencing
        spl_id_map = tactic6(spl_data)
        ### spl_id_map
        # {'23_1_2_5_6_8': 'pcr-8', '6_3_4_5_8': 'pcr-8', '13_1_2_5_7': 'pcr-7', ..., '34_1_2_4': 'pcr-4', '43_1_8': 'pcr-8'}
        ### len(spl_id_map)
        # num_reads_for_sequencing
        trees = spl_data[:, 0].ravel().tolist()
        # print(trees)
        ### trees
        # ['23_1_2_5_6_8', '6_3_4_5_8', '13_1_2_5_7', ... '34_1_2_4', '43_1_8']
        ### len(trees)
        # num_reads_for_sequencing
        self.console.print('=========>Sampled reads contain {} unrepeated PCR amplified molecules'.format(np.unique(trees).size))

        res_data = []
        mol_map_to_all_its_pcr_trees = {}
        for tree in trees:
            mol_map_to_all_its_pcr_trees[tree.split('_')[0]] = []
            # print(tree, tree.split('_')[0])
        for tree in trees:
            mol_map_to_all_its_pcr_trees[tree.split('_')[0]].append(tree.split('_')[1:])
        ### mol_map_to_all_its_pcr_trees
        # {'23': [['1', '2', '5', '6', '8'], ['1', '3', '4', '8', '9', '10'],
        # ['1', '4', '7', '9'], ['2', '3', '4', '5', '7', '9'], ['1', '3', '6'],
        # ['1', '2', '4', '7', '8']], '6': [['3', '4', '5', '8'],
        # ['2', '3', '4', '5', '6', '10'], ['5', '6', '7'], ['1', '2', '4', '6', '7', '9'],
        # ['1', '2', '9']], '13': [['1', '2', '5', '7'], ['1', '4', '7', '8', '10'],
        # ['1', '4', '8', '9'], ['6', '7', '9']],
        # ...
        # '48': [['1'], ['3', '4', '5', '8', '9', '10'], ['1', '9', '10'],
        # ['3', '4', '5', '7', '9', '10']], '11': [['1', '3', '5', '6', '7'],
        # ['2', '3', '4', '7', '8', '10'], ['1', '2', '3', '6', '9']],
        # '29': [['2', '5', '9', '10'], ['1', '2', '3', '5', '6', '8', '10'],
        # ['3', '4', '7', '9', '10']], '39': [['6', '7', '9']]}
        self.console.print('=========>Sampled reads cover {} initial molecules'.format(len(mol_map_to_all_its_pcr_trees)))
        ### len(mol_map_to_all_its_pcr_trees)
        # 50

        for k, trees in mol_map_to_all_its_pcr_trees.items():
            # k, trees
            # 23 [['1', '2', '5', '6', '8'], ['1', '3', '4', '8', '9', '10'], ['1', '4', '7', '9'],
            # ['2', '3', '4', '5', '7', '9'], ['1', '3', '6'], ['1', '2', '4', '7', '8']]
            # 6 [['3', '4', '5', '8'], ['2', '3', '4', '5', '6', '10'], ['5', '6', '7'],
            # ['1', '2', '4', '6', '7', '9'], ['1', '2', '9']]
            # 13 [['1', '2', '5', '7'], ['1', '4', '7', '8', '10'], ['1', '4', '8', '9'], ['6', '7', '9']
            # ...
            # 29 [['2', '5', '9', '10'], ['1', '2', '3', '5', '6', '8', '10'], ['3', '4', '7', '9', '10']]
            # 39 [['6', '7', '9']]
            read = umi_map[str(k)]
            # read
            # CCCCCCGGGGGGAAAAAATTTCCCTTTGGGAAACCC
            read_cache_table = [[] for _ in range(len(trees))]
            bool_flag_table = [[True] for _ in range(len(trees))]
            # realvalued_flag_table = [[] for _ in range(len(trees))]
            trees_ori = [[] for _ in range(len(trees))]
            max_len_pcr_tree = max([len(tree) for tree in trees])
            # print(max_len_pcr_tree)
            for _, tree in enumerate(trees):
                if len(tree) < max_len_pcr_tree:
                    tree += ['0' for _ in range(max_len_pcr_tree - len(tree))]
                # k, tree
                # 23 ['1', '2', '5', '6', '8', '0']
                # 23 ['1', '3', '4', '8', '9', '10']
                # 23 ['1', '4', '7', '9', '0', '0']
                # 23 ['2', '3', '4', '5', '7', '9']
                # 23 ['1', '3', '6', '0', '0', '0']
                # 23 ['1', '2', '4', '7', '8', '0']
                # 6 ['3', '4', '5', '8', '0', '0']
                # 6 ['2', '3', '4', '5', '6', '10']
                # 6 ['5', '6', '7', '0', '0', '0']
                # 6 ['1', '2', '4', '6', '7', '9']
                # 6 ['1', '2', '9', '0', '0', '0']
                # 13 ['1', '2', '5', '7', '0']
                # 13 ['1', '4', '7', '8', '10']
                # 13 ['1', '4', '8', '9', '0']
                # 13 ['6', '7', '9', '0', '0']
                # ...
                # 29 ['2', '5', '9', '10', '0', '0', '0']
                # 29 ['1', '2', '3', '5', '6', '8', '10']
                # 29 ['3', '4', '7', '9', '10', '0', '0']
                # 39 ['6', '7', '9']
            trees_np = np.array(trees)
            ### k, trees_np
            # 23 [['1' '2' '5' '6' '8' '0']
            #  ['1' '3' '4' '8' '9' '10']
            #  ['1' '4' '7' '9' '0' '0']
            #  ['2' '3' '4' '5' '7' '9']
            #  ['1' '3' '6' '0' '0' '0']
            #  ['1' '2' '4' '7' '8' '0']]
            # 6 [['3' '4' '5' '8' '0' '0']
            #  ['2' '3' '4' '5' '6' '10']
            #  ['5' '6' '7' '0' '0' '0']
            #  ['1' '2' '4' '6' '7' '9']
            #  ['1' '2' '9' '0' '0' '0']]
            # 13 [['1' '2' '5' '7' '0']
            #  ['1' '4' '7' '8' '10']
            #  ['1' '4' '8' '9' '0']
            #  ['6' '7' '9' '0' '0']]
            # ...
            # 29 [['2' '5' '9' '10' '0' '0' '0']
            #  ['1' '2' '3' '5' '6' '8' '10']
            #  ['3' '4' '7' '9' '10' '0' '0']]
            # 39 [['6' '7' '9']]

            ### +++++++++++++++ block: Construct bool and real-valued flag tables +++++++++++++++
            #         trees (i.e., PCR tree)
            #             [['2', '5', '6', '12', '13', '15'],
            #              ['3', '4', '7', '9', '12'],
            #              ['2', '6', '9', '10', '11', '12', '13', '15'],
            #              ['2', '4', '6', '8', '11'],
            #              ['3', '4', '7', '9', '13', '14']]
            #         trees_np
            #             [['2' '5' '6' '12' '13' '15' '0' '0']
            #              ['3' '4' '7' '9' '12' '0' '0' '0']
            #              ['2' '6' '9' '10' '11' '12' '13' '15']
            #              ['2' '4' '6' '8' '11' '0' '0' '0']
            #              ['3' '4' '7' '9' '13' '14' '0' '0']]
            #         bool_flag_table
            #             [[ True False False False False False False False]
            #              [ True  True  True  True False False False False]
            #              [ True False False False False False False False]
            #              [ True  True False False False False False False]
            #              [ True  True  True  True False False False False]]
            #         realvalued_flag_table
            #             [['2' '-1' '-1' '-1' '-1' '-1' '-1' '-1']
            #              ['3' '4' '7' '9' '-1' '-1' '-1' '-1']
            #              ['2' '-1' '-1' '-1' '-1' '-1' '-1' '-1']
            #              ['2' '4' '-1' '-1' '-1' '-1' '-1' '-1']
            #              ['3' '4' '7' '9' '-1' '-1' '-1' '-1']]
            ### +++++++++++++++ block: construct bool and real-valued flag tables +++++++++++++++
            for id_horiz_in_a_tree in range(trees_np.shape[1]):
                repeat_in_a_col = self.findListDuplicates(trees_np[:, id_horiz_in_a_tree])
                # print('repeated in a col', repeat_in_a_col)
                # print(trees_np[:, id_horiz_in_a_tree])
                for id_vert_across_trees, ele_in_a_col in enumerate(trees_np[:, id_horiz_in_a_tree]):
                    # print('ele_in_a_col {}'.format(ele_in_a_col))
                    if ele_in_a_col in repeat_in_a_col:
                        ids_repeat_in_a_col = [i for i, value in enumerate(trees_np[:, id_horiz_in_a_tree]) if value == ele_in_a_col]
                        # print('asds {}'.format(ids_repeat_in_a_col))
                        if bool_flag_table[id_vert_across_trees][id_horiz_in_a_tree] is False:
                            # print('repeated ele: {}'.format(ele_in_a_col))
                            bool_flag_table[id_vert_across_trees].append(False)
                            # realvalued_flag_table[id_vert_across_trees].append(-1)
                        else:
                            # @@ any positions that is overlay
                            # print('non-repeated ele: {}'.format(ele_in_a_col))
                            inspector_flags = [1 if bool_flag_table[i][id_horiz_in_a_tree] is True else 0 for i in ids_repeat_in_a_col]
                            # print('inspector_flags {}'.format(inspector_flags))
                            if sum(inspector_flags) > 1:
                                bool_flag_table[id_vert_across_trees].append(True)
                                # realvalued_flag_table[id_vert_across_trees].append(ele_in_a_col)
                            else:
                                bool_flag_table[id_vert_across_trees].append(False)
                                # realvalued_flag_table[id_vert_across_trees].append(-1)
                    else:
                        bool_flag_table[id_vert_across_trees].append(False)
                        # realvalued_flag_table[id_vert_across_trees].append(-1)

            bool_flag_table = np.array(bool_flag_table)[:, 1:]
            ### trees_np
            ### bool_flag_table
            # 23
            # [['1' '2' '5' '6' '8' '0']
            #  ['1' '3' '4' '8' '9' '10']
            #  ['1' '4' '7' '9' '0' '0']
            #  ['2' '3' '4' '5' '7' '9']
            #  ['1' '3' '6' '0' '0' '0']
            #  ['1' '2' '4' '7' '8' '0']]
            # [[ True  True False False False False]
            #  [ True  True  True False False False]
            #  [ True False False False False False]
            #  [False False False False False False]
            #  [ True  True False False False False]
            #  [ True  True  True False False False]]
            # 6 [['3' '4' '5' '8' '0' '0']
            #  ['2' '3' '4' '5' '6' '10']
            #  ['5' '6' '7' '0' '0' '0']
            #  ['1' '2' '4' '6' '7' '9']
            #  ['1' '2' '9' '0' '0' '0']]
            # [[False False False False False False]
            #  [False False False False False False]
            #  [False False False False False False]
            #  [ True  True False False False False]
            #  [ True  True False False False False]]
            # 13 [['1' '2' '5' '7' '0']
            #  ['1' '4' '7' '8' '10']
            #  ['1' '4' '8' '9' '0']
            #  ['6' '7' '9' '0' '0']]
            # [[ True False False False False]
            #  [ True  True False False False]
            #  [ True  True False False False]
            #  [False False False False False]]
            # ...
            # 29 [['2' '5' '9' '10' '0' '0' '0']
            #  ['1' '2' '3' '5' '6' '8' '10']
            #  ['3' '4' '7' '9' '10' '0' '0']]
            # [[False False False False False False False]
            #  [False False False False False False False]
            #  [False False False False False False False]]
            # 39 [['6' '7' '9']]
            # [[False False False]]

            for jj in range(trees_np.shape[1]):
                read_for_repeat_tmp_per_col = {}
                for ii, val_in_a_col in enumerate(trees_np[:, jj]):
                    if val_in_a_col != '0':
                        if jj == 0:
                            # print(val_in_a_col, bool_flag_table[ii][jj])
                            if bool_flag_table[ii][jj] == True:
                                if val_in_a_col not in [*read_for_repeat_tmp_per_col.keys()]:
                                    r1 = self.errlib.mutated(
                                        read=read,
                                        pcr_error=pcr_dict['pcr_error'],
                                    )
                                    if pcr_dict['pcr_deletion']:
                                        r1 = self.errlib.deletion(read=r1, del_rate=pcr_dict['pcr_del_rate'])
                                    if pcr_dict['pcr_insertion']:
                                        r1 = self.errlib.insertion(read=r1, ins_rate=pcr_dict['pcr_ins_rate'])
                                    read_for_repeat_tmp_per_col[val_in_a_col] = r1
                                    read_cache_table[ii].append(r1)
                                else:
                                    read_cache_table[ii].append(read_for_repeat_tmp_per_col[val_in_a_col])
                            else:
                                r1 = self.errlib.mutated(
                                    read=read,
                                    pcr_error=pcr_dict['pcr_error'],
                                )
                                if pcr_dict['pcr_deletion']:
                                    r1 = self.errlib.deletion(read=r1, del_rate=pcr_dict['pcr_del_rate'])
                                if pcr_dict['pcr_insertion']:
                                    r1 = self.errlib.insertion(read=r1, ins_rate=pcr_dict['pcr_ins_rate'])
                                read_cache_table[ii].append(r1)
                        if jj > 0:
                            if bool_flag_table[ii][jj] == True:
                                if val_in_a_col + '_' + '_'.join(list(trees_np[ii][:jj])) not in [*read_for_repeat_tmp_per_col.keys()]:
                                    r1 = self.errlib.mutated(
                                        read=read_cache_table[ii][jj - 1],
                                        pcr_error=pcr_dict['pcr_error'],
                                    )
                                    if pcr_dict['pcr_deletion']:
                                        r1 = self.errlib.deletion(read=r1, del_rate=pcr_dict['pcr_del_rate'])
                                    if pcr_dict['pcr_insertion']:
                                        r1 = self.errlib.insertion(read=r1, ins_rate=pcr_dict['pcr_ins_rate'])
                                    read_for_repeat_tmp_per_col[val_in_a_col + '_' + '_'.join(list(trees_np[ii][:jj]))] = r1
                                    read_cache_table[ii].append(r1)
                                else:
                                    read_cache_table[ii].append(read_for_repeat_tmp_per_col[val_in_a_col + '_' + '_'.join(list(trees_np[ii][:jj]))])
                                # print('id', [i for i, value in enumerate(trees_np[:, jj]) if value == val_in_a_col])
                            else:
                                r1 = self.errlib.mutated(
                                    read=read_cache_table[ii][jj - 1],
                                    pcr_error=pcr_dict['pcr_error'],
                                )
                                if pcr_dict['pcr_deletion']:
                                    r1 = self.errlib.deletion(read=r1, del_rate=pcr_dict['pcr_del_rate'])
                                if pcr_dict['pcr_insertion']:
                                    r1 = self.errlib.insertion(read=r1, ins_rate=pcr_dict['pcr_ins_rate'])
                                read_cache_table[ii].append(r1)
            ### k, read_cache_table
            # 23 [['CCCCCCGGGGGGAAAAAATTTCCCTTTGGGAAACCC', 'CCCCCCGGGGGGAAAAAATTTCCCTTTGGGAAACCC', 'CCCCCCGGGGGGAAAAAATTTCCCTTTGGGAAACCC', 'CCCCCCGGGGGGAAAAAATTTCCCTTTGGGAAACCC', 'CCCCCCGGGGGGAAAAAATTTCCCTTTGGGAAACCC'], ['CCCCCCGGGGGGAAAAAATTTCCCTTTGGGAAACCC', 'CCCCCCGGGGGGAAAAAATTTCCCTTTGGGAAACCC', 'CCCCCCGGGGGGAAAAAATTTCCCTTTGGGAAACCC', 'CCCCCCGGGGGGAAAAAATTTCCCTTTGGGAAACCC', 'CCCCCCGGGGGGAAAAAATTTCCCTTTGGGAAACCC', 'CCCCCCGGGGGGAAAAAATTTCCCTTTGGGAAACCC'], ['CCCCCCGGGGGGAAAAAATTTCCCTTTGGGAAACCC', 'CCCCCCGGGGGGAAAAAATTTCCCTTTGGGAAACCC', 'CCCCCCGGGGGGAAAAAATTTCCCTTTGGGAAACCC', 'CCCCCCGGGGGGAAAAAATTTCCCTTTGGGAAACCC'], ['CCCCCCGGGGGGAAAAAATTTCCCTTTGGGAAACCC', 'CCCCCCGGGGGGAAAAAATTTCCCTTTGGGAAACCC', 'CCCCCCGGGGGGAAAAAATTTCCCTTTGGGAAACCC', 'CCCCCCGGGGGGAAAAAATTTCCCTTTGGGAAACCC', 'CCCCCCGGGGGGAAAAAATTTCCCTTTGGGAAACCC', 'CCCCCCGGGGGGAAAAAATTTCCCTTTGGGAAACCC'], ['CCCCCCGGGGGGAAAAAATTTCCCTTTGGGAAACCC', 'CCCCCCGGGGGGAAAAAATTTCCCTTTGGGAAACCC', 'CCCCCCGGGGGGAAAAAATTTCCCTTTGGGAAACCC'], ['CCCCCCGGGGGGAAAAAATTTCCCTTTGGGAAACCC', 'CCCCCCGGGGGGAAAAAATTTCCCTTTGGGAAACCC', 'CCCCCCGGGGGGAAAAAATTTCCCTTTGGGAAACCC', 'CCCCCCGGGGGGAAAAAATTTCCCTTTGGGAAACCC', 'CCCCCCGGGGGGAAAAAATTTCCCTTTGGGAAACCC']]
            # 6 [['AAACCCAAAAAACCCTTTCCCAAACCCGGGGGGTTT', 'AAACCCAAAAAACCCTTTCCCAAACCCGGGGGGTTT', 'AAACCCAAAAAACCCTTTCCCAAACCCGGGGGGTTT', 'AAACCCAAAAAACCCTTTCCCAAACCCGGGGGGTTT'], ['AAACCCAAAAAACCCTTTCCCAAACCCGGGGGGTTT', 'AAACCCAAAAAACCCTTTCCCAAACCCGGGGGGTTT', 'AAACCCAAAAAACCCTTTCCCAAACCCGGGGGGTTT', 'AAACCCAAAAAACCCTTTCCCAAACCCGGGGGGTTT', 'AAACCCAAAAAACCCTTTCCCAAACCCGGGGGGTTT', 'AAACCCAAAAAACCCTTTCCCAAACCCGGGGGGTTT'], ['AAACCCAAAAAACCCTTTCCCAAACCCGGGGGGTTT', 'AAACCCAAAAAACCCTTTCCCAAACCCGGGGGGTTT', 'AAACCCAAAAAACCCTTTCCCAAACCCGGGGGGTTT'], ['AAACCCAAAAAACCCTTTCCCAAACCCGGGGGGTTT', 'AAACCCAAAAAACCCTTTCCCAAACCCGGGGGGTTT', 'AAACCCAAAAAACCCTTTCCCAAACCCGGGGGGTTT', 'AAACCCAAAAAACCCTTTCCCAAACCCGGGGGGTTT', 'AAACCCAAAAAACCCTTTCCCAAACCCGGGGGGTTT', 'AAACCCAAAAAACCCTTTCCCAAACCCGGGGGGTTT'], ['AAACCCAAAAAACCCTTTCCCAAACCCGGGGGGTTT', 'AAACCCAAAAAACCCTTTCCCAAACCCGGGGGGTTT', 'AAACCCAAAAAACCCTTTCCCAAACCCGGGGGGTTT']]
            # 13 [['CCCAAATTTTTTAAAGGGTTTGGGCCCGGGGGGGGG', 'CCCAAATTTTTTAAAGGGTTTGGGCCCGGGGGGGGG', 'CCCAAATTTTTTAAAGGGTTTGGGCCCGGGGGGGGG', 'CCCAAATTTTTTAAAGGGTTTGGGCCCGGGGGGGGG'], ['CCCAAATTTTTTAAAGGGTTTGGGCCCGGGGGGGGG', 'CCCAAATTTTTTAAAGGGTTTGGGCCCGGGGGGGGG', 'CCCAAATTTTTTAAAGGGTTTGGGCCCGGGGGGGGG', 'CCCAAATTTTTTAAAGGGTTTGGGCCCGGGGGGGGG', 'CCCAAATTTTTTAAAGGGTTTGGGCCCGGGGGGGGG'], ['CCCAAATTTTTTAAAGGGTTTGGGCCCGGGGGGGGG', 'CCCAAATTTTTTAAAGGGTTTGGGCCCGGGGGGGGG', 'CCCAAATTTTTTAAAGGGTTTGGGCCCGGGGGGGGG', 'CCCAAATTTTTTAAAGGGTTTGGGCCCGGGGGGGGG'], ['CCCAAATTTTTTAAAGGGTTTGGGCCCGGGGGGGGG', 'CCCAAATTTTTTAAAGGGTTTGGGCCCGGGGGGGGG', 'CCCAAATTTTTTAAAGGGTTTGGGCCCGGGGGGGGG']]
            # ...
            # 29 [['TTTCCCCCCAAACCCTTTGGGGGGGGGAAAGGGTTT', 'TTTCCCCCCAAACCCTTTGGGGGGGGGAAAGGGTTT', 'TTTCCCCCCAAACCCTTTGGGGGGGGGAAAGGGTTT', 'TTTCCCCCCAAACCCTTTGGGGGGGGGAAAGGGTTT'], ['TTTCCCCCCAAACCCTTTGGGGGGGGGAAAGGGTTT', 'TTTCCCCCCAAACCCTTTGGGGGGGGGAAAGGGTTT', 'TTTCCCCCCAAACCCTTTGGGGGGGGGAAAGGGTTT', 'TTTCCCCCCAAACCCTTTGGGGGGGGGAAAGGGTTT', 'TTTCCCCCCAAACCCTTTGGGGGGGGGAAAGGGTTT', 'TTTCCCCCCAAACCCTTTGGGGGGGGGAAAGGGTTT', 'TTTCCCCCCAAACCCTTTGGGGGGGGGAAAGGGTTT'], ['TTTCCCCCCAAACCCTTTGGGGGGGGGAAAGGGTTT', 'TTTCCCCCCAAACCCTTTGGGGGGGGGAAAGGGTTT', 'TTTCCCCCCAAACCCTTTGGGGGGGGGAAAGGGTTT', 'TTTCCCCCCAAACCCTTTGGGGGGGGGAAAGGGTTT', 'TTTCCCCCCAAACCCTTTGGGGGGGGGAAAGGGTTT']]
            # 39 [['AAAGGGTTTAAACCCAAAGGGTTTGGGGGGAAAAAA', 'AAAGGGTTTAAACCCAAAGGGTTTGGGGGGAAAAAA', 'AAAGGGTTTAAACCCAAAGGGTTTGGGGGGAAAAAA']]

            for ii, u in enumerate(trees_np):
                for jj, v in enumerate(u):
                    if v != '0':
                        trees_ori[ii].append(trees_np[ii][jj])
            ### k, trees_ori
            # 23 [['1', '2', '5', '6', '8'], ['1', '3', '4', '8', '9', '10'], ['1', '4', '7', '9'], ['2', '3', '4', '5', '7', '9'], ['1', '3', '6'], ['1', '2', '4', '7', '8']]
            # 6 [['3', '4', '5', '8'], ['2', '3', '4', '5', '6', '10'], ['5', '6', '7'], ['1', '2', '4', '6', '7', '9'], ['1', '2', '9']]
            # 13 [['1', '2', '5', '7'], ['1', '4', '7', '8', '10'], ['1', '4', '8', '9'], ['6', '7', '9']]
            # ...
            # 29 [['2', '5', '9', '10'], ['1', '2', '3', '5', '6', '8', '10'], ['3', '4', '7', '9', '10']]
            # 39 [['6', '7', '9']]

            for i, tree in enumerate(trees_ori):
                res_data.append([
                    read_cache_table[i][-1] if read_cache_table[i] != [] else read,  # read
                    str(k) + '_' + '_'.join(tree) if read_cache_table[i] != [] else str(k),  # sample id
                    spl_id_map[str(k) + '_' + '_'.join(tree)] if read_cache_table[i] != [] else 'init',  # source
                ])
            # print('res_data {}'.format(res_data))
        # print(res_data)
        # print(len(res_data))
        return np.array(res_data)

    def findListDuplicates(
            self,
            l
    ):
        seen = set()
        seen_add = seen.add
        seen_twice = set(x for x in l if x in seen or seen_add(x))
        return list(seen_twice)