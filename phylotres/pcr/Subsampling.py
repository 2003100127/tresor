__version__ = "v1.0"
__copyright__ = "Copyright 2024"
__license__ = "MIT"
__developer__ = "Jianfeng Sun"
__maintainer__ = "Jianfeng Sun"
__email__="jianfeng.sunmt@gmail.com"
__lab__ = "Cribbslab"

import numpy as np
from phylotres.pcr.Amplify import Amplify as pcr
from phylotres.seq.Calling import Calling as seq
from phylotres.util.sequence.fastq.Write import write as wfastq
from phylotres.util.file.read.Reader import Reader as pfreader
from phylotres.util.random.Number import Number as rannum
from phylotres.util.sequence.symbol.Single import Single as dnasgl
from phylotres.util.Console import Console
from phylotres.util.Kit import tactic6


class Subsampling:

    def __init__(
            self,

            # verbose=False
    ):
        self.pcr = pcr
        self.seq = seq
        self.wfastq = wfastq
        self.console = Console()
        # self.console.verbose = verbose
    def minnow(self, pcr_dict):
        umi_map = pfreader().generic(pcr_dict['read_lib_fpn'])[0].to_dict()
        # print(umi_map)
        nn = pcr_dict['data'].shape[0]
        spl_ids = rannum().uniform(
            low=0, high=nn, num=250, use_seed=False, seed=1
        )
        # print(spl_ids)
        # print(pcr_dict['data'])
        spl_id_map = tactic6(pcr_dict['data'][:, [1, 2]])
        # print(spl_id_map)
        # print(len(spl_id_map))
        spl_mut_info = pcr_dict['mut_info'][spl_ids]
        print(pcr_dict['mut_info'].shape)
        print(spl_mut_info)
        print(len(spl_mut_info))
        keys = spl_mut_info[:, 2]
        print(len(keys))
        pos_dict = tactic6(pcr_dict['mut_info'][:, [2, 0]])
        base_dict = tactic6(pcr_dict['mut_info'][:, [2, 1]])
        # print(pos_dict)
        # print(base_dict)
        # print(tactic6(base_np))
        res_data = []
        for key in keys:
            mol_id = key.split('_')[0]
            k = key.split('_')[1:]
            print('kkk', key, k)
            read = umi_map[int(mol_id)]
            print(read)
            for i in range(len(k)):
                print('id', i)
                sub_k = mol_id + '_' + '_'.join(k[: i+1]) if k != [] else mol_id
                print(sub_k)
                print(pos_dict[sub_k], base_dict[sub_k])
                read = self.change(read, pos_list=pos_dict[sub_k], base_list=base_dict[sub_k])
            #     print(read)
            # print(read)
            res_data.append([
                read,  # read
                str(mol_id) + '_' + '_'.join(k) if k != [] else str(mol_id),  # sam id
                spl_id_map[str(mol_id) + '_' + '_'.join(k)] if k != [] else 'init',  # source
            ])
            # print(read)
        print(np.array(res_data).shape)
        return np.array(res_data)

    def pcrtree(self, pcr_dict):
        """

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
        umi_map = pfreader().generic(pcr_dict['read_lib_fpn'])[0].to_dict()
        # print(umi_map)

        nn = pcr_dict['data'].shape[0]
        # print(nn)
        spl_ids = rannum().uniform(
            low=0, high=nn, num=9000, use_seed=False, seed=1
        )
        # print(spl_ids)
        spl_data = pcr_dict['data'][spl_ids]
        # print(spl_data)
        spl_id_map = tactic6(spl_data)
        # print(spl_id_map)

        trees = spl_data[:, 0].ravel().tolist()
        # print(trees)
        # print(len(trees))

        res_data = []
        mol_map_to_all_its_pcr_trees = {}
        for tree in trees:
            mol_map_to_all_its_pcr_trees[tree.split('_')[0]] = []
        for tree in trees:
            mol_map_to_all_its_pcr_trees[tree.split('_')[0]].append(tree.split('_')[1:])
        # print(mol_map_to_all_its_pcr_trees)
        # print(len(mol_map_to_all_its_pcr_trees))

        for k, trees in mol_map_to_all_its_pcr_trees.items():
            # print(k, trees)
            read = umi_map[int(k)]
            # print(read)
            read_cache_table = [[] for _ in range(len(trees))]
            bool_flag_table = [[True] for _ in range(len(trees))]
            # print('asdasdasdsa {}'.format(bool_flag_table))
            realvalued_flag_table = [[] for _ in range(len(trees))]
            trees_ori = [[] for _ in range(len(trees))]
            max_len_pcr_tree = max([len(tree) for tree in trees])
            # print(max_len_pcr_tree)
            for _, tree in enumerate(trees):
                if len(tree) < max_len_pcr_tree:
                    tree += ['0' for _ in range(max_len_pcr_tree - len(tree))]
                # print(k, tree)
            trees_np = np.array(trees)

            # ### /*** block. construct bool and real-valued flag tables ***/
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
            # ### /*** block. construct bool and real-valued flag tables ***/
            for id_horiz_in_a_tree in range(trees_np.shape[1]):
                repeat_in_a_col = self.findListDuplicates(trees_np[:, id_horiz_in_a_tree])
                # print('repeated in a col', repeat_in_a_col)
                # print(trees_np[:, id_horiz_in_a_tree])
                for id_vert_across_trees, ele_in_a_col in enumerate(trees_np[:, id_horiz_in_a_tree]):
                    # print('ele_in_a_col {}'.format(ele_in_a_col))
                    if ele_in_a_col in repeat_in_a_col:
                        ids_repeat_in_a_col = [i for i, value in
                                               enumerate(trees_np[:, id_horiz_in_a_tree]) if
                                               value == ele_in_a_col]
                        # print('asds {}'.format(ids_repeat_in_a_col))
                        if bool_flag_table[id_vert_across_trees][id_horiz_in_a_tree] is False:
                            # print('repeated ele: {}'.format(ele_in_a_col))
                            bool_flag_table[id_vert_across_trees].append(False)
                            realvalued_flag_table[id_vert_across_trees].append(-1)
                        else:
                            # print('non-repeated ele: {}'.format(ele_in_a_col))
                            inspector_flags = [1 if bool_flag_table[i][id_horiz_in_a_tree] is True else 0 for i in
                                               ids_repeat_in_a_col]
                            # print('inspector_flags {}'.format(inspector_flags))
                            if sum(inspector_flags) > 1:
                                bool_flag_table[id_vert_across_trees].append(True)
                                realvalued_flag_table[id_vert_across_trees].append(ele_in_a_col)
                            else:
                                bool_flag_table[id_vert_across_trees].append(False)
                                realvalued_flag_table[id_vert_across_trees].append(-1)
                    else:
                        bool_flag_table[id_vert_across_trees].append(False)
                        realvalued_flag_table[id_vert_across_trees].append(-1)

            bool_flag_table = np.array(bool_flag_table)[:, 1:]
            realvalued_flag_table = np.array(realvalued_flag_table)
            # print(trees_np)
            # print(bool_flag_table)
            # print(realvalued_flag_table)

            for jj in range(trees_np.shape[1]):
                read_for_repeat_tmp_per_col = {}
                for ii, val_in_a_col in enumerate(trees_np[:, jj]):
                    if val_in_a_col != '0':
                        if jj == 0:
                            # print(val_in_a_col, bool_flag_table[ii][jj])
                            if bool_flag_table[ii][jj] == True:
                                if val_in_a_col not in [*read_for_repeat_tmp_per_col.keys()]:
                                    r1 = self.mutated(
                                        read=read,
                                        pcr_error=pcr_dict['pcr_error'],
                                    )
                                    read_for_repeat_tmp_per_col[val_in_a_col] = r1
                                    read_cache_table[ii].append(r1)
                                else:
                                    read_cache_table[ii].append(read_for_repeat_tmp_per_col[val_in_a_col])
                            else:
                                r1 = self.mutated(
                                    read=read,
                                    pcr_error=pcr_dict['pcr_error'],
                                )
                                read_cache_table[ii].append(r1)
                        if jj > 0:
                            if bool_flag_table[ii][jj] == True:
                                if val_in_a_col + '_' + '_'.join(list(trees_np[ii][:jj])) not in [*read_for_repeat_tmp_per_col.keys()]:
                                    r1 = self.mutated(
                                        read=read_cache_table[ii][jj - 1],
                                        pcr_error=pcr_dict['pcr_error'],
                                    )
                                    read_for_repeat_tmp_per_col[val_in_a_col + '_' + '_'.join(list(trees_np[ii][:jj]))] = r1
                                    read_cache_table[ii].append(r1)
                                else:
                                    read_cache_table[ii].append(read_for_repeat_tmp_per_col[val_in_a_col + '_' + '_'.join(list(trees_np[ii][:jj]))])
                                # print('id', [i for i, value in enumerate(trees_np[:, jj]) if value == val_in_a_col])
                            else:
                                r1 = self.mutated(
                                    read=read_cache_table[ii][jj - 1],
                                    pcr_error=pcr_dict['pcr_error'],
                                )
                                read_cache_table[ii].append(r1)
            print('read_cache_table {}'.format(read_cache_table))

            for ii, u in enumerate(trees_np):
                for jj, v in enumerate(u):
                    if v != '0':
                        trees_ori[ii].append(trees_np[ii][jj])
            # print(trees_ori)
            print('trees_ori {}'.format(trees_ori))


            for i, tree in enumerate(trees_ori):
                # print(read_cache_table[i])
                # print(read_cache_table[i][-1])
                res_data.append([
                    read_cache_table[i][-1] if read_cache_table[i] != [] else read,  # read
                    str(k) + '_' + '_'.join(tree) if read_cache_table[i] != [] else str(k),  # sam id
                    spl_id_map[str(k) + '_' + '_'.join(tree)] if read_cache_table[i] != [] else 'init',  # source
                ])
            # print('res_data {}'.format(res_data))
        # print(res_data)
        # print(len(res_data))
        return np.array(res_data)

    def change(self, read, pos_list, base_list):
        read_l = list(read)
        # print(read_l)
        for i, pos in enumerate(pos_list):
            dna_map = dnasgl().todict(
                nucleotides=dnasgl().getEleTrimmed(
                    ele_loo=read_l[pos],
                    universal=True,
                ),
                reverse=True,
            )
            read_l[pos] = dna_map[base_list[i]]
        return ''.join(read_l)

    def mutated(self, read, pcr_error):
        num_err_per_read = rannum().binomial(
            n=len(read), p=pcr_error, use_seed=False, seed=False
        )
        pos_list = rannum().uniform(
            low=0, high=len(read), num=num_err_per_read, use_seed=False, seed=False
        )
        base_list = rannum().uniform(
            low=0, high=3, num=num_err_per_read, use_seed=False
        )
        read_l = list(read)
        for i, pos in enumerate(pos_list):
            dna_map = dnasgl().todict(
                nucleotides=dnasgl().getEleTrimmed(
                    ele_loo=read_l[pos],
                    universal=True,
                ),
                reverse=True,
            )
            read_l[pos] = dna_map[base_list[i]]
        return ''.join(read_l)

    def findListDuplicates(self, l):
        seen = set()
        seen_add = seen.add
        seen_twice = set(x for x in l if x in seen or seen_add(x))
        return list(seen_twice)