__version__ = "0.0.1"
__copyright__ = "Copyright 2025"
__license__ = "GPL-3.0"
__developer__ = "Jianfeng Sun"
__maintainer__ = "Jianfeng Sun"
__email__ = "jianfeng.sunmt@gmail.com"



def tactic6(arr_2d):
    result = {}
    len_arr = len(arr_2d[0])
    if len_arr == 2:
        for item in arr_2d:
            result[item[0]] = item[1]
    else:
        for item in arr_2d:
            result[item[0]] = item[1:]
    return result
