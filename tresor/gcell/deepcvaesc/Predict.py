__version__ = "v1.0"
__copyright__ = "Copyright 2024"
__license__ = "MIT"
__developer__ = "Jianfeng Sun"
__maintainer__ = "Jianfeng Sun"
__email__ = "jianfeng.sunmt@gmail.com"
__lab__ = "Cribbslab"

import numpy as np
from tresor.util.plot.DimensionReduction import DimensionReduction as drplot


def sample(
        model_fpn,
        num_cells_per_cluster=100,
        image_size=144,
        num_labels=11,
        batch_size=32,
):
    import tensorflow as tf
    decoder = tf.keras.models.load_model(model_fpn)
    print(decoder)
    xs = []
    ys = []
    for ilabel in range(num_labels):
        print('===>label {}'.format(ilabel))
        digit_size = image_size
        z_combo_spl = np.random.normal(loc=1, scale=4, size=[num_cells_per_cluster, 2])
        # z_combo_spl = np.random.uniform(0, 4, size=[num_, 2])
        ilabel_ = np.eye(num_labels)[np.array([ilabel])]
        for i, zcspl in enumerate(z_combo_spl):
            z_sample = np.array([zcspl])
            x_decoded = decoder.get_layer("decoder").predict([z_sample, ilabel_])
            # print(x_decoded.shape)
            # print(x_decoded[0].shape)
            digit = x_decoded[0].reshape(digit_size, digit_size)
            # print(digit)
            xs.append(digit.reshape(digit_size * digit_size))
            t = digit.reshape(digit_size * digit_size)*260
            t = np.floor(t).astype(int)
            print(t[t>1])
            print(t[t>1].shape)
            # xs.append(x_decoded[0])
            ys.append(ilabel)
            # sds = digit * 255
            # print(np.where(sds > 1, sds, 0))
    drplot().single(X=xs, y=ys, tech='TSNE', marker_size=3, cmap='tab20b', title='CondiCVAE')
    drplot().single(X=xs, y=ys, tech='PCA', marker_size=3, cmap='tab20b', title='CondiCVAE')
    drplot().single(X=xs, y=ys, tech='UMAP', marker_size=3, cmap='tab20b', title='CondiCVAE')