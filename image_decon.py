import histomicstk as htk
import numpy as np
import scipy as sp
import skimage.io
import skimage.measure
import skimage.color
import matplotlib.patches as mpatches
from skimage.transform import resize
from matplotlib.colors import ListedColormap
from histomicstk.preprocessing.color_conversion import lab_mean_std
from histomicstk.preprocessing.color_normalization import reinhard
from histomicstk.saliency.tissue_detection import (
    get_slide_thumbnail, get_tissue_mask)
from histomicstk.annotations_and_masks.annotation_and_mask_utils import (
    get_image_from_htk_response)
from histomicstk.preprocessing.color_normalization.\
    deconvolution_based_normalization import deconvolution_based_normalization
from histomicstk.preprocessing.color_deconvolution.\
    color_deconvolution import color_deconvolution_routine, stain_unmixing_routine
from histomicstk.preprocessing.augmentation.\
    color_augmentation import rgb_perturb_stain_concentration, perturb_stain_concentration
import os
import shutil
from PIL import Image


# this is an example function for color normalization using HistomicsTK
# both the color normalization and target image properties are defined from a TCGA sample and are the same as used
#   in Amgad et al 2019 publication for reproducibility 

def colornorm(self,image_location,save_location,img_file):

    # color norm. standard (from TCGA-A2-A3XS-DX1, Amgad et al, 2019)
    cnorm = {
        'mu': np.array([8.74108109, -0.12440419,  0.0444982]),
        'sigma': np.array([0.6135447, 0.10989545, 0.0286032]),
    }

    # TCGA-A2-A3XS-DX1_xmin21421_ymin37486_.png, Amgad et al, 2019)
    # for macenco (obtained using rgb_separate_stains_macenko_pca()
    # and reordered such that columns are the order:
    # Hamtoxylin, Eosin, Null
    W_target = np.array([
        [0.5807549,  0.08314027,  0.08213795],
        [0.71681094,  0.90081588,  0.41999816],
        [0.38588316,  0.42616716, -0.90380025]
    ])

    # specify stains of input image
    stains = ['hematoxylin',  # nuclei stain
              'eosin',  # cytoplasm stain
              'null']  # set to null if input contains only two stains

    imInput = skimage.io.imread(os.path.join(image_location,img_file))[:, :, :3]

    stain_unmixing_routine_params = {
        'stains': ['hematoxylin', 'eosin'],
        'stain_unmixing_method': 'macenko_pca',
    }
    im_norm = deconvolution_based_normalization(
        imInput, W_target=W_target,
        stain_unmixing_routine_params=stain_unmixing_routine_params)

    stain_color_map = htk.preprocessing.color_deconvolution.stain_color_map

    # create stain matrix
    W = np.array([stain_color_map[st] for st in stains]).T

    # perform standard color deconvolution
    imDeconvolved = htk.preprocessing.color_deconvolution.color_deconvolution(im_norm, W)

    norm = Image.fromarray(im_norm)
    norm_bw = norm.convert('L')
    imDeconvolved.Stains[:, :, 2] = np.array(norm_bw)
    decon = Image.fromarray(imDeconvolved.Stains)

    norm.save(os.path.join(save_location+'_norm',img_file))
    norm_bw.save(os.path.join(save_location+'_bw',img_file))
    decon.save(os.path.join(save_location+'_decon',img_file))