import sys
import os
import numpy as np
from PIL import Image
from PIL import ImageOps
import cv2
import openslide
from openslide import open_slide
from openslide.deepzoom import DeepZoomGenerator
import xml.etree.ElementTree as ET
from xml.dom import minidom
import pandas as pd
from skimage import draw
import matplotlib.pyplot as plt
from shapely.geometry import Polygon, Point, MultiPoint, MultiPolygon, box
import argparse
from utils import *

class makeMap:

    def __init__(self):
        # self.image_locations = args.image_locations
        # self.has_xml = args.has_xml
        # self.save_location = args.save_location
        # self.master_pred = args.master_pred
        # self.coded_filenames = args.coded_filenames
        self.image_locations = '/path/to/svs/imgs'
        self.has_xml = True
        self.save_location = 'path/to/test_maps'
        self.master_pred = 'path/to/prediction_file_example.csv'
        # THIS PREDICTION FILE IS A CSV WITH THE FOLLOWING FORMAT
        # mag, output/file/location/
        # example for 3 res:
        # 10.0, /my/output/file/for/10x/inference.csv
        # 5.0, /my/output/file/for/5x/inference.csv
        # 2.5, /my/output/file/for/2.5x/inference.csv
        self.coded_filenames = True # for TCGA work this is FALSE because data was already anon
        # this is not working right now, will produce ground truth map of annotations
        self.stride_ratio = 2
        self.writemask = False
        self.simplifyroi = True  # default = False, set to true for complex shapes
        self.nolabel = False  # if all regions in an annotation file belong to the same class, they are labeled as 'tumor'
        #      nolabel=FALSE should only be used if the "Text" attribute in xml corresponds to label

    def parsePredictions(self):
        all_preds = self.read_preds()
        if self.coded_filenames == True:
            imglist = findCodedFiles(image_locations=self.image_locations,has_xml=self.has_xml)
        else:
            imglist = findImageFiles(image_locations=self.image_locations, has_xml=self.has_xml)

        #find the base layer and define output map from here
        lowest_mag = max(list(all_preds.keys())) #every base mag prediction will fill ONE voxel
        lowest_preds = all_preds[lowest_mag]
        other_levels = dict(filter(lambda elem: elem[0] != lowest_mag, all_preds.items()))

        for pt in list(set(lowest_preds['patient'])):
            print(pt)
            pt_preds = lowest_preds[lowest_preds['patient'] == pt]

            for block in list(set(pt_preds['block'])):
                print(block)
                block_preds = pt_preds[pt_preds['block'] == block]
                if self.coded_filenames == True:
                    image_name = pt+'_'+block
                else:
                    image_name = block # THIS SHOULD BE CHANGED DEPENDING ON HOW YOU DID THINGS

                #find actual svs image file
                block_img = imglist[imglist['savename']==image_name]
                # find the image size that was actually extracted on prediction files
                phys_size = max([int(i) for i in list(set(block_preds['size']))[0].split('-')])

                # now we utilize the original image to create a probability map
                # first grab data from digital header
                oslide = openslide.OpenSlide(list(block_img['image'])[0])
                # this is physical microns per pixel
                acq_mag = 10.0 / float(oslide.properties[openslide.PROPERTY_NAME_MPP_X])
                # this is nearest multiple of 20 for base layer
                base_mag = int(20 * round(float(acq_mag) / 20))
                base_dim = oslide.dimensions

                base_img = self.make_level_img(base_dim=base_dim, phys_size=phys_size, block_preds=block_preds, image_name=image_name, level_name = lowest_mag)

                for lvl in list(other_levels.keys()):
                    lvl_preds = other_levels[lvl]
                    lvl_pt = lvl_preds[lvl_preds['patient'] == pt]
                    lvl_block = lvl_pt[lvl_pt['block'] == block]
                    lvl_img = self.make_level_img(base_dim=base_dim, phys_size=phys_size, block_preds=lvl_block, image_name=image_name, level_name = lvl)
                    base_img +=lvl_img

                base_img = base_img/3
                slideimg = Image.fromarray(np.uint8(base_img * 255))
                slideimg = slideimg.convert('L')
                slideimg.save(os.path.join(self.save_location, image_name + '_avg.jpeg'))
        return all_preds, imglist

    def read_preds(self):
        ''' create a dictionary of predictions from all magnifications provided in master csv file '''
        df_csv = pd.read_csv(self.master_pred,sep=',',header=0)
        all_preds = {}
        for index, pred_i in df_csv.iterrows():
            mag = pred_i['mag']
            print(str(mag)+' - '+pred_i['file'])
            df_i = pd.read_csv(pred_i['file'],sep=',',header=0)
            # YOU WILL NEED TO CHANGE THIS IF YOU CHANGE YOUR NAMING CONVENTION
            # df_i['RP'], df_i['block'], df_i['mag'], df_i['loc'], df_i['size'], df_i['ws'], df_i['label'] = df_i[
            #     'img_name'].str.split('_', expand=True)
            new_df = df_i['img_name'].str.split('_', expand=True)
            new_df = new_df.rename(columns={0: 'patient', 1: 'block', 2: 'mag', 3: 'loc', 4: 'size', 5: 'ws', 6: 'label'})
            df_i = pd.merge(df_i,new_df,left_index=True,right_index=True)
            all_preds[mag] = df_i
        return all_preds

    def make_level_img(self,base_dim,phys_size,block_preds,image_name,level_name):
        x = np.zeros((round(base_dim[1]*self.stride_ratio/phys_size),round(base_dim[0]*self.stride_ratio/phys_size)), np.float)
        x_count = np.ones((round(base_dim[1]*self.stride_ratio/phys_size),round(base_dim[0]*self.stride_ratio/phys_size)), np.float)
        x_mask = np.zeros((round(base_dim[1] * self.stride_ratio / phys_size),round(base_dim[0] * self.stride_ratio / phys_size)), np.float)

        for index, patch in block_preds.iterrows():
            patch_start = [round(int(i)*self.stride_ratio/phys_size) for i in patch['loc'].split('-')]
            x_count[patch_start[1]:patch_start[1]+self.stride_ratio,patch_start[0]:patch_start[0]+self.stride_ratio] += 1
            x_mask[patch_start[1]:patch_start[1] + self.stride_ratio,patch_start[0]:patch_start[0] + self.stride_ratio] = 1
            x[patch_start[1]:patch_start[1]+self.stride_ratio,patch_start[0]:patch_start[0]+self.stride_ratio] += patch['HR']

        x = x/x_count
        slideimg = Image.fromarray(np.uint8(x*255))
        slideimg = slideimg.convert('L')
        slideimg.save(os.path.join(self.save_location,image_name+'_'+str(level_name)+'.jpeg'))
        slidemask = Image.fromarray(np.uint8(x_mask*255))
        slidemask = slidemask.convert('L')
        slidemask.save(os.path.join(self.save_location,image_name+'_'+str(level_name)+'_mask.jpeg'))
        return x

if __name__ == '__main__':
    # parser = argparse.ArgumentParser()
    # parser.add_argument('--image_locations')
    # parser.add_argument('--save_location')
    # parser.add_argument('--has_xml')
    # parser.add_argument('--master_pred')
    # parser.add_argument('--coded_filenames')
    # args = parser.parse_args()
    c = makeMap()
    all_preds, imglist = c.parsePredictions()
