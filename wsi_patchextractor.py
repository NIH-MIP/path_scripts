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

class extractPatch:

    def __init__(self):
        self.file_location = 'V:/TCGABLCA/cohort_final/TCGA-2F-A9KQ'
        self.image_file = 'TCGA-2F-A9KQ-01Z-00-DX1.1C8CB2DD-5CC6-4E99-A0F9-32A0F598F5F9.svs'
        self.xml_file = 'TCGA-2F-A9KQ-01Z-00-DX1.1C8CB2DD-5CC6-4E99-A0F9-32A0F598F5F9.xml'
        self.save_location = 'V:/dummy_save'
        self.save_name = 'TCGA-2F-A9KQ-01Z-00-DX1'
        self.mag_extract = [5,10,20] # specify which magnifications you wish to pull images from
        self.save_image_size = 300   # specify image size to be saved (note this is the same for all magnifications)
        self.pixel_overlap = 0       # specify the level of pixel overlap, note this will operate at highest magnification
        self.limit_bounds = True     # this is weird, dont change it

    def parseMeta_and_pullTiles(self):
        # first grab data from digital header
        oslide = openslide.OpenSlide(os.path.join(self.file_location,self.image_file))

        # this is physical microns per pixel
        acq_mag = 10.0/float(oslide.properties[openslide.PROPERTY_NAME_MPP_X])

        # this is nearest multiple of 20 for base layer
        base_mag = int(20 * round(float(acq_mag) / 20))

        # this is how much we need to resample our physical patches for uniformity across studies
        physSize = round(self.save_image_size*acq_mag/base_mag)

        # grab tiles accounting for the physical size we need to pull for standardized tile size across studies
        tiles = DeepZoomGenerator(oslide, tile_size=physSize, overlap=self.pixel_overlap, limit_bounds=self.limit_bounds)

        # calculate the effective magnification at each level of tiles, determined from base magnification
        tile_lvls = tuple(base_mag/(tiles._l_z_downsamples[i]*tiles._l0_l_downsamples[tiles._slide_from_dz_level[i]]) for i in range(0,tiles.level_count))

        # pull tiles from levels specified by self.mag_extract
        for lvl in self.mag_extract:
            if lvl in tile_lvls:
                if not os.path.exists(os.path.join(self.save_location, str(lvl) + "x")):
                    os.mkdir(os.path.join(self.save_location, str(lvl) + "x"))
                x_tiles, y_tiles = tiles.level_tiles[tile_lvls.index(lvl)]
                for y in range(0,y_tiles): #
                    for x in range(0,x_tiles): #
                        tile_coords = tiles.get_tile_coordinates(tile_lvls.index(lvl), (x, y))
                        # note to self, call a function here for labeling based on xml

                        tile_size = tiles.get_tile_dimensions(tile_lvls.index(lvl), (x, y))
                        tile_pull = tiles.get_tile(tile_lvls.index(lvl), (x, y))

                        # check size and pad with data if too small
                        if tile_size != (physSize,physSize):
                            # one way to get desired size is to just pad with zeros using ImageOps
                            #new_tile = ImageOps.expand(tile_pull,(0,0,physSize - tile_size[0],physSize - tile_size[1]))

                            # or we could just reflect data from border to fill image using cv2 (i like this better)
                            tile_pull = Image.fromarray(cv2.copyMakeBorder(np.array(tile_pull),0,physSize-tile_size[1],0,physSize-tile_size[0],cv2.BORDER_REFLECT))

                        # remember we pulled a physical size that allows us to resample to true size
                        tile_pull = tile_pull.resize(size=(self.save_image_size, self.save_image_size), resample=Image.ANTIALIAS)

                        # check whitespace amount
                        ws = self.whitespace_check(im=tile_pull)

                        tile_savename = self.save_name + "_" + str(lvl) + "_" \
                                        + str(tile_coords[0][0]) + "-" + str(tile_coords[0][1]) + "_" \
                                        + '%.0f'%(tiles._l0_l_downsamples[tile_coords[1]]*tile_coords[2][0]) + "-" + '%.0f'%(tiles._l0_l_downsamples[tile_coords[1]]*tile_coords[2][1]) \
                                        + "_ws-" + '%.2f'%(ws)

                        tile_pull.save(os.path.join(self.save_location, str(lvl) + "x", tile_savename + ".png"))
                        # print(tile_savename)

            else:
                print("WARNING: YOU ENTERED AN INCORRECT MAGNIFICATION LEVEL")

        return

    def whitespace_check(self,im):
        bw = im.convert('L')
        bw = np.array(bw)
        bw = bw.astype('float')
        bw=bw/255
        prop_ws = (bw > 0.8).sum()/(bw>0).sum()
        return prop_ws

if __name__ == '__main__':
    c = extractPatch()
    c.parseMeta_and_pullTiles()