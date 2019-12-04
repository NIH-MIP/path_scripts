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
from shapely.geometry import Polygon, Point, MultiPoint, MultiPolygon

class extractPatch:

    def __init__(self):
        self.file_location = 'V:/TCGABLCA/cohort_final/TCGA-2F-A9KQ'
        self.image_file = 'TCGA-2F-A9KQ-01Z-00-DX1.1C8CB2DD-5CC6-4E99-A0F9-32A0F598F5F9.svs'
        self.xml_file = 'TCGA-2F-A9KQ-01Z-00-DX1.1C8CB2DD-5CC6-4E99-A0F9-32A0F598F5F9.xml' #if you dont have an xml file specify 'none'
        self.save_location = 'V:/dummy_save'
        self.save_name = 'TCGA-2F-A9KQ-01Z-00-DX1'
        self.mag_extract = [5,10,20] # specify which magnifications you wish to pull images from
        self.save_image_size = 300   # specify image size to be saved (note this is the same for all magnifications)
        self.pixel_overlap = 0       # specify the level of pixel overlap in your saved images
        self.limit_bounds = True     # this is weird, dont change it
        self.write_all = False       # default is to only write patches that overlap with xml regions (if no xml provided, all patches written)
        self.nolabel = True          # if all regions in an annotation file belong to the same class, they are labeled as 'tumor'
                                     #      nolabel=FALSE should only be used if the "Text" attribute in xml corresponds to label

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
        tiles = DeepZoomGenerator(oslide, tile_size=physSize-round(self.pixel_overlap*acq_mag/base_mag), overlap=round(self.pixel_overlap*acq_mag/base_mag/2), limit_bounds=self.limit_bounds)

        # calculate the effective magnification at each level of tiles, determined from base magnification
        tile_lvls = tuple(base_mag/(tiles._l_z_downsamples[i]*tiles._l0_l_downsamples[tiles._slide_from_dz_level[i]]) for i in range(0,tiles.level_count))

        #only do this if xml exists, if it does spit a warning that all tiles will be saved
        if os.path.exists(os.path.join(self.file_location,self.xml_file)):
            # pull xml data - remember vertices always correspond to lowest level
            regions, region_labels = self.get_regions(os.path.join(self.file_location,self.xml_file),nolabel=self.nolabel)
        else:
            print("WARNING! xml file not found. extracting all tiles and assigning label = NA")

        # pull tiles from levels specified by self.mag_extract
        for lvl in self.mag_extract:
            if lvl in tile_lvls:
                if not os.path.exists(os.path.join(self.save_location, str(lvl) + "x")):
                    os.mkdir(os.path.join(self.save_location, str(lvl) + "x"))
                # pull tile info for level
                x_tiles, y_tiles = tiles.level_tiles[tile_lvls.index(lvl)]

                # note to self, we have to iterate b/c deepzoom does not allow casting all at once at list (??)
                for y in range(0,y_tiles):
                    for x in range(0,x_tiles):

                        # grab tile coordinates
                        tile_coords = tiles.get_tile_coordinates(tile_lvls.index(lvl), (x, y))
                        save_coords = str(tile_coords[0][0]) + "-" + str(tile_coords[0][1]) + "_" + '%.0f'%(tiles._l0_l_downsamples[tile_coords[1]]*tile_coords[2][0]) + "-" + '%.0f'%(tiles._l0_l_downsamples[tile_coords[1]]*tile_coords[2][1])

                        # label tile based on xml region membership
                        if os.path.exists(os.path.join(self.file_location, self.xml_file)):
                            tile_ends = (int(tile_coords[0][0] + tiles._l0_l_downsamples[tile_coords[1]] * tile_coords[2][0]),int(tile_coords[0][1] + tiles._l0_l_downsamples[tile_coords[1]] * tile_coords[2][1]))
                            tile_label = self.assign_label(tile_starts=tile_coords[0],tile_ends=tile_ends,regions=regions,region_labels=region_labels)
                        else:
                            tile_label = {}
                            tile_label['NA']=1

                        if self.write_all == True:
                            if len(tile_label) == 0:
                                tile_label['NA'] = 1
                            tile_size = tiles.get_tile_dimensions(tile_lvls.index(lvl), (x, y))
                            tile_pull = tiles.get_tile(tile_lvls.index(lvl), (x, y))
                            self.save_to_disk(tile_pull=tile_pull, save_coords=save_coords, lvl=lvl, tile_label=tile_label, tile_size=tile_size, physSize=physSize)
                        else:
                            if len(tile_label)>0:
                                tile_size = tiles.get_tile_dimensions(tile_lvls.index(lvl), (x, y))
                                tile_pull = tiles.get_tile(tile_lvls.index(lvl), (x, y))
                                self.save_to_disk(tile_pull=tile_pull, save_coords=save_coords, lvl=lvl, tile_label=tile_label, tile_size=tile_size, physSize=physSize)

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

    def get_regions(self, path, nolabel):
        ''' Parses the xml at the given path, assuming annotation format importable by ImageScope. '''
        ''' nolabel=True identifies instances where any annotation region corresponds to single class (i.e. cancer) '''
        ''' nolabel=False assumes labels are given per region per "Text" attribute '''

        xml = minidom.parse(path)

        regions_ = xml.getElementsByTagName("Region")
        regions, region_labels = [], []
        for region in regions_:
            r_label = region.getAttribute('Text')
            region_labels.append(r_label)
            vertices = region.getElementsByTagName("Vertex")

            # Store x, y coordinates into a 2D array in format [x1, y1], [x2, y2], ...
            coords = np.zeros((len(vertices), 2))

            for i, vertex in enumerate(vertices):
                coords[i][0] = vertex.attributes['X'].value
                coords[i][1] = vertex.attributes['Y'].value
            regions.append(coords)

        if(nolabel == True):
            region_labels = [label.replace('', 'tumor') for label in region_labels]

        return regions, region_labels


    def assign_label(self,tile_starts,tile_ends,regions,region_labels):
        ''' calculates overlap of tile with xml regions and creates dictionary based on unique labels '''

        tile_box = [tile_starts[0],tile_starts[1]],[tile_starts[0],tile_ends[1]],[tile_ends[0],tile_starts[1]],[tile_ends[0],tile_ends[1]]
        tile_box = list(tile_box)
        tile_box = MultiPoint(tile_box).convex_hull

        tile_label = {}
        # create a dictionary of label/value pairs: returns percent of tile containing unique
        for label in set(region_labels):

            # grab regions that correspond to this label
            label_list = [i for i, e in enumerate(region_labels) if e == label]
            labels = tuple(regions[i] for i in label_list)

            # loop over every region associated with a given label, sum the overlap
            box_label = False # initialize
            ov = 0 # initialize
            for reg in labels:
                poly = Polygon(reg)
                if poly.is_valid == False:
                    poly = poly.buffer(0)
                poly_label = tile_box.intersects(poly)
                if poly_label == True:
                    box_label = True
                    ov_reg = tile_box.intersection(poly)
                    ov += ov_reg.area/tile_box.area

            if box_label == True:
                tile_label[label] = ov

        # p.s. if you are curious, you can plot the polygons by the following
        #   plt.plot(*poly.exterior.xy) and plt.plot(*tile_box.exterior.xy)
        return tile_label

    def save_to_disk(self, tile_pull, save_coords, lvl, tile_label, tile_size, physSize):
        # edge tiles will not be correct size (too small), so we reflect image data until correct size
        if tile_size != (physSize, physSize):
            tile_pull = Image.fromarray(cv2.copyMakeBorder(np.array(tile_pull), 0, physSize - tile_size[1], 0, physSize - tile_size[0],cv2.BORDER_REFLECT))
        # check whitespace amount
        ws = self.whitespace_check(im=tile_pull)

        label_text = 'label'
        for key, value in tile_label.items():
            label_text += '-' + key + "-" + '%.2f' % (value)

        tile_savename = self.save_name + "_" + str(lvl) + "_" \
                        + save_coords + "_" \
                        + "ws-" + '%.2f' % (ws) \
                        + "_" + label_text
        tile_pull = tile_pull.resize(size=(self.save_image_size, self.save_image_size), resample=Image.ANTIALIAS)
        tile_pull.save(os.path.join(self.save_location, str(lvl) + "x", tile_savename + ".png"))
        # print(tile_savename)
        return


if __name__ == '__main__':
    c = extractPatch()
    c.parseMeta_and_pullTiles()
