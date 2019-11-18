import os
import numpy as np
from PIL import Image
import openslide
from openslide import open_slide
from openslide.deepzoom import DeepZoomGenerator
import xml.etree.ElementTree as ET
from xml.dom import minidom
from shapely.geometry import Polygon, Point

#author @t_sanf, modifying code from py_wsi by @ysbecca

class MakeBox:

    def __init__(self):
        self.filedir = '/Users/sanforth/Desktop/files_xml'
        self.save_loc='/Users/sanforth/Desktop/save_tiles'

    def sample_and_store_patches(self,
                                 file_name='TCGA-HQ-A2OF-01A-01-TSA_died.svs',
                                 file_dir='/Users/sanforth/Desktop/files_xml/',
                                 pixel_overlap=0,
                                 patch_size=400,
                                 xml_dir=True,
                                 label_map={'1':1,'2':1},
                                 limit_bounds=True,
                                 rows_per_txn=20,
                                 db_location='/Users/sanforth/Desktop/db/',
                                 ):
        ''' Sample patches of specified size from .svs file.
            - file_name             name of whole slide image to sample from
            - file_dir              directory file is located in
            - pixel_overlap         pixels overlap on each side
            - level                 0 is lowest resolution; level_count - 1 is highest
            - xml_dir               directory containing annotation XML files
            - label_map             dictionary mapping string labels to integers
            - rows_per_txn          how many patches to load into memory at once
            - storage_option        the patch storage option

            Note: patch_size is the dimension of the sampled patches, NOT equivalent to openslide's definition
            of tile_size. This implementation was chosen to allow for more intuitive usage.
        '''

        tile_size = self.patch_to_tile_size(patch_size, pixel_overlap)
        slide = open_slide(os.path.join(file_dir,file_name))
        mag_0 = int(10.0 / float(slide.properties[openslide.PROPERTY_NAME_MPP_X])) #get magnification of - level
        tiles = DeepZoomGenerator(slide,tile_size=tile_size,overlap=pixel_overlap,limit_bounds=limit_bounds)
        tiles_10 = DeepZoomGenerator(slide, tile_size=(tile_size), overlap=pixel_overlap, limit_bounds=limit_bounds)
        tiles_5 = DeepZoomGenerator(slide, tile_size=(tile_size), overlap=pixel_overlap, limit_bounds=limit_bounds)

        # Make sure to select the 20x dimension as base level
        if 39<mag_0<41: level_20=tiles.level_count-2
        if 19<mag_0<21: level_20=tiles.level_count-1 #note - minus 1 due to 0 indexing

        # Select 5x and 10x levels
        level_dimension=tiles.level_dimensions
        dim_20=level_dimension[level_20]
        n_level=0; level_10=0; level_5=0
        for dim in level_dimension:
            n_level+=1
            if 3.5<(dim_20[0]/dim[0])<4.5: level_10=n_level
            if 15.5<(dim_20[0] / dim[0])<16.5: level_5=n_level

        level_dict={'20':(tiles,level_20),'10':(tiles_10,level_10),'5':(tiles_5,level_5)}
        print(dict)

        if xml_dir:
            # Expect filename of XML annotations to match SVS file name
            regions, region_labels = self.get_regions(os.path.join(file_dir,file_name[:-4] + ".xml"))

        if level_20 >= tiles.level_count:
            print("[py-wsi error]: requested level does not exist. Number of slide levels: " + str(tiles.level_count))
            return 0
        x_tiles, y_tiles = tiles.level_tiles[level_20]
        print("num x tiles is {}".format(x_tiles))

        x, y = 0, 0
        count, batch_count = 0, 0
        patches, coords, labels = [], [], []

        for val in level_dict.values():
            print(val)
            while y < y_tiles:
                while x < x_tiles:

                    new_tile= np.array(val[0].get_tile(val[1], (x, y)), dtype=np.uint8)

                    # OpenSlide calculates overlap in such a way that sometimes depending on the dimensions, edge
                    # patches are smaller than the others. We will ignore such patches.
                    if np.shape(new_tile) == (patch_size, patch_size, 3):
                        print((new_tile.shape))
                        print(x*tile_size); print(y*tile_size)

                        patches.append(new_tile)
                        coords.append(np.array([x, y]))
                        count += 1

                        # Calculate the patch label based on centre point.
                        if xml_dir:
                            converted_coords = tiles.get_tile_coordinates(val[1], (x, y))[0]
                            labels.append(self.generate_label(regions, region_labels, converted_coords, label_map))
                    x += 1

                # To save memory, we will save data into the dbs every rows_per_txn rows. i.e., each transaction will commit
                # rows_per_txn rows of patches. Write after last row regardless.
                if (y % rows_per_txn == 0 and y != 0) or y == y_tiles-1:
                    self.save_to_disk(db_location, patches, coords, file_name[:-4], labels)

                y += 1
                x = 0

        return count


    def get_regions(self, path):
        ''' Parses the xml at the given path, assuming annotation format importable by ImageScope. '''

        xml = minidom.parse(path)

        # The first region marked is always the tumour delineation
        regions_ = xml.getElementsByTagName("Region")
        regions, region_labels = [], []
        reg_num=1
        for region in regions_:
            r_label = str(reg_num)
            region_labels.append(r_label)
            vertices = region.getElementsByTagName("Vertex")

            # Store x, y coordinates into a 2D array in format [x1, y1], [x2, y2], ...
            coords = np.zeros((len(vertices), 2))

            for i, vertex in enumerate(vertices):
                coords[i][0] = vertex.attributes['X'].value
                coords[i][1] = vertex.attributes['Y'].value
            reg_num+=1
            regions.append(coords)
        return regions, region_labels

    def save_to_disk(self, db_location, patches, coords, file_name, labels):
        """ Saves numpy patches to .png files (full resolution).
            Meta data is saved in the file name.
            - db_location       folder to save images in
            - patches           numpy images
            - coords            x, y tile coordinates
            - file_name         original source WSI name
            - labels            patch labels (opt)
        """
        save_labels = len(labels)
        for i, patch in enumerate(patches):
            # Construct the new PNG filename
            patch_fname = file_name + "_" + str(coords[i][0]) + "_" + str(coords[i][1]) + "_"

            if save_labels:
                patch_fname += str(labels[i])

            # Save the image.
            Image.fromarray(patch).save(db_location + patch_fname + ".png")

    def generate_label(self,regions, region_labels, point, label_map):
        ''' Generates a label given an array of regions.
            - regions               array of vertices
            - region_labels         corresponding labels for the regions
            - point                 x, y tuple
            - label_map             the label dictionary mapping string labels to integer labels
        '''
        for i in range(len(region_labels)):
            poly = Polygon(regions[i])
            if poly.contains(Point(point[0], point[1])):
                return 1
            else:
                return 0


    def patch_to_tile_size(self,patch_size, overlap):
        return patch_size - overlap * 2

    def check_label_exists(self,label, label_map):
        ''' Checking if a label is a valid label.
        '''
        if label in label_map:
            return True
        else:
            print("py_wsi error: provided label " + str(label) + " not present in label map.")
            print("Setting label as -1 for UNRECOGNISED LABEL.")
            print(label_map)
            return False

    def whitespace_check(self,im):
        bw = im.convert('L')
        bw = np.array(bw)
        bw = bw.astype('float')
        bw=bw/255
        prop_ws = (bw > 0.8).sum()/(bw>0).sum()
        return prop_ws





if __name__=='__main__':
    MakeBox().sample_and_store_patches()
