import sys
import os
import pandas as pd


def findCodedFiles(image_locations,has_xml):
    dirlist = sorted(os.listdir(image_locations))
    data = []
    for subi in dirlist:
        if subi != 'Thumbs.db':
            subsearch = sorted(os.listdir(os.path.join(image_locations, subi)))
            block = 1
            RPname = [i for i in subsearch if i.startswith("RP")]

            if has_xml == True:
                sublist = [i for i in subsearch if i.endswith(".xml")]
                for subli in sublist:
                    if (os.path.exists(
                            os.path.join(image_locations, subi, subli.replace(".xml", ".svs")))):
                        xml_file = os.path.join(image_locations, subi, subli)
                        image_file = os.path.join(image_locations, subi, subli.replace(".xml", ".svs"))
                        save_name = RPname[0].replace(".txt", "") + "_block" + str(block)
                        data.append([save_name,image_file,xml_file])
                        block += 1

            else:
                sublist = [i for i in subsearch if i.endswith(".svs")]
                for subli in sublist:
                    xml_file = 'NA'
                    image_file = os.path.join(image_locations, subi, subli)
                    save_name = RPname[0].replace(".txt", "") + "_block" + str(block)
                    data.append([save_name,image_file,xml_file])
                    block += 1
    col_names = ['savename','image','xml']
    imglist = pd.DataFrame(data, columns = col_names)
    return imglist


def findImageFiles(image_locations,has_xml):
    dirlist = sorted(os.listdir(image_locations))
    data = []
    for subi in dirlist:
        if subi != 'Thumbs.db':
            subsearch = sorted(os.listdir(os.path.join(image_locations, subi)))

            if has_xml == True:
                sublist = [i for i in subsearch if i.endswith(".xml")]
                for subli in sublist:
                    if (os.path.exists(
                            os.path.join(image_locations, subi, subli.replace(".xml", ".svs")))):
                        xml_file = os.path.join(image_locations, subi, subli)
                        image_file = os.path.join(image_locations, subi, subli.replace(".xml", ".svs"))
                        save_name = subli.replace(".xml", "") #this is the same savename initialized in wsi_region_extractor, so if you changed it there you'll need to change it here
                        data.append([save_name,image_file,xml_file])

            else:
                sublist = [i for i in subsearch if i.endswith(".svs")]
                for subli in sublist:
                    xml_file = 'NA'
                    image_file = os.path.join(image_locations, subi, subli)
                    save_name = subli.replace(".xml", "") #this is the same savename initialized in wsi_region_extractor, so if you changed it there you'll need to change it here
                    data.append([save_name,image_file,xml_file])

    col_names = ['savename','image','xml']
    imglist = pd.DataFrame(data, columns = col_names)
    return imglist