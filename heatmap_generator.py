import os
import numpy as np
from PIL import Image
from matplotlib import cm


class GenHeatmap:

    def __init__(self):
        self.basePATH='/Users/sanforth/Desktop/save_tiles'

    def multires_heatmap(self):
        '''create multiresolution heatmap
        :param: path(str) path to the
        '''

        path=self.basePATH #right now, only works for one patient.  Will need to add loop for multiple patients

        #create dictionary with heatmap for every resolution
        cor_dict = {'5x': 16, '10x': 8, '20x': 4}
        array_dict={}
        for res in os.listdir(path):
            if res.startswith('.'): #for all the stoopid hidden filenames
                pass
            else:
                correction=cor_dict[res]
                array_dict[res]=self.heatmap(res,correction)

        #combine the multiple resolutions
        array=np.zeros(array_dict['5x'].shape)
        for res in array_dict.keys():
            print(res)
            print(array_dict[res].shape)
            self.display(array_dict[res])
            array+=array_dict[res]

        self.display(array)

        return(array)

    def heatmap(self,res='5x',correction=4):
        '''
        Create a nxm array that assigns predictions to
        :return:
        '''

        #find all files, patch size, and unit size
        files=os.listdir(os.path.join(self.basePATH,res))
        size_patch=int(files[0].split('_')[3].split('-')[0])
        us=int(size_patch/correction) #us=unit size

        #initialize empty numpy array
        file=self.find_largest_file(file_list=files)
        xy_coord=file.split('_')[2]; box_size=file.split('_')[3]
        x_coord=int(xy_coord.split('-')[0])+int(box_size.split('-')[0])
        y_coord = int(xy_coord.split('-')[1]) + int(box_size.split('-')[1])
        array=np.zeros((int(x_coord/us),int(y_coord/us)))

        #iterate over files and
        for file in files:
            pred=float(file.split('_')[4].split('-')[1].split('.png')[0]) #will need to change this line for prediction --> right now just using white space predictions
            xy_coord=file.split('_')[2]; xy_size=file.split('_')[3]

            #select out the corredinates of the submatrix you are interested in
            topL_coord_x=int(int(xy_coord.split('-')[0])/us)
            topL_coord_y=int(int(xy_coord.split('-')[1])/us)
            botR_coord_x= int(topL_coord_x+ int(xy_size.split('-')[0])/us)
            botR_coord_y=int(topL_coord_y+int(xy_size.split('-')[1])/us)

            #assign the prediction to the array
            array[topL_coord_x:botR_coord_x,topL_coord_y:botR_coord_y]=pred

        return array


    def find_largest_file(self,file_list):
        '''
        helper function to find the largest x and y position among filenames in a give directory
        :return:
        '''
        largest_sum=0
        largest_file=''
        for file in file_list:
            xy_coord=file.split('_')[2]
            sum=int(xy_coord.split('-')[0])+int(xy_coord.split('-')[1])
            if sum>largest_sum:
                largest_sum=sum
                largest_file=file
        return largest_file

    def display(self,array):
        '''quick helper function to display image array'''
        im = Image.fromarray(np.uint8(cm.gist_earth(array) * 255))
        im=im.rotate(270,expand=True)
        im=im.transpose(Image.FLIP_LEFT_RIGHT)
        im.show()


if __name__=='__main__':
    c=GenHeatmap()
    c.multires_heatmap()
