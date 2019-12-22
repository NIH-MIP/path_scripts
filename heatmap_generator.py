import os
import numpy as np
from PIL import Image


class GenHeatmap:

    def __init__(self):
        self.basePATH='/Users/sanforth/Desktop/save_tiles'

    def heatmap(self):
        '''
        trying to make a heatmap at one resoultion first, then combining them should be simple.
        :return:
        '''

        res='5x'
        files=os.listdir(os.path.join(self.basePATH,res))
        size_patch=int(files[0].split('_')[3].split('-')[0])
        us=int(size_patch/16) #us=unit size

        #initialize empty numpy array
        file=self.find_largest_file(file_list=files)
        xy_coord=file.split('_')[2];box_size=file.split('_')[3]
        x_coord=int(xy_coord.split('-')[0])+int(box_size.split('-')[0])
        y_coord = int(xy_coord.split('-')[1]) + int(box_size.split('-')[1])
        x_coord_unit=int(x_coord/us); y_coord_unit=int(y_coord/us)
        array = np.array([range(0,1)] * x_coord_unit*y_coord_unit).reshape((x_coord_unit, y_coord_unit))
        print(array.shape)

        for file in files:
            pred=float(file.split('_')[4].split('-')[1].split('.png')[0])
            pred=pred*255 if pred<1 else pred
            xy_coord=file.split('_')[2]; xy_size=file.split('_')[3]

            #select out the corredinates of the submatrix you are interested in
            topL_coord_x=int(int(xy_coord.split('-')[0])/us)
            topL_coord_y=int(int(xy_coord.split('-')[1])/us)
            botR_coord_x= int(topL_coord_x+ int(xy_size.split('-')[0])/us)
            botR_coord_y=int(topL_coord_y+int(xy_size.split('-')[1])/us)

            #assign the
            array_sm=array[topL_coord_x:botR_coord_x,topL_coord_y:botR_coord_y]

            array[topL_coord_x:botR_coord_x,topL_coord_y:botR_coord_y]=pred

        print(array.shape)
        print(type(array))
        img = Image.fromarray(array,mode='L')
        img.show()
        return array


    def find_largest_file(self,file_list):
        '''
        helper function to find the largest x and y position
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





if __name__=='__main__':
    c=GenHeatmap()
    c.heatmap()

