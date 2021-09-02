import numpy as np
from scipy import ndimage
import numpy.ma as ma
import math

from matplotlib import pyplot as plt

""" matrice to ground coordinates"""
def matrix_to_ground_coordinates(x,y, nRows,nCols, cellSize):
	#1rst step :correction of axes such that 0 is the center
	i = x - (nRows-1)/2.0
	j = y - (nCols-1)/2.0
	#2nd step : matrix coordinates -> ground coordinates
	X = i * cellSize
	Y = j * cellSize
	return X,Y
    
def point_rotation(i, j, cols, rows, angle) :
   ic = i - (rows -1)/2
   jc = jc = j - (cols -1)/2
   icrot = jc*math.sin(angle) + ic*math.cos(angle)
   jcrot = jc*math.cos(angle) - ic*math.sin(angle)
   irot = icrot + (rows -1)/2
   jrot = jcrot + (cols -1)/2
   return irot,jrot


    
class terrain(object):

    def __init__(self,filename, cellSize = 2.5, noValue = -9999 ):
        self._filename = filename
        if type(filename) == str:
            self._matrix = np.loadtxt(filename,skiprows=6)
            self._orig_matrix = np.copy(self._matrix)
        elif type(filename) == np.ndarray:
            self._matrix = filename
        else :
            raise TypeError
        n_rows, n_cols = np.shape(self._matrix)
        self._n_rows = n_rows
        self._n_cols = n_cols
        self._cell_size = float(cellSize)
        self._noValue = noValue
        self._angle = 0
        self._dimension = 0
        self.dim()
        self.bInfI = 40
        self.bMaxI = 240
        self.bInfJ = 95
        self.bMaxJ = 150
        self._zoom = 1

    def dim(self):
        self._dimension = ((self.n_rows-1)*self._cell_size, (self.n_cols-1)*self._cell_size)

    def zoom(self,zoom):
        if zoom == 0:
            return None
        self._matrix = ndimage.zoom(self._matrix,zoom,order=1)
        n_rows, n_cols = np.shape(self._matrix)
        self._n_rows = n_rows
        self._n_cols = n_cols
        self._cell_size = self._dimension[0]/float((self._n_rows-1))
        self._zoom = zoom

    def change_min(self,val = 1100):
	    self._matrix[self._matrix == self._noValue]=val
	    self._noValue = val

    def plot2d(self, show=True, newfig=True):
        if newfig: plt.figure()
        plt.imshow(self._matrix.T - self._orig_matrix.T)
        if show: plt.show()


    def cut(self,bInfI = 40,bMaxI = 240,bInfJ = 95,bMaxJ = 150):
        self._matrix = self._matrix[bInfI:bMaxI,bInfJ:bMaxJ]
        self.bInfI = bInfI
        self.bMaxI = bMaxI
        self.bInfJ = bInfJ
        self.bMaxJ = bMaxJ
        n_rows, n_cols = np.shape(self._matrix)
        self._n_rows = n_rows
        self._n_cols = n_cols
        self.dim()

    def rotation(self,angle):
        self.change_min()
        self._matrix = ndimage.rotate(self._matrix,angle,mode='constant',cval=self._matrix.min(),order=1)
        n_rows, n_cols = np.shape(self._matrix)
        self._n_rows = n_rows
        self._n_cols = n_cols
        self._angle = angle
        self.dim()
        self.cut()


    def uniform_noise(self, noise_amp=1):
        bruit = np.random.uniform(0,noise_amp,np.shape(self._matrix))
        self._matrix =  self._matrix + bruit

    def sinusoidal_noise(self, noise_amp=1, noise_freq=1, phase = 0):
        i,j = np.shape(self._matrix)
        b,c = np.meshgrid(np.linspace(0,self._dimension[0], i),
                          np.linspace(0,self._dimension[1], j))
        bruit1 = noise_amp*np.cos(b*noise_freq*2*math.pi+phase)
        bruit2 = noise_amp*np.sin(c*noise_freq*2*math.pi+phase)
        self._matrix =  self._matrix + bruit1.T + bruit2.T

    def initial_position(self,filenamePosInit):

        if type(filenamePosInit) == np.ndarray or type(filenamePosInit) == list:
            assert len(filenamePosInit)==3
            return filenamePosInit[0],filenamePosInit[1],filenamePosInit[2]
        
        posInit = np.loadtxt(filenamePosInit,skiprows=6)
        rowsPI,colsPI = np.shape(posInit)
        idxFlat = np.argmax(posInit) # By default, the index is into the flattened array, otherwise along the specified axis.
        (i,j)=np.unravel_index(idxFlat, (rowsPI,colsPI))
        val = posInit[i,j]
        if self._angle != 0:
            angle = -(self._angle)*math.pi/180
            irotsg, jrotsg = point_rotation (0, 0, colsPI, rowsPI, angle)
            irotig, jrotig = point_rotation (rowsPI, 0, colsPI, rowsPI, angle)
            irot, jrot = point_rotation (i, j, colsPI, rowsPI, angle)
            icut = irot -(irotsg+self.bInfI)
            jcut = jrot - (jrotig +self.bInfJ)
            X,Y = matrix_to_ground_coordinates(icut*self._zoom,jcut*self._zoom,self._n_rows,self._n_cols, self._cell_size)
            Z = self._matrix[int(icut)*self._zoom,int(jcut)*self._zoom] + val
        else:
            X,Y = matrix_to_ground_coordinates (i*self._zoom,j*self._zoom, self._n_rows,self._n_cols, self._cell_size)
            Z = np.loadtxt(self._filename,skiprows=6)[i,j] + val
        return X,Y,Z


    def _get_matrix(self):
        return self._matrix

    def _get_n_rows(self):
        return self._n_rows

    def _get_n_cols(self):
        return self._n_cols

    def _get_cell_size(self):
        return self._cell_size

    def _get_noValue(self):
        return self._noValue

    def _get_dimension(self):
        return self._dimension

    matrix = property(_get_matrix)
    n_rows = property(_get_n_rows)
    n_cols = property(_get_n_cols)
    cellSize  = property(_get_cell_size)
    noValue = property(_get_noValue)
    dimension = property(_get_dimension)




if __name__ == "__main__":
    Terrain = terrain('./data/dem.asc')
    print(Terrain)
    print(Terrain.matrix)
    print(Terrain.dimension)
