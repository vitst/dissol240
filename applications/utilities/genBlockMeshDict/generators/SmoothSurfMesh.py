import genBlockMesh as gbm
import os, sys, getopt
import matplotlib.pylab as plt
import numpy as np
import scipy.stats as sp
import scipy.signal as signal
import matplotlib.mlab as mlab
import math

from scipy.interpolate import UnivariateSpline

__author__ = 'cmarraj'

class SmoothSurfMesh(gbm.AbstractBaseGenerator):
  '''
  This generator creates a random smooth curve at the inlet.
  '''
  def __init__(self):
    self.__genName__ = 'smoothsurfmesh'

    # an array of parameters for current generator
    self.parameters = []

    # the format for parameters: [name, initial_value, description
    self.parameters.append(['Lx',  100.0, 'size of the fracture in X direction'])
    self.parameters.append(['Ly',    1.0, 'size of the fracture in Y direction'])
    self.parameters.append(['Lz',  500.0, 'size of the fracture in Z direction'])
    self.parameters.append(['x_res', 800, 'number of cells in X direction'])
    self.parameters.append(['y_res',   8, 'number of cells in Y direction'])
    self.parameters.append(['z_res', 256, 'number of cells in Z direction'])
    self.parameters.append(['x_G',   1.0, 'grading in X direction'])
    self.parameters.append(['y_G',   1.0, 'grading in Y direction'])
    self.parameters.append(['z_G',   1.0, 'grading in Z direction'])
    self.parameters.append(['Nx',   10000, 'number of points used to generate curve'])
    self.parameters.append(['Nseeds',100, 'number of random seeds'])
    self.parameters.append(['CVar',10, 'coefficient of variance for generated curve'])
    self.parameters.append(['K',1.0, 'seed height scalar'])
    self.parameters.append(['block_z',5.0, 'size of the first block in Z direction'])
    self.parameters.append(['RandomSeed',0, 'call a specific random geometry, 0 if none'])

  def createBlockMeshDict(self, dictFileName):
    lines = self.read_dict(dictFileName)
    empty, Lx = self.check_par('Lx', lines)
    empty, Ly = self.check_par('Ly', lines)
    empty, Lz = self.check_par('Lz', lines)
    empty, x_res = self.check_par('x_res', lines)
    empty, y_res = self.check_par('y_res', lines)
    empty, z_res = self.check_par('z_res', lines)
    empty, x_G = self.check_par('x_G', lines)
    empty, y_G = self.check_par('y_G', lines)
    empty, z_G = self.check_par('z_G', lines)
    empty, Nx = self.check_par('Nx', lines)
    empty, Nseeds = self.check_par('Nseeds', lines)
    empty, CVar = self.check_par('CVar', lines)
    empty, K = self.check_par('K', lines)
    empty, block_z = self.check_par('block_z', lines)
    empty, RandomSeed = self.check_par('RandomSeed', lines)

    # file content
    fileCont = '\n'
    fileCont += 'convertToMeters 1.0; \n'
    fileCont += '\n'

    ######################################################
    # add vertices                                       #
    ######################################################
    fileCont += 'vertices\n'
    fileCont += '(\n'
    #Lowest plane
    fileCont += '    ({0:8g} {1:8g} {2:8g})\n'.format(-Lx/2., -Ly/2., 0.0)   # 0
    fileCont += '    ({0:8g} {1:8g} {2:8g})\n'.format( Lx/2., -Ly/2., 0.0)   # 1
    fileCont += '    ({0:8g} {1:8g} {2:8g})\n'.format( Lx/2.,  Ly/2., 0.0)   # 2
    fileCont += '    ({0:8g} {1:8g} {2:8g})\n'.format(-Lx/2.,  Ly/2., 0.0)   # 3

    fileCont += '    ({0:8g} {1:8g} {2:8g})\n'.format(-Lx/2., -Ly/2., block_z)     # 4
    fileCont += '    ({0:8g} {1:8g} {2:8g})\n'.format( Lx/2., -Ly/2., block_z)     # 5
    fileCont += '    ({0:8g} {1:8g} {2:8g})\n'.format( Lx/2.,  Ly/2., block_z)     # 6
    fileCont += '    ({0:8g} {1:8g} {2:8g})\n'.format(-Lx/2.,  Ly/2., block_z)     # 7

    #Highest plane
    fileCont += '    ({0:8g} {1:8g} {2:8g})\n'.format(-Lx/2., -Ly/2., Lz)     # 8
    fileCont += '    ({0:8g} {1:8g} {2:8g})\n'.format( Lx/2., -Ly/2., Lz)     # 9
    fileCont += '    ({0:8g} {1:8g} {2:8g})\n'.format( Lx/2.,  Ly/2., Lz)     # 10
    fileCont += '    ({0:8g} {1:8g} {2:8g})\n'.format(-Lx/2.,  Ly/2., Lz)     # 11
    fileCont += ');\n'

    ######################################################
    # add edges                                          #
    ######################################################
    fileCont += '\n'
    fileCont += 'edges\n'
    fileCont += '(\n'


    fileCont += '    spline 0 1 (\n'
    while True:
      x1 = np.linspace( -Lx, 0.0, Nx)
      x = np.linspace( 0.0, Lx, Nx)
      x2 = np.linspace( Lx, 2*Lx, Nx)
      x_mesh = np.linspace( -Lx/2,Lx/2, Nx)
    
      if RandomSeed == 0:
        gen_Seed = np.random.random_integers(1,1e7)
        np.random.seed(gen_Seed)
      else:
        np.random.seed(RandomSeed)
        
      xs = Lx * np.random.random(Nseeds)
      variance = CVar * np.random.random(Nseeds)
      sigma = np.sqrt( variance )
      
      y = mlab.normpdf(x,xs[0],sigma[0])
      for i in range(1,Nseeds):
        y += mlab.normpdf(x,xs[i],sigma[i]) + \
             mlab.normpdf(x1,xs[i],sigma[i]) + \
             mlab.normpdf(x2,xs[i],sigma[i])
      
      Scale = 0.5*float(Nseeds)/y[x[0]]
      y = Scale*y/float(Nseeds)
      
      if y[0] == 0.5 and y[Nx-1] == 0.5:
          break
    
    for i in range(len(y)):
        y[i] = y[i] - 0.5
        y[i] = y[i] * K
        y[i] = y[i] + 0.5
    
    if RandomSeed == 0:
        print("Random Seed Value:") 
        print(gen_Seed)
    
    for i in range(Nx):
      fileCont += '        ( {0:f}\t{1:f}\t{2:f} )\n'.format(x_mesh[i], y[i]-1, 0.0)
    fileCont += '    )\n'

    fileCont += '    spline 3 2 (\n'
    for i in range(Nx):
      fileCont += '        ( {0:f}\t{1:f}\t{2:f} )\n'.format(x_mesh[i], y[i], 0.0)
    fileCont += '    )\n'
    

    fileCont += ');\n'

    ######################################################
    # add blocks                                         #
    ######################################################
    fileCont += '\n'
    fileCont += 'blocks\n'
    fileCont += '(\n'
    
    #####################################################################################
    # advanced grading
    # calculate the size of the lowest cell
    H = Lz - block_z
    lam = pow(z_G, 1.0/float(z_res-1))
    sum_lam = 0.0
    for i in range(z_res):
      sum_lam += pow(lam, i)

    z1 = H / sum_lam

    block_res = int(np.ceil(block_z/z1))
    
    fileCont += '    hex (0  1  2  3  4  5  6  7) ({0:d} {1:d} {2:d})\n' \
                '      simpleGrading (\n' \
                '        {3:g}\n' \
                '        {4:g}\n' \
                '        {5:g}\n' \
                '      )\n'. \
                format(x_res, y_res, block_res, x_G, y_G, 1.0)

    fileCont += '    hex (4  5  6  7 8 9 10 11) ({0:d} {1:d} {2:d})\n' \
                '      simpleGrading (\n' \
                '        {3:g}\n' \
                '        {4:g}\n' \
                '        {5:g}\n' \
                '      )\n'. \
                format(x_res, y_res, z_res, x_G, y_G, z_G)

    fileCont += ');\n'

    ######################################################
    # add boundaries                                     #
    ######################################################
    fileCont += '\n'
    fileCont += 'boundary\n'
    fileCont += '(\n'

    fileCont += '    outlet\n'
    fileCont += '    {\n'
    fileCont += '      type patch;\n'
    fileCont += '      faces\n'
    fileCont += '      (\n'
    fileCont += '        (8  9  10  11)\n'
    fileCont += '      );\n'
    fileCont += '    }\n'
    fileCont += '\n'

    fileCont += '    inlet\n'
    fileCont += '    {\n'
    fileCont += '      type patch;\n'
    fileCont += '      faces\n'
    fileCont += '      (\n'
    fileCont += '        (0   1   2   3)\n'
    fileCont += '      );\n'
    fileCont += '    }\n'
    fileCont += '\n'

    fileCont += '    walls\n'
    fileCont += '    {\n'
    fileCont += '      type wall;\n'
    fileCont += '      faces\n'
    fileCont += '      (\n'
    fileCont += '        ( 0   1   5   4)\n'
    fileCont += '        ( 4   5   9   8)\n'
    fileCont += '        ( 2   6   7   3)\n'
    fileCont += '        ( 6  10  11   7)\n'
    fileCont += '      );\n'
    fileCont += '    }\n'
    fileCont += '\n'

    fileCont += '    periodicx1\n'
    fileCont += '    {\n'
    fileCont += '      type cyclic;\n'
    fileCont += '      neighbourPatch periodicx2;\n'
    fileCont += '      faces\n'
    fileCont += '      (\n'
    fileCont += '        ( 0   3   7   4)\n'
    fileCont += '        ( 4   7  11   8)\n'
    fileCont += '      );\n'
    fileCont += '    }\n'
    fileCont += '\n'

    fileCont += '    periodicx2\n'
    fileCont += '    {\n'
    fileCont += '      type cyclic;\n'
    fileCont += '      neighbourPatch periodicx1;\n'
    fileCont += '      faces\n'
    fileCont += '      (\n'
    fileCont += '        ( 1   2   6   5)\n'
    fileCont += '        ( 5   6  10   9)\n'
    fileCont += '      );\n'
    fileCont += '    }\n'
    fileCont += '\n'

    fileCont += ');\n'

    ######################################################
    # add boundaries                                     #
    ######################################################
    fileCont += '\n'
    fileCont += 'mergePatchPairs\n'
    fileCont += '(\n'
    fileCont += ');\n\n'

    fileCont += '// ********************************************************************* //\n'

    # writing the file
    self.writeBlockMesh( fileCont )
    return True

