import os, sys
import numpy as np
import genBlockMesh as gbm

__author__ = 'vstar'

class Seeded2blocks(gbm.AbstractBaseGenerator):
  '''
  This generator creates a seeded fracture. Cosine like seed is at the inlet at Lx/2.
  blockMeshDict consist of two blocks. First block contains the seed and the second
  is a flat fracture.
  '''
  def __init__(self):
    self.__genName__ = 'seeded2blocks'

    # an array of parameters for current generator
    self.parameters = []

    # the format for parameters: [name, initial_value, description
    self.parameters.append(['Lx',  100.0, 'size of the fracture in X direction'])
    self.parameters.append(['Ly',    1.0, 'size of the fracture in Y direction'])
    self.parameters.append(['Lz',  500.0, 'size of the fracture in Z direction'])
    self.parameters.append(['x_res', 200, 'number of cells in X direction'])
    self.parameters.append(['y_res',   8, 'number of cells in Y direction'])
    self.parameters.append(['z_res', 512, 'number of cells in Z direction'])
    self.parameters.append(['x_G',   5.0, 'grading in X direction'])
    self.parameters.append(['y_G',   2.0, 'grading in Y direction'])
    self.parameters.append(['z_G',   1.0, 'grading in Z direction'])
    self.parameters.append(['nseed',    400, 'number of points for the spline'])
    self.parameters.append(['seed_x',   5.0, 'size of the seed in X direction'])
    self.parameters.append(['seed_y',   0.5, 'size of the seed in Y direction'])
    self.parameters.append(['seed_z',   5.0, 'size of the seed in Z direction'])

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
    empty, n = self.check_par('nseed', lines)
    empty, seed_x = self.check_par('seed_x', lines)
    empty, seed_y = self.check_par('seed_y', lines)
    empty, seed_z = self.check_par('seed_z', lines)

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

    fileCont += '    ({0:8g} {1:8g} {2:8g})\n'.format(-Lx/2., -Ly/2., seed_z)     # 4
    fileCont += '    ({0:8g} {1:8g} {2:8g})\n'.format( Lx/2., -Ly/2., seed_z)     # 5
    fileCont += '    ({0:8g} {1:8g} {2:8g})\n'.format( Lx/2.,  Ly/2., seed_z)     # 6
    fileCont += '    ({0:8g} {1:8g} {2:8g})\n'.format(-Lx/2.,  Ly/2., seed_z)     # 7

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

    xx = np.linspace(-Lx/2.,Lx/2.,n+1)
    #yi = Ly/2. + seed_y*(1+np.cos(2*np.pi*x/W))

    fileCont += '    spline 0 1 (\n'
    for xi in xx:
      yi = -Ly/2.
      if xi>-seed_x/2. and xi<seed_x/2.:
        yi -= seed_y/2.*(1+np.cos(2.*np.pi*xi/seed_x))
      fileCont += '        ( {0:f}\t{1:f}\t{2:f} )\n'.format(xi, yi, 0.0)
    fileCont += '    )\n'

    fileCont += '    spline 3 2 (\n'
    for xi in xx:
      yi = Ly/2.
      if xi>-seed_x/2. and xi<seed_x/2.:
        yi += seed_y/2.*(1+np.cos(2.*np.pi*xi/seed_x))
      fileCont += '        ( {0:f}\t{1:f}\t{2:f} )\n'.format(xi, yi, 0.0)
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
    H = Lz - seed_z
    lam = pow(z_G, 1.0/float(z_res-1))
    sum_lam = 0.0
    for i in range(z_res):
      sum_lam += pow(lam, i)

    z1 = H / sum_lam

    seed_res = int(np.ceil(seed_z/z1))

    # by default symmetric grading across the fracture
    y_Gs = '((0.5 0.5 {0:g}) (0.5 0.5 {1:g}))'.format(y_G, 1.0/y_G)

    # Grading in X direction. We want the center be more detaled
    Lseed = seed_x * 3.2
    Gside = x_G
    Lside = (Lx - Lseed)/2.

    cur_side_res=int((x_res-2)/2.)
    cur_seed_res=x_res-2*cur_side_res

    while True:
      lam = pow(Gside, 1.0/float(cur_side_res-1))
      sum_lam = 0.0
      for i in range(cur_side_res):
        sum_lam += pow(lam, i)

      x1 = Lside / sum_lam
      x2 = Lseed / float(cur_seed_res)

      if(x2<x1):
        cur_side_res-=1
        cur_seed_res+=2
        break
      else:
        cur_side_res-=1
        cur_seed_res+=2

    if cur_seed_res>x_res:
      print('cur_seed_res={0:d} cur_side_res={1:d}\n'.format(cur_seed_res,cur_side_res))
      print('*** Error. Change your X grading.\n')
      return False

    side_cell_part = cur_side_res/float(x_res)
    seed_cell_part = cur_seed_res/float(x_res)

    x_Gs = '(({0:g} {1:g} {2:g}) ({3:g} {4:g} {5:g}) ({6:g} {7:g} {8:g}))'.\
      format(Lside/Lx, side_cell_part, 1/Gside,
             Lseed/Lx, seed_cell_part, 1.0,
             Lside/Lx, side_cell_part, Gside)
    #####################################################################################


    fileCont += '    hex (0  1  2  3  4  5  6  7) ({0:d} {1:d} {2:d})\n' \
                '      simpleGrading (\n' \
                '        {3:s}\n' \
                '        {4:s}\n' \
                '        {5:g}\n' \
                '      )\n'. \
                format(x_res, y_res, seed_res, x_Gs, y_Gs, 1.0)

    fileCont += '    hex (4  5  6  7 8 9 10 11) ({0:d} {1:d} {2:d})\n' \
                '      simpleGrading (\n' \
                '        {3:s}\n' \
                '        {4:s}\n' \
                '        {5:g}\n' \
                '      )\n'. \
                format(x_res, y_res, z_res, x_Gs, y_Gs, z_G)

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

