import os, sys
import numpy as np
import genBlockMesh as gbm

__author__ = 'sdale'

class Curved(gbm.AbstractBaseGenerator):
  '''
  This is the generator for curved meshes.
  The mesh may be generated with symmetric or asymmetric curved boundaries.
  		f = -1 ********* f = 1
  
  		|__|  ********  |__|
  		(__(  ******** (____)
 		)__)  ********  )__(     ^
 		(__(  ******** (____)    |
 		|__|  ********  |__|     z --- y->
  '''

  def __init__(self):
    self.__genName__ = 'curved'

    # an array of parameters for current generator
    self.parameters = []

    # the format for parameters: [name, initial_value, description]
    self.parameters.append(['Lz',  500.0, 'size of the fracture in Z direction'])
    self.parameters.append(['Ly',   1.00, 'width of the fracture (actual width = 1/2 of value)'])
    self.parameters.append(['Lx',   2.00, 'depth of the fracture (fracture begins at -x and ends at +x)'])
    self.parameters.append(['x_res', 4, 'number of cells in X direction'])
    self.parameters.append(['y_res', 8, 'number of cells in Y direction'])
    self.parameters.append(['z_res', 4096, 'number of cells in Z direction'])
    self.parameters.append(['x_G',   1.0, 'grading in X direction'])
    self.parameters.append(['y_G',   1.0, 'grading in Y direction'])
    self.parameters.append(['z_G',   1.0, 'grading in Z direction'])
    self.parameters.append(['n_',   5000, 'number of points on the curved boundaries'])
    self.parameters.append(['amp',   0.25, 'wave amplitude factor'])
    self.parameters.append(['wave_',   50, 'number of waves on the curved boundary'])
    self.parameters.append(['symm',   -1, 'symmetry factor, f= -1 for asymmetric waves, f= 1 for symmetric waves'])
    
  def curved_function(self, f, delta):
    amp = self.parameters[10][1]
    Lx = self.parameters[2][1]
    wave_ = self.parameters[11][1]
    n_ = self.parameters[9][1]
    return amp*(1 - f * Lx*np.cos(2*wave_*np.pi*delta/float(n_)))
    
  def createBlockMeshDict(self, dictFileName):
  	 # read parameters
    lines = self.read_dict(dictFileName)
    empty, Lz = self.check_par('Lz', lines)
    empty, Ly = self.check_par('Ly', lines)
    empty, Lx = self.check_par('Lx', lines)
    empty, x_res = self.check_par('x_res', lines)
    empty, y_res = self.check_par('y_res', lines)
    empty, z_res = self.check_par('z_res', lines)
    empty, x_G = self.check_par('x_G', lines)
    empty, y_G = self.check_par('y_G', lines)
    empty, z_G = self.check_par('z_G', lines)
    empty, n_ = self.check_par('n_', lines)
    empty, amp = self.check_par('amp', lines)
    empty, wave_ = self.check_par('wave_', lines)
    empty, symm = self.check_par('symm', lines)

    fileCont = '\n'
    fileCont += 'convertToMeters 1.0; \n'
    fileCont += '\n'

    ######################################################
    # add vertices                                       #
    ######################################################
    fileCont += 'vertices\n'
    fileCont += '(\n'
    
    ome = np.pi * 2 * wave_ / Lz
    dphi = Lz / float(n_)
    phi = np.pi / 2
    
    #Lowest plane
    fileCont += '    ({0:8g} {1:8g} {2:8g})\n'.format(-Lx, -self.curved_function(symm, 0), 0.0)   # 0
    fileCont += '    ({0:8g} {1:8g} {2:8g})\n'.format( Lx, -self.curved_function(symm, 0), 0.0)   # 1
    fileCont += '    ({0:8g} {1:8g} {2:8g})\n'.format( Lx,  self.curved_function(1, 0),    0.0)   # 2
    fileCont += '    ({0:8g} {1:8g} {2:8g})\n'.format(-Lx,  self.curved_function(1, 0),    0.0)   # 3
   
    #Highest plane
    fileCont += '    ({0:8g} {1:8g} {2:8g})\n'.format(-Lx, -self.curved_function(symm, 0),  Lz)   # 4
    fileCont += '    ({0:8g} {1:8g} {2:8g})\n'.format( Lx, -self.curved_function(symm, 0),  Lz)   # 5
    fileCont += '    ({0:8g} {1:8g} {2:8g})\n'.format( Lx,  self.curved_function(1, 0),     Lz)   # 6
    fileCont += '    ({0:8g} {1:8g} {2:8g})\n'.format(-Lx,  self.curved_function(1, 0),     Lz)   # 7

    fileCont += ');\n'
    
    ######################################################
    # add edges                                          #
    ######################################################
    fileCont += '\n'
    fileCont += 'edges\n'
    fileCont += '(\n'


    fileCont += '    spline 0 4 (\n'
    for i in range(0, n_+1, 1):
      currentX = -Lx
      currentY = -self.curved_function(symm, i)
      currentZ =  dphi * i
      fileCont += '        ( {0:f}\t{1:f}\t{2:f} )\n'.format(currentX, currentY, currentZ)
    fileCont += '    )\n'

    fileCont += '    spline 1 5 (\n'
    for i in range(0, n_+1, 1):
      currentX =  Lx
      currentY = -self.curved_function(symm, i)
      currentZ =  dphi * i
      fileCont += '        ( {0:f}\t{1:f}\t{2:f} )\n'.format(currentX, currentY, currentZ)
    fileCont += '    )\n'

    fileCont += '    spline 2 6 (\n'
    for i in range(0, n_+1, 1):
      currentX =  Lx
      currentY =  self.curved_function(1, i)
      currentZ =  dphi * i
      fileCont += '        ( {0:f}\t{1:f}\t{2:f} )\n'.format(currentX, currentY, currentZ)
    fileCont += '    )\n'

    fileCont += '    spline 3 7 (\n'
    for i in range(0, n_+1, 1):
      currentX = -Lx
      currentY =  self.curved_function(1, i)
      currentZ =  dphi * i
      fileCont += '        ( {0:f}\t{1:f}\t{2:f} )\n'.format(currentX, currentY, currentZ)
    fileCont += '    )\n'
  
  
    fileCont += ');\n'


    ######################################################
    # add blocks                                         #
    ######################################################
    fileCont += '\n'
    fileCont += 'blocks\n'
    fileCont += '(\n'

    #####################################################################################
    #block 0
    fileCont += '    hex (0  1  2  3  4  5  6  7) ({0:d} {1:d} {2:d})\n' \
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
    fileCont += '        ( 4  5  6  7)\n'
    fileCont += '      );\n'
    fileCont += '    }\n'
    fileCont += '\n'

    fileCont += '    inlet\n'
    fileCont += '    {\n'
    fileCont += '      type patch;\n'
    fileCont += '      faces\n'
    fileCont += '      (\n'
    fileCont += '        ( 0   3   2   1)\n'
    fileCont += '      );\n'
    fileCont += '    }\n'
    fileCont += '\n'

    fileCont += '    walls\n'
    fileCont += '    {\n'
    fileCont += '      type wall;\n'
    fileCont += '      faces\n'
    fileCont += '      (\n'
    fileCont += '        ( 0   1   5   4)\n'
    
    fileCont += '        ( 2   3   7   6)\n'
    fileCont += '      );\n'
    fileCont += '    }\n'
    fileCont += '\n'
  
    fileCont += '    periodicx1\n'
    fileCont += '    {\n'
    fileCont += '      type cyclic;\n'
    fileCont += '      neighbourPatch periodicx2;\n'
    fileCont += '      faces\n'
    fileCont += '      (\n'
    fileCont += '        ( 4   7   3   0)\n'
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
    fileCont += '      );\n'
    fileCont += '    }\n'
    fileCont += '\n'

    fileCont += ');\n'
  # *****************************************************************************
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








    
