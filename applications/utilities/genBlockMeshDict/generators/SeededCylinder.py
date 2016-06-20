import os, sys
import numpy as np
import genBlockMesh as gbm

__author__ = 'sdale'

class SeededCylinder(gbm.AbstractBaseGenerator):
  '''
  This is the generator for seeded cylindrical meshes.
  The mesh is created with an inner rectangular prism, defined by x_rect, y_rect.
  '''

  def __init__(self):
    self.__genName__ = 'seededCylinder'

    # an array of parameters for current generator
    self.parameters = []

    # the format for parameters: [name, initial_value, description]
    self.parameters.append(['x_rect',  0.75, 'size of the rectangular prism in x direction'])
    self.parameters.append(['y_rect',  0.25, 'size of the rectangular prism in y direction'])
    self.parameters.append(['Lz',  500.0, 'size of the fracture in Z direction'])
    self.parameters.append(['x_res', 16, 'number of cells in X direction'])
    self.parameters.append(['y_res', 8, 'number of cells in Y direction'])
    self.parameters.append(['z_res', 512, 'number of cells in Z direction'])
    self.parameters.append(['x_G',   1.0, 'grading in X direction'])
    self.parameters.append(['y_G',   10.0, 'grading in Y direction'])
    self.parameters.append(['z_G',   10.0, 'grading in Z direction'])
    self.parameters.append(['n_surf', 128, 'number of points on the edges of the cylinder'])
    self.parameters.append(['beta',   1.0, 'Y direction elliptical radius'])
    self.parameters.append(['alpha',   1.5, 'X direction elliptical radius'])
    self.parameters.append(['nseed',    400, 'number of points for the spline'])
    self.parameters.append(['seed_x',   0.5, 'size of the seed in X direction'])
    self.parameters.append(['seed_y',   0.5, 'size of the seed in Y direction'])
    self.parameters.append(['seed_z',   5.0, 'size of the seed in Z direction'])
  def ellipsoid_function(self, theta, delt, axis):
    alpha = self.parameters[11][1]
    beta = self.parameters[10][1]
    n_surf = self.parameters[9][1]
    #axis 0 = x, 1 = y
    if axis == 0:
    	func = alpha * np.cos(theta * np.pi/4 + delt * np.pi/2/n_surf)
    if axis == 1:
    	func = beta * np.sin(theta * np.pi/4 + delt * np.pi/2/n_surf)
    return func
  def spline_function(self, type_, point, theta, H):
    #valid inputs for type_ : ellipse, seed
    #point should be a string of two numbers with at least one space separating them
    x_rect = self.parameters[0][1]
    y_rect = self.parameters[1][1]
    Lz = self.parameters[2][1]
    n_surf = self.parameters[9][1]
    beta = self.parameters[10][1]
    alpha = self.parameters[11][1]
    nseed = self.parameters[12][1]
    seed_x = self.parameters[13][1]
    seed_y = self.parameters[14][1]
    seed_z = self.parameters[15][1]
    fileCont2 = '\n'
    if type_ == "ellipse":
      if (theta == 5 or theta == 7):
        range1 = 0
        range2 = n_surf+1
        step = 1
      if (theta == 1 or theta == 3):
        range1 = n_surf
        range2 = -1
        step = -1
      fileCont2 += '    spline %s (\n' % (point)
      for i in range(range1, range2, step):
        currentX = self.ellipsoid_function(theta, i, 0)
        currentY = self.ellipsoid_function(theta, i, 1)
        currentZ = H
        fileCont2 += '        ( {0:f}\t{1:f}\t{2:f} )\n'.format(currentX, currentY, currentZ)
      fileCont2 += '    )\n'
    if type_ == "seed":
      seed_start = -seed_x/2
      seed_end = seed_x/2
      coef = 1
      if (theta == 5 or theta == 7):
        coef = -1
      fileCont2 += '    spline %s (\n' % (point)
      for i in range (0, n_surf+1, 1):
        currentX = self.ellipsoid_function(theta, i, 0)
        currentY = self.ellipsoid_function(theta, i, 1)
        if (currentX>seed_start and currentX<seed_end):
            currentY = currentY + coef*seed_y/2.*(1+np.cos(2.*np.pi*currentX/seed_x))
        fileCont2 += '        ( {0:f}\t{1:f}\t{2:f} )\n'.format(currentX, currentY, H)
      fileCont2 += '    )\n'
      
    return fileCont2

  def createBlockMeshDict(self, dictFileName):
  	 # read parameters
    lines = self.read_dict(dictFileName)
    empty, x_rect = self.check_par('x_rect', lines)
    empty, y_rect = self.check_par('y_rect', lines)
    empty, Lz = self.check_par('Lz', lines)
    empty, x_res = self.check_par('x_res', lines)
    empty, y_res = self.check_par('y_res', lines)
    empty, z_res = self.check_par('z_res', lines)
    empty, x_G = self.check_par('x_G', lines)
    empty, y_G = self.check_par('y_G', lines)
    empty, z_G = self.check_par('z_G', lines)
    empty, n_surf = self.check_par('n_surf', lines)
    empty, beta = self.check_par('beta', lines)
    empty, alpha = self.check_par('alpha', lines)
    empty, nseed = self.check_par('nseed', lines)
    empty, seed_x = self.check_par('seed_x', lines)
    empty, seed_y = self.check_par('seed_y', lines)
    empty, seed_z = self.check_par('seed_z', lines) 
    
    fileCont = '\n'
    fileCont += 'convertToMeters 1.0; \n'
    fileCont += '\n'

    ######################################################
    # add vertices                                       #
    ######################################################
    fileCont += 'vertices\n'
    fileCont += '(\n'

    #Lowest plane
    fileCont += '    ({0:8g} {1:8g} {2:8g})\n'.format(self.ellipsoid_function(5, 0, 0), self.ellipsoid_function(5, 0, 1), 0.0)   # 0
    fileCont += '    ({0:8g} {1:8g} {2:8g})\n'.format(self.ellipsoid_function(7, 0, 0), self.ellipsoid_function(7, 0, 1), 0.0)   # 1
    fileCont += '    ({0:8g} {1:8g} {2:8g})\n'.format( 								 -x_rect,  								  -y_rect, 0.0)   # 2
    fileCont += '    ({0:8g} {1:8g} {2:8g})\n'.format( 								  x_rect,  								  -y_rect, 0.0)   # 3
    fileCont += '    ({0:8g} {1:8g} {2:8g})\n'.format( 								 -x_rect,  								   y_rect, 0.0)   # 4
    fileCont += '    ({0:8g} {1:8g} {2:8g})\n'.format( 								  x_rect,  								   y_rect, 0.0)   # 5
    fileCont += '    ({0:8g} {1:8g} {2:8g})\n'.format(self.ellipsoid_function(3, 0, 0), self.ellipsoid_function(3, 0, 1), 0.0)   # 6
    fileCont += '    ({0:8g} {1:8g} {2:8g})\n'.format(self.ellipsoid_function(1, 0, 0), self.ellipsoid_function(1, 0, 1), 0.0)   # 7
    
    #End of Seed plane
    fileCont += '    ({0:8g} {1:8g} {2:8g})\n'.format(self.ellipsoid_function(5, 0, 0), self.ellipsoid_function(5, 0, 1),  seed_z)   # 8
    fileCont += '    ({0:8g} {1:8g} {2:8g})\n'.format(self.ellipsoid_function(7, 0, 0), self.ellipsoid_function(7, 0, 1),  seed_z)   # 9
    fileCont += '    ({0:8g} {1:8g} {2:8g})\n'.format( 								 -x_rect,  								  -y_rect,  seed_z)   #10
    fileCont += '    ({0:8g} {1:8g} {2:8g})\n'.format( 								  x_rect,  							 	  -y_rect,  seed_z)   #11
    fileCont += '    ({0:8g} {1:8g} {2:8g})\n'.format( 								 -x_rect,  								   y_rect,  seed_z)   #12
    fileCont += '    ({0:8g} {1:8g} {2:8g})\n'.format( 								  x_rect,  								   y_rect,  seed_z)   #13
    fileCont += '    ({0:8g} {1:8g} {2:8g})\n'.format(self.ellipsoid_function(3, 0, 0), self.ellipsoid_function(3, 0, 1),  seed_z)   #14
    fileCont += '    ({0:8g} {1:8g} {2:8g})\n'.format(self.ellipsoid_function(1, 0, 0), self.ellipsoid_function(1, 0, 1),  seed_z)   #15

    #Highest plane
    fileCont += '    ({0:8g} {1:8g} {2:8g})\n'.format(self.ellipsoid_function(5, 0, 0), self.ellipsoid_function(5, 0, 1),  Lz)   #16
    fileCont += '    ({0:8g} {1:8g} {2:8g})\n'.format(self.ellipsoid_function(7, 0, 0), self.ellipsoid_function(7, 0, 1),  Lz)   #17
    fileCont += '    ({0:8g} {1:8g} {2:8g})\n'.format( 								 -x_rect,  								  -y_rect,  Lz)   #18
    fileCont += '    ({0:8g} {1:8g} {2:8g})\n'.format( 								  x_rect,  							 	  -y_rect,  Lz)   #19
    fileCont += '    ({0:8g} {1:8g} {2:8g})\n'.format( 								 -x_rect,  								   y_rect,  Lz)   #20
    fileCont += '    ({0:8g} {1:8g} {2:8g})\n'.format( 								  x_rect,  								   y_rect,  Lz)   #21
    fileCont += '    ({0:8g} {1:8g} {2:8g})\n'.format(self.ellipsoid_function(3, 0, 0), self.ellipsoid_function(3, 0, 1),  Lz)   #22
    fileCont += '    ({0:8g} {1:8g} {2:8g})\n'.format(self.ellipsoid_function(1, 0, 0), self.ellipsoid_function(1, 0, 1),  Lz)   #23

    fileCont += ');\n'

    ######################################################
    # add edges                                          #
    ######################################################
    fileCont += '\n'
    fileCont += 'edges\n'
    fileCont += '(\n'

    #******Cylinder******
    #Inlet
    fileCont += self.spline_function("seed", "0 1", 5, 0)
    fileCont += self.spline_function("seed", "7 6", 1, 0)
    fileCont += self.spline_function("ellipse", "0 6", 3, 0)
    fileCont += self.spline_function("ellipse", "1 7", 7, 0)
    #Middle
    fileCont += self.spline_function("ellipse", " 8  9", 5, seed_z)
    fileCont += self.spline_function("ellipse", "14 15", 1, seed_z)
    fileCont += self.spline_function("ellipse", " 8 14", 3, seed_z)
    fileCont += self.spline_function("ellipse", " 9 15", 7, seed_z)
    #Outlet
    fileCont += self.spline_function("ellipse", "16 17", 5, Lz)
    fileCont += self.spline_function("ellipse", "22 23", 1, Lz)
    fileCont += self.spline_function("ellipse", "16 22", 3, Lz)
    fileCont += self.spline_function("ellipse", "17 23", 7, Lz)
    #********************

    fileCont += ');\n'

    ######################################################
    # add blocks                                         #
    ######################################################
    fileCont += '\n'
    fileCont += 'blocks\n'
    fileCont += '(\n'

    #####################################################################################
    z_G_seed = 1
    H = Lz - seed_z
    lam = pow(z_G, 1.0/float(z_res-1))
    sum_lam = 0.0
    for i in range(z_res):
      sum_lam += pow(lam, i)

    z1 = H / sum_lam

    seed_res = int(np.ceil(seed_z/z1))
    x_res1 = int(2 * x_res * alpha / beta)
    #block 0
    fileCont += '    hex (0  1  3  2   8   9  11  10) ({0:d} {1:d} {2:d})\n' \
                '      simpleGrading (\n' \
                '        {3:g}\n' \
                '        {4:g}\n' \
                '        {5:g}\n' \
                '      )\n'. \
                format(x_res1, y_res, seed_res, x_G, y_G, z_G_seed)
    #block 1            
    fileCont += '    hex (7  6  4  5  15  14  12  13) ({0:d} {1:d} {2:d})\n' \
                '      simpleGrading (\n' \
                '        {3:g}\n' \
                '        {4:g}\n' \
                '        {5:g}\n' \
                '      )\n'. \
                format(x_res1, y_res, seed_res, x_G, y_G, z_G_seed)
    #block 2            
    fileCont += '    hex (6  0  2  4  14   8  10  12) ({0:d} {1:d} {2:d})\n' \
                '      simpleGrading (\n' \
                '        {3:g}\n' \
                '        {4:g}\n' \
                '        {5:g}\n' \
                '      )\n'. \
                format(x_res, y_res, seed_res, x_G, y_G, z_G_seed)
    #block 3            
    fileCont += '    hex (1  7  5  3   9  15  13  11) ({0:d} {1:d} {2:d})\n' \
                '      simpleGrading (\n' \
                '        {3:g}\n' \
                '        {4:g}\n' \
                '        {5:g}\n' \
                '      )\n'. \
                format(x_res, y_res, seed_res, x_G, y_G, z_G_seed)
    #block 4            
    fileCont += '    hex (2  3  5  4  10  11  13  12) ({0:d} {1:d} {2:d})\n' \
                '      simpleGrading (\n' \
                '        {3:g}\n' \
                '        {4:g}\n' \
                '        {5:g}\n' \
                '      )\n'. \
                format(x_res1, x_res, seed_res, 1, 1, z_G_seed)
    #block 5
    fileCont += '    hex (8  9 11 10  16  17  19  18) ({0:d} {1:d} {2:d})\n' \
                '      simpleGrading (\n' \
                '        {3:g}\n' \
                '        {4:g}\n' \
                '        {5:g}\n' \
                '      )\n'. \
                format(x_res1, y_res, z_res, x_G, y_G, z_G)
    #block 6            
    fileCont += '    hex (15 14 12 13  23  22  20  21) ({0:d} {1:d} {2:d})\n' \
                '      simpleGrading (\n' \
                '        {3:g}\n' \
                '        {4:g}\n' \
                '        {5:g}\n' \
                '      )\n'. \
                format(x_res1, y_res, z_res, x_G, y_G, z_G)

    #block 7            
    fileCont += '    hex (14  8 10 12  22  16  18  20) ({0:d} {1:d} {2:d})\n' \
                '      simpleGrading (\n' \
                '        {3:g}\n' \
                '        {4:g}\n' \
                '        {5:g}\n' \
                '      )\n'. \
                format(x_res, y_res, z_res, x_G, y_G, z_G)
                
    #block 8            
    fileCont += '    hex (9 15 13 11  17  23  21  19) ({0:d} {1:d} {2:d})\n' \
                '      simpleGrading (\n' \
                '        {3:g}\n' \
                '        {4:g}\n' \
                '        {5:g}\n' \
                '      )\n'. \
                format(x_res, y_res, z_res, x_G, y_G, z_G)
                
    #block 9            
    fileCont += '    hex (10 11 13 12  18  19  21  20) ({0:d} {1:d} {2:d})\n' \
                '      simpleGrading (\n' \
                '        {3:g}\n' \
                '        {4:g}\n' \
                '        {5:g}\n' \
                '      )\n'. \
                format(x_res1, x_res, z_res, 1, 1, z_G)
                
    fileCont += ');\n'

    ######################################################
    # add boundaries                                     #TODO
    ######################################################
    fileCont += '\n'
    fileCont += 'boundary\n'
    fileCont += '(\n'

    fileCont += '    outlet\n'
    fileCont += '    {\n'
    fileCont += '      type patch;\n'
    fileCont += '      faces\n'
    fileCont += '      (\n'
    fileCont += '        (16  17  19  18)\n'
    fileCont += '        (22  20  21  23)\n'
    fileCont += '        (16  18  20  22)\n'
    fileCont += '        (17  23  21  19)\n'
    fileCont += '        (18  19  21  20)\n'
    fileCont += '      );\n'
    fileCont += '    }\n'
    fileCont += '\n'

    fileCont += '    inlet\n'
    fileCont += '    {\n'
    fileCont += '      type patch;\n'
    fileCont += '      faces\n'
    fileCont += '      (\n'
    fileCont += '        ( 0   1   3   2)\n'
    fileCont += '        ( 6   4   5   7)\n'
    fileCont += '        ( 0   2   4   6)\n'
    fileCont += '        ( 1   7   5   3)\n'
    fileCont += '        ( 2   3   5   4)\n'
    fileCont += '      );\n'
    fileCont += '    }\n'
    fileCont += '\n'

    fileCont += '    walls\n'
    fileCont += '    {\n'
    fileCont += '      type wall;\n'
    fileCont += '      faces\n'
    fileCont += '      (\n'
    fileCont += '        ( 0   1   9   8)\n'
    fileCont += '        ( 6   7  15  14)\n'
    fileCont += '        ( 6   0   8  14)\n'
    fileCont += '        ( 1   7  15   9)\n'
    fileCont += '        ( 8   9  17  16)\n'
    fileCont += '        (14  15  23  22)\n'
    fileCont += '        (14   8  16  22)\n'
    fileCont += '        ( 9  15  23  17)\n'
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
