import os, sys
import numpy as np
import genBlockMesh as gbm

__author__ = 'sdale'

class SimpleOneBumpCyl(gbm.AbstractBaseGenerator):
  '''
  This is the generator for cylindrical meshes.
  The mesh is created with an inner rectangular prism, defined by x_rect, y_rect.
  '''

  def __init__(self):
    self.__genName__ = 'simpleOneBumpCyl'

    # an array of parameters for current generator
    self.parameters = []

    # the format for parameters: [name, initial_value, description]
    self.parameters.append(['x_rect',  1.00, 'size of the rectangular prism in X direction'])
    self.parameters.append(['y_rect',  0.25, 'size of the rectangular prism in Y direction'])
    self.parameters.append(['Lz',  500.0, 'size of the fracture in Z direction'])
    self.parameters.append(['x_res', 8, 'number of cells in X direction'])
    self.parameters.append(['y_res', 8, 'number of cells in Y direction'])
    self.parameters.append(['z_res', 2048, 'number of cells in Z direction'])
    self.parameters.append(['x_G',   1.0, 'grading in X direction'])
    self.parameters.append(['y_G',   10.0, 'grading in Y direction'])
    self.parameters.append(['z_G',   10.0, 'grading in Z direction'])
    self.parameters.append(['n_surf',   32, 'number of points on the surface of the ellipse'])
    self.parameters.append(['beta',   1.0, 'Y direction elliptical radius'])
    self.parameters.append(['alpha',   2.0, 'X direction elliptical radius'])
    self.parameters.append(['n_spline', 128, 'number of points on the bump spline'])
    self.parameters.append(['bump_z', 25.0, 'position of the bump'])
    self.parameters.append(['bump_A',  0.5, 'amplitude of the wave'])
    self.parameters.append(['bump_w',  2.0, 'wave length'])
    self.parameters.append(['bump_T',  np.pi/2, 'angle that determines the direction of the bump (in radians). Ex: 0 => bump in x, np.pi/2 => bump in y'])

  def ellipsoid_function(self, theta, delt, axis):
    alpha = self.parameters[11][1]
    beta = self.parameters[10][1]
    n_surf = self.parameters[9][1]
    #axis 0 = x, 1 = y
    if axis == 0:
    	func = alpha * np.cos(theta * np.pi/4 + delt * np.pi/2/n_surf)
    if axis == 1:
    	func =beta * np.sin(theta * np.pi/4 + delt * np.pi/2/n_surf)
    return func 
    
  def bump_function(self, delt, theta, axis):
    n_spline = self.parameters[12][1]
    bump_z = self.parameters[13][1]
    bump_A = self.parameters[14][1]
    bump_w = self.parameters[15][1]
    bump_T = self.parameters[16][1]
    #axis 0 = x, 1 = y
    if axis == 1:
    	func = np.cos(bump_T)*bump_A/2.*(1 - np.cos(2*np.pi*delt/n_spline))+self.ellipsoid_function(theta, 0, axis)
    if axis == 0:
    	func = np.sin(bump_T)*bump_A/2.*(1 - np.cos(2*np.pi*delt/n_spline))+self.ellipsoid_function(theta, 0, axis)
    return func
    
  def inside_bump_function(self, delt, sign, axis):
    n_spline = self.parameters[12][1]
    bump_z = self.parameters[13][1]
    bump_A = self.parameters[14][1]
    bump_w = self.parameters[15][1]
    x_rect = self.parameters[0][1]
    y_rect = self.parameters[1][1]
    bump_T = self.parameters[16][1]
    #axis 0 = x, 1 = y
    if axis == 0:
    	func = np.sin(bump_T)*bump_A/2.*(1-np.cos(2*np.pi*delt/n_spline))+sign*x_rect
    if axis == 1:
    	func = np.cos(bump_T)*bump_A/2.*(1-np.cos(2*np.pi*delt/n_spline))+sign*y_rect
    return func
    
  def spline_function(self, type_, point, theta, H):
    #valid inputs for type_ : ellipse, bump, inside
    #point should be a string of two numbers with at least one space separating them
    x_rect = self.parameters[0][1]
    y_rect = self.parameters[1][1]
    Lz = self.parameters[2][1]
    n_surf = self.parameters[9][1]
    n_spline = self.parameters[12][1]
    bump_z = self.parameters[13][1]
    bump_w = self.parameters[15][1]
    fileCont2 = '\n'
    n_flat = int(n_spline * (Lz - bump_w)/bump_w)
    delta_z = bump_w / float(n_spline)
    delta_z1 = (bump_z-bump_w/2) / float(n_flat)
    delta_z2 = (Lz - (bump_z+bump_w/2))/ float(n_flat)
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
    if type_ == "bump":
      fileCont2 += '    spline %s (\n' % (point)
      for i in range(0, n_flat, 1):
        currentX = self.ellipsoid_function(theta, 0, 0)
        currentY = self.ellipsoid_function(theta, 0, 1)
        currentZ = delta_z1*i
        fileCont2 += '        ( {0:f}\t{1:f}\t{2:f} )\n'.format(currentX, currentY, currentZ)
      for i in range(0, n_spline, 1):
        currentX = self.bump_function(i, theta, 0)
        currentY = self.bump_function(i, theta, 1)
        currentZ = bump_z-bump_w/2+(delta_z * i)
        fileCont2 += '        ( {0:f}\t{1:f}\t{2:f} )\n'.format(currentX, currentY, currentZ)
      for i in range(0, n_flat+1, 1):
        currentX = self.ellipsoid_function(theta, 0, 0)
        currentY = self.ellipsoid_function(theta, 0, 1)
        currentZ = bump_z+bump_w/2 + (delta_z2 * i)
        fileCont2 += '        ( {0:f}\t{1:f}\t{2:f} )\n'.format(currentX, currentY, currentZ)
      fileCont2 += '    )\n'
    if type_ == "inside":
      x_sign = 1
      y_sign = 1
      if (theta == 5 or theta == 3):
        x_sign = -1
      if (theta == 5 or theta == 7):
        y_sign = -1
      fileCont2 += '    spline %s (\n' % (point)
      for i in range(0, n_flat, 1):
        currentX =  x_sign * x_rect
        currentY =  y_sign * y_rect
        currentZ = delta_z1*i
        fileCont2 += '        ( {0:f}\t{1:f}\t{2:f} )\n'.format(currentX, currentY, currentZ)
      for i in range(0, n_spline, 1):
        currentX = self.inside_bump_function(i, x_sign, 0)
        currentY = self.inside_bump_function(i, y_sign, 1)
        currentZ = bump_z-bump_w/2+(delta_z * i)
        fileCont2 += '        ( {0:f}\t{1:f}\t{2:f} )\n'.format(currentX, currentY, currentZ)
      for i in range(0, n_flat+1, 1):
        currentX = x_sign * x_rect
        currentY = y_sign * y_rect
        currentZ = bump_z+bump_w/2 + (delta_z2 * i)
        fileCont2 += '        ( {0:f}\t{1:f}\t{2:f} )\n'.format(currentX, currentY, currentZ)
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
    empty, n_spline = self.check_par('n_spline', lines)
    empty, bump_z = self.check_par('bump_z', lines)
    empty, bump_A = self.check_par('bump_A', lines)
    empty, bump_w = self.check_par('bump_w', lines)
    empty, bump_T = self.check_par('bump_T', lines)
            
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
    
    #Highest plane
    fileCont += '    ({0:8g} {1:8g} {2:8g})\n'.format(self.ellipsoid_function(5, 0, 0), self.ellipsoid_function(5, 0, 1), Lz)   # 8
    fileCont += '    ({0:8g} {1:8g} {2:8g})\n'.format(self.ellipsoid_function(7, 0, 0), self.ellipsoid_function(7, 0, 1), Lz)   # 9
    fileCont += '    ({0:8g} {1:8g} {2:8g})\n'.format( 								 -x_rect,  								  -y_rect, Lz)   #10
    fileCont += '    ({0:8g} {1:8g} {2:8g})\n'.format( 								  x_rect,  							 	  -y_rect, Lz)   #11
    fileCont += '    ({0:8g} {1:8g} {2:8g})\n'.format( 								 -x_rect,  								   y_rect, Lz)   #12
    fileCont += '    ({0:8g} {1:8g} {2:8g})\n'.format( 								  x_rect,  								   y_rect, Lz)   #13
    fileCont += '    ({0:8g} {1:8g} {2:8g})\n'.format(self.ellipsoid_function(3, 0, 0), self.ellipsoid_function(3, 0, 1), Lz)   #14
    fileCont += '    ({0:8g} {1:8g} {2:8g})\n'.format(self.ellipsoid_function(1, 0, 0), self.ellipsoid_function(1, 0, 1), Lz)   #15

    fileCont += ');\n'

    ######################################################
    # add edges                                          #
    ######################################################
    fileCont += '\n'
    fileCont += 'edges\n'
    fileCont += '(\n'

    #******Cylinder******
    #Inlet
    fileCont += self.spline_function("ellipse", "0 1", 5, 0)
    fileCont += self.spline_function("ellipse", "6 7", 1, 0)
    fileCont += self.spline_function("ellipse", "0 6", 3, 0)
    fileCont += self.spline_function("ellipse", "1 7", 7, 0)
    #Outlet    
    fileCont += self.spline_function("ellipse", " 8  9", 5, Lz)
    fileCont += self.spline_function("ellipse", "14 15", 1, Lz)
    fileCont += self.spline_function("ellipse", " 8 14", 3, Lz)
    fileCont += self.spline_function("ellipse", " 9 15", 7, Lz)
    #********************

    #********Bump********
    #Outside
    fileCont += self.spline_function("bump", " 0  8", 5, 0)
    fileCont += self.spline_function("bump", " 7 15", 1, 0)
    fileCont += self.spline_function("bump", " 6 14", 3, 0)
    fileCont += self.spline_function("bump", " 1  9", 7, 0)
    #Inside
    fileCont += self.spline_function("inside", " 2 10", 5, 0)
    fileCont += self.spline_function("inside", " 5 13", 1, 0)
    fileCont += self.spline_function("inside", " 4 12", 3, 0)
    fileCont += self.spline_function("inside", " 3 11", 7, 0)
    #********************

    fileCont += ');\n'

    ######################################################
    # add blocks                                         #
    ######################################################
    fileCont += '\n'
    fileCont += 'blocks\n'
    fileCont += '(\n'

    #####################################################################################
    #block 0
    fileCont += '    hex (0  1  3  2   8   9  11  10) ({0:d} {1:d} {2:d})\n' \
                '      simpleGrading (\n' \
                '        {3:g}\n' \
                '        {4:g}\n' \
                '        {5:g}\n' \
                '      )\n'. \
                format(x_res, y_res, z_res, x_G, y_G, z_G)
    #block 1            
    fileCont += '    hex (7  6  4  5  15  14  12  13) ({0:d} {1:d} {2:d})\n' \
                '      simpleGrading (\n' \
                '        {3:g}\n' \
                '        {4:g}\n' \
                '        {5:g}\n' \
                '      )\n'. \
                format(x_res, y_res, z_res, x_G, y_G, z_G)

    #block 2            
    fileCont += '    hex (6  0  2  4  14   8  10  12) ({0:d} {1:d} {2:d})\n' \
                '      simpleGrading (\n' \
                '        {3:g}\n' \
                '        {4:g}\n' \
                '        {5:g}\n' \
                '      )\n'. \
                format(x_res, y_res, z_res, x_G, y_G, z_G)
                
    #block 3            
    fileCont += '    hex (1  7  5  3   9  15  13  11) ({0:d} {1:d} {2:d})\n' \
                '      simpleGrading (\n' \
                '        {3:g}\n' \
                '        {4:g}\n' \
                '        {5:g}\n' \
                '      )\n'. \
                format(x_res, y_res, z_res, x_G, y_G, z_G)
                
    #block 3            
    fileCont += '    hex (2  3  5  4  10  11  13  12) ({0:d} {1:d} {2:d})\n' \
                '      simpleGrading (\n' \
                '        {3:g}\n' \
                '        {4:g}\n' \
                '        {5:g}\n' \
                '      )\n'. \
                format(x_res, x_res, z_res, 1, 1, z_G)
                
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
    fileCont += '        ( 8   9  11  10)\n'
    fileCont += '        (14  12  13  15)\n'
    fileCont += '        ( 8  10  12  14)\n'
    fileCont += '        ( 9  15  13  11)\n'
    fileCont += '        (10  11  13  12)\n'
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
