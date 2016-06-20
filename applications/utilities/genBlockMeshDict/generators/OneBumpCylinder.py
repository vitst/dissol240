import os, sys
import numpy as np
import genBlockMesh as gbm

__author__ = 'sdale'

class OneBumpCylinder(gbm.AbstractBaseGenerator):
  '''
  This is the generator for cylindrical meshes.
  The mesh is created with an inner rectangular prism, defined by x_rect, y_rect.
  '''

  def __init__(self):
    self.__genName__ = 'oneBumpCylinder'

    # an array of parameters for current generator
    self.parameters = []

    # the format for parameters: [name, initial_value, description]
    self.parameters.append(['x_rect',  1.00, 'size of the rectangular prism in X direction'])
    self.parameters.append(['y_rect',  0.25, 'size of the rectangular prism in Y direction'])
    self.parameters.append(['Lz',  500.0, 'size of the fracture in Z direction'])
    self.parameters.append(['x_res', 8, 'number of cells in X direction'])
    self.parameters.append(['y_res', 8, 'number of cells in Y direction'])
    self.parameters.append(['z_res', 1024, 'Baseline for determining number of cells in Z direction (actual number will be larger).'])
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
    n_flat = (int)((n_spline * (Lz - bump_w)/bump_w)/2)
    Lz_bump = (int)(7 * bump_w)
    Lz_above = (int)(bump_z + Lz_bump/2)
    Lz_below = (int)(bump_z - Lz_bump/2)
    delta_z = bump_w / float(n_spline)
    delta_z1 = (bump_z-bump_w/2 - Lz_below) / float(n_flat)
    delta_z2 = (Lz_above - (bump_z+bump_w/2))/ float(n_flat)
    # to revert: Lz_ above in z2 => Lz ; remove lz_below from z1
    
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
        currentZ = delta_z1*i + Lz_below
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
        currentZ = delta_z1*i + Lz_below
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
    
    LzBelow = (int)(bump_z - 3.5 * bump_w)
    LzAbove = (int)(bump_z + 3.5 * bump_w)
    #Lowest plane
    fileCont += '    ({0:8g} {1:8g} {2:8g})\n'.format(self.ellipsoid_function(5, 0, 0), self.ellipsoid_function(5, 0, 1), 0.0)   # 0
    fileCont += '    ({0:8g} {1:8g} {2:8g})\n'.format(self.ellipsoid_function(7, 0, 0), self.ellipsoid_function(7, 0, 1), 0.0)   # 1
    fileCont += '    ({0:8g} {1:8g} {2:8g})\n'.format( 								 -x_rect,  								  -y_rect, 0.0)   # 2
    fileCont += '    ({0:8g} {1:8g} {2:8g})\n'.format( 								  x_rect,  								  -y_rect, 0.0)   # 3
    fileCont += '    ({0:8g} {1:8g} {2:8g})\n'.format( 								 -x_rect,  								   y_rect, 0.0)   # 4
    fileCont += '    ({0:8g} {1:8g} {2:8g})\n'.format( 								  x_rect,  								   y_rect, 0.0)   # 5
    fileCont += '    ({0:8g} {1:8g} {2:8g})\n'.format(self.ellipsoid_function(3, 0, 0), self.ellipsoid_function(3, 0, 1), 0.0)   # 6
    fileCont += '    ({0:8g} {1:8g} {2:8g})\n'.format(self.ellipsoid_function(1, 0, 0), self.ellipsoid_function(1, 0, 1), 0.0)   # 7
    
    #Below Bump
    fileCont += '    ({0:8g} {1:8g} {2:8g})\n'.format(self.ellipsoid_function(5, 0, 0), self.ellipsoid_function(5, 0, 1), LzBelow)   # 8
    fileCont += '    ({0:8g} {1:8g} {2:8g})\n'.format(self.ellipsoid_function(7, 0, 0), self.ellipsoid_function(7, 0, 1), LzBelow)   # 9
    fileCont += '    ({0:8g} {1:8g} {2:8g})\n'.format( 								 -x_rect,  								  -y_rect, LzBelow)   #10
    fileCont += '    ({0:8g} {1:8g} {2:8g})\n'.format( 								  x_rect,  							 	  -y_rect, LzBelow)   #11
    fileCont += '    ({0:8g} {1:8g} {2:8g})\n'.format( 								 -x_rect,  								   y_rect, LzBelow)   #12
    fileCont += '    ({0:8g} {1:8g} {2:8g})\n'.format( 								  x_rect,  								   y_rect, LzBelow)   #13
    fileCont += '    ({0:8g} {1:8g} {2:8g})\n'.format(self.ellipsoid_function(3, 0, 0), self.ellipsoid_function(3, 0, 1), LzBelow)   #14
    fileCont += '    ({0:8g} {1:8g} {2:8g})\n'.format(self.ellipsoid_function(1, 0, 0), self.ellipsoid_function(1, 0, 1), LzBelow)   #15
    
    #Above Bump
    fileCont += '    ({0:8g} {1:8g} {2:8g})\n'.format(self.ellipsoid_function(5, 0, 0), self.ellipsoid_function(5, 0, 1), LzAbove)   #16
    fileCont += '    ({0:8g} {1:8g} {2:8g})\n'.format(self.ellipsoid_function(7, 0, 0), self.ellipsoid_function(7, 0, 1), LzAbove)   #17
    fileCont += '    ({0:8g} {1:8g} {2:8g})\n'.format( 								 -x_rect,  								  -y_rect, LzAbove)   #18
    fileCont += '    ({0:8g} {1:8g} {2:8g})\n'.format( 								  x_rect,  							 	  -y_rect, LzAbove)   #19
    fileCont += '    ({0:8g} {1:8g} {2:8g})\n'.format( 								 -x_rect,  								   y_rect, LzAbove)   #20
    fileCont += '    ({0:8g} {1:8g} {2:8g})\n'.format( 								  x_rect,  								   y_rect, LzAbove)   #21
    fileCont += '    ({0:8g} {1:8g} {2:8g})\n'.format(self.ellipsoid_function(3, 0, 0), self.ellipsoid_function(3, 0, 1), LzAbove)   #22
    fileCont += '    ({0:8g} {1:8g} {2:8g})\n'.format(self.ellipsoid_function(1, 0, 0), self.ellipsoid_function(1, 0, 1), LzAbove)   #23
    
    #Highest plane
    fileCont += '    ({0:8g} {1:8g} {2:8g})\n'.format(self.ellipsoid_function(5, 0, 0), self.ellipsoid_function(5, 0, 1), Lz)   #24
    fileCont += '    ({0:8g} {1:8g} {2:8g})\n'.format(self.ellipsoid_function(7, 0, 0), self.ellipsoid_function(7, 0, 1), Lz)   #25
    fileCont += '    ({0:8g} {1:8g} {2:8g})\n'.format( 								 -x_rect,  								  -y_rect, Lz)   #26
    fileCont += '    ({0:8g} {1:8g} {2:8g})\n'.format( 								  x_rect,  							 	  -y_rect, Lz)   #27
    fileCont += '    ({0:8g} {1:8g} {2:8g})\n'.format( 								 -x_rect,  								   y_rect, Lz)   #28
    fileCont += '    ({0:8g} {1:8g} {2:8g})\n'.format( 								  x_rect,  								   y_rect, Lz)   #29
    fileCont += '    ({0:8g} {1:8g} {2:8g})\n'.format(self.ellipsoid_function(3, 0, 0), self.ellipsoid_function(3, 0, 1), Lz)   #30
    fileCont += '    ({0:8g} {1:8g} {2:8g})\n'.format(self.ellipsoid_function(1, 0, 0), self.ellipsoid_function(1, 0, 1), Lz)   #31

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
    #BelowBump
    fileCont += self.spline_function("ellipse", " 8  9", 5, LzBelow)
    fileCont += self.spline_function("ellipse", "14 15", 1, LzBelow)
    fileCont += self.spline_function("ellipse", " 8 14", 3, LzBelow)
    fileCont += self.spline_function("ellipse", " 9 15", 7, LzBelow)
    #AboveBump
    fileCont += self.spline_function("ellipse", "16 17", 5, LzAbove)
    fileCont += self.spline_function("ellipse", "22 23", 1, LzAbove)
    fileCont += self.spline_function("ellipse", "16 22", 3, LzAbove)
    fileCont += self.spline_function("ellipse", "17 23", 7, LzAbove)
    #Outlet
    fileCont += self.spline_function("ellipse", "24 25", 5, Lz)
    fileCont += self.spline_function("ellipse", "30 31", 1, Lz)
    fileCont += self.spline_function("ellipse", "24 30", 3, Lz)
    fileCont += self.spline_function("ellipse", "25 31", 7, Lz)
    #********************

    #********Bump********
    #Outside
    fileCont += self.spline_function("bump", " 8 16", 5, 0)
    fileCont += self.spline_function("bump", "15 23", 1, 0)
    fileCont += self.spline_function("bump", "14 22", 3, 0)
    fileCont += self.spline_function("bump", " 9 17", 7, 0)
    #Inside
    fileCont += self.spline_function("inside", "10 18", 5, 0)
    fileCont += self.spline_function("inside", "13 21", 1, 0)
    fileCont += self.spline_function("inside", "12 20", 3, 0)
    fileCont += self.spline_function("inside", "11 19", 7, 0)
    #Outside
    #fileCont += self.spline_function("bump", " 0 24", 5, 0)
    #fileCont += self.spline_function("bump", " 7 31", 1, 0)
    #fileCont += self.spline_function("bump", " 6 30", 3, 0)
    #fileCont += self.spline_function("bump", " 1 25", 7, 0)
    #Inside
    #fileCont += self.spline_function("inside", " 2 26", 5, 0)
    #fileCont += self.spline_function("inside", " 5 29", 1, 0)
    #fileCont += self.spline_function("inside", " 4 28", 3, 0)
    #fileCont += self.spline_function("inside", " 3 27", 7, 0)
    #********************

    fileCont += ');\n'

    ######################################################
    # add blocks                                         #
    ######################################################
    fileCont += '\n'
    fileCont += 'blocks\n'
    fileCont += '(\n'

    #####################################################################################
    
    #Grading Calculations
    Lz_bump = 7 * bump_w
    Lz_above = Lz - bump_z - Lz_bump/2
    Lz_below = Lz - Lz_bump - Lz_above
    
    z_res_bump = int(Lz_bump * 7 * z_res / Lz)
    z_G_bump = 1
    cell_width_bump = Lz_bump / z_res_bump
    z_G_above = 1 / cell_width_bump
    z_G_below = cell_width_bump
    
    
    dz_bump = Lz_bump / z_res_bump
    Hx1 = Lz_above/dz_bump
    z_res_above = 1
    sumG = 0.0
    while sumG<Hx1:
      z_res_above += 1
      lam = pow(z_G, 1.0/float(z_res_above-1))
      sumG = 0.0
      for i in range(z_res_above):
        sumG += pow(lam, i)
    z_res_above = int(z_res_above)
    
    Hx1 = Lz_below/dz_bump
    z_res_below = 1
    sumG = 0.0
    while sumG<Hx1:
      z_res_below += 1
      lam = pow(z_G, 1.0/float(z_res_below-1))
      sumG = 0.0
      for i in range(z_res_below):
        sumG += pow(lam, i)
    z_res_below = int(z_res_below)
    
    
    
   

    #block 0
    fileCont += '    hex (0  1  3  2   8   9  11  10) ({0:d} {1:d} {2:d})\n' \
                '      simpleGrading (\n' \
                '        {3:g}\n' \
                '        {4:g}\n' \
                '        {5:g}\n' \
                '      )\n'. \
                format(x_res, y_res, z_res_below, x_G, y_G, z_G_below)
    #block 1            
    fileCont += '    hex (7  6  4  5  15  14  12  13) ({0:d} {1:d} {2:d})\n' \
                '      simpleGrading (\n' \
                '        {3:g}\n' \
                '        {4:g}\n' \
                '        {5:g}\n' \
                '      )\n'. \
                format(x_res, y_res, z_res_below, x_G, y_G, z_G_below)

    #block 2            
    fileCont += '    hex (6  0  2  4  14   8  10  12) ({0:d} {1:d} {2:d})\n' \
                '      simpleGrading (\n' \
                '        {3:g}\n' \
                '        {4:g}\n' \
                '        {5:g}\n' \
                '      )\n'. \
                format(x_res, y_res, z_res_below, x_G, y_G, z_G_below)
                
    #block 3            
    fileCont += '    hex (1  7  5  3   9  15  13  11) ({0:d} {1:d} {2:d})\n' \
                '      simpleGrading (\n' \
                '        {3:g}\n' \
                '        {4:g}\n' \
                '        {5:g}\n' \
                '      )\n'. \
                format(x_res, y_res, z_res_below, x_G, y_G, z_G_below)
                
    #block 4            
    fileCont += '    hex (2  3  5  4  10  11  13  12) ({0:d} {1:d} {2:d})\n' \
                '      simpleGrading (\n' \
                '        {3:g}\n' \
                '        {4:g}\n' \
                '        {5:g}\n' \
                '      )\n'. \
                format(x_res, x_res, z_res_below, 1, 1, z_G_below)
    #block 5
    fileCont += '    hex (8  9 11 10  16  17  19  18) ({0:d} {1:d} {2:d})\n' \
                '      simpleGrading (\n' \
                '        {3:g}\n' \
                '        {4:g}\n' \
                '        {5:g}\n' \
                '      )\n'. \
                format(x_res, y_res, z_res_bump, x_G, y_G, z_G_bump)
    #block 6            
    fileCont += '    hex (15 14 12 13  23  22  20  21) ({0:d} {1:d} {2:d})\n' \
                '      simpleGrading (\n' \
                '        {3:g}\n' \
                '        {4:g}\n' \
                '        {5:g}\n' \
                '      )\n'. \
                format(x_res, y_res, z_res_bump, x_G, y_G, z_G_bump)

    #block 7            
    fileCont += '    hex (14  8 10 12  22  16  18  20) ({0:d} {1:d} {2:d})\n' \
                '      simpleGrading (\n' \
                '        {3:g}\n' \
                '        {4:g}\n' \
                '        {5:g}\n' \
                '      )\n'. \
                format(x_res, y_res, z_res_bump, x_G, y_G, z_G_bump)
                
    #block 8            
    fileCont += '    hex (9 15 13 11  17  23  21  19) ({0:d} {1:d} {2:d})\n' \
                '      simpleGrading (\n' \
                '        {3:g}\n' \
                '        {4:g}\n' \
                '        {5:g}\n' \
                '      )\n'. \
                format(x_res, y_res, z_res_bump, x_G, y_G, z_G_bump)
                
    #block 9            
    fileCont += '    hex (10 11 13 12  18  19  21  20) ({0:d} {1:d} {2:d})\n' \
                '      simpleGrading (\n' \
                '        {3:g}\n' \
                '        {4:g}\n' \
                '        {5:g}\n' \
                '      )\n'. \
                format(x_res, x_res, z_res_bump, 1, 1, z_G_bump)
    #block 10
    fileCont += '    hex (16 17 19 18  24  25  27  26) ({0:d} {1:d} {2:d})\n' \
                '      simpleGrading (\n' \
                '        {3:g}\n' \
                '        {4:g}\n' \
                '        {5:g}\n' \
                '      )\n'. \
                format(x_res, y_res, z_res_above, x_G, y_G, z_G_above)
    #block 11            
    fileCont += '    hex (23 22 20 21 31  30  28  29) ({0:d} {1:d} {2:d})\n' \
                '      simpleGrading (\n' \
                '        {3:g}\n' \
                '        {4:g}\n' \
                '        {5:g}\n' \
                '      )\n'. \
                format(x_res, y_res, z_res_above, x_G, y_G, z_G_above)

    #block 12            
    fileCont += '    hex (22 16 18 20  30  24  26  28) ({0:d} {1:d} {2:d})\n' \
                '      simpleGrading (\n' \
                '        {3:g}\n' \
                '        {4:g}\n' \
                '        {5:g}\n' \
                '      )\n'. \
                format(x_res, y_res, z_res_above, x_G, y_G, z_G_above)
                
    #block 13            
    fileCont += '    hex (17 23 21 19  25  31  29  27) ({0:d} {1:d} {2:d})\n' \
                '      simpleGrading (\n' \
                '        {3:g}\n' \
                '        {4:g}\n' \
                '        {5:g}\n' \
                '      )\n'. \
                format(x_res, y_res, z_res_above, x_G, y_G, z_G_above)
                
    #block 14            
    fileCont += '    hex (18 19 21 20  26  27  29  28) ({0:d} {1:d} {2:d})\n' \
                '      simpleGrading (\n' \
                '        {3:g}\n' \
                '        {4:g}\n' \
                '        {5:g}\n' \
                '      )\n'. \
                format(x_res, x_res, z_res_above, 1, 1, z_G_above)
                
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
    fileCont += '        (24  25  27  26)\n'
    fileCont += '        (30  28  29  31)\n'
    fileCont += '        (24  26  28  30)\n'
    fileCont += '        (25  31  29  27)\n'
    fileCont += '        (26  27  29  28)\n'
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
    
    fileCont += '        (16  17  25  24)\n'
    fileCont += '        (22  23  31  30)\n'
    fileCont += '        (22  16  24  30)\n'
    fileCont += '        (17  23  31  25)\n'
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
