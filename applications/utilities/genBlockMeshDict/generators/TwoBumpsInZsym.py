import os, sys
import numpy as np
import genBlockMesh as gbm

__author__ = 'vstar'

class TwoBumpsInZsym(gbm.AbstractBaseGenerator):
  '''
  This is the generator for the flat fracture with two waves in Z direction
  It is a symmetric case where the bumps are on the opposite walls and directed
  inside of the fracture.
  '''

  def __init__(self):
    self.__genName__ = 'twoBumpsSym'

    # an array of parameters for current generator
    self.parameters = []

    # the format for parameters: [name, initial_value, description
    self.parameters.append(['Lx',  4.0, 'size of the fracture in X direction'])
    self.parameters.append(['Ly',    1.0, 'size of the fracture in Y direction'])
    self.parameters.append(['Lz',  500.0, 'size of the fracture in Z direction'])
    self.parameters.append(['x_res', 4, 'number of cells in X direction'])
    self.parameters.append(['y_res',   8, 'number of cells in Y direction'])
    self.parameters.append(['z_res', 512, 'number of cells in Z direction'])
    self.parameters.append(['x_G',   1.0, 'grading in X direction'])
    self.parameters.append(['y_G',   2.0, 'variable is used for symmetric grading in Y direction'])
    self.parameters.append(['z_G',   10.0, 'variable is used for advanced grading in Z direction'])
    self.parameters.append(['n_spline',    4000, 'number of points for the spline'])
    self.parameters.append(['bump1_z', 100.0, 'position of the first seed'])
    self.parameters.append(['bump2_z', 200.0, 'position of the second seed'])
    self.parameters.append(['bump_A',    0.5, 'amplitude of the wave'])
    self.parameters.append(['bump_w',    5.0, 'wave length'])


  def read_spline_points(self):
    '''
    Reads the edge points and returns them as an array. Returns coordinates for the
    inlet and outlet as well.
    '''
    try:
      file = open('edges', 'r')
    except IOError:
      print('\nError. Cannot open `edges`.')
      sys.exit(2)

    header = file.readline()
    lines = file.readlines()
    posy = []
    negy = []
    for line in lines:
      clrline = line.replace('(', '')
      clrline = clrline.replace(')', '')
      spl_line_fl = [float(x) for x in clrline.split()]

      posy.append([spl_line_fl[0], spl_line_fl[1], spl_line_fl[2]])
      negy.append([spl_line_fl[3], spl_line_fl[4], spl_line_fl[5]])

    inletPpos = posy[0]
    outletPpos = posy[-1]
    inletPneg = negy[0]
    outletPneg = negy[-1]

    return posy, negy, inletPpos, outletPpos, inletPneg, outletPneg

  def bump_function(self, y, z):
    bump1_z = self.get_parameter('bump1_z')
    bump2_z = self.get_parameter('bump2_z')
    A = self.get_parameter('bump_A')
    wl = self.get_parameter('bump_w')

    newY = y
    modifyY = False

    if( z > bump1_z-wl/2. and z < bump1_z+wl/2.):
      rel_z = z - (bump1_z-wl/2.)
      modifyY = True
    elif(z > bump2_z-wl/2. and z < bump2_z+wl/2.):
      rel_z = z - (bump2_z-wl/2.)
      modifyY = True

    if modifyY:
      newY = y - np.sign(y) * A/2.*(1-np.cos(2*np.pi*rel_z/wl))

    return newY

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
    empty, n_spline = self.check_par('n_spline', lines)
    empty, bump1_z = self.check_par('bump1_z', lines)
    empty, bump2_z = self.check_par('bump2_z', lines)
    empty, bump_A = self.check_par('bump_A', lines)
    empty, bump_w = self.check_par('bump_w', lines)

    # points to define inlet and outlet planes
    inletP0 = [-Lx/2., -Ly/2., 0.0]
    inletP1 = [ Lx/2., -Ly/2., 0.0]
    inletP2 = [ Lx/2.,  Ly/2., 0.0]
    inletP3 = [-Lx/2.,  Ly/2., 0.0]
    outletP0 = [-Lx/2., -Ly/2., Lz]
    outletP1 = [ Lx/2., -Ly/2., Lz]
    outletP2 = [ Lx/2.,  Ly/2., Lz]
    outletP3 = [-Lx/2.,  Ly/2., Lz]
    # create edges as arrays
    slp04 = []
    slp15 = []
    slp26 = []
    slp37 = []

    dz = Lz / float(n_spline-1)
    if self.__remesh__:
      posy, negy, inletPpos, outletPpos, inletPneg, outletPneg = self.read_spline_points()
      inletP0 = [-inletPneg[0], inletPneg[1], inletPneg[2]]
      inletP1 = [ inletPneg[0], inletPneg[1], inletPneg[2]]
      inletP2 = [ inletPpos[0], inletPpos[1], inletPpos[2]]
      inletP3 = [-inletPpos[0], inletPpos[1], inletPpos[2]]
      outletP0 = [-outletPneg[0], outletPneg[1], outletPneg[2]]
      outletP1 = [ outletPneg[0], outletPneg[1], outletPneg[2]]
      outletP2 = [ outletPpos[0], outletPpos[1], outletPpos[2]]
      outletP3 = [-outletPpos[0], outletPpos[1], outletPpos[2]]
      for i in range(len(posy)):
        slp04.append([-negy[i][0], negy[i][1], negy[i][2]])
        slp15.append([ negy[i][0], negy[i][1], negy[i][2]])
        slp26.append([ posy[i][0], posy[i][1], posy[i][2]])
        slp37.append([-posy[i][0], posy[i][1], posy[i][2]])
    else:
      for i in range(n_spline):
        currentX = -Lx/2.
        currentY = -Ly/2.
        currentZ = i*dz
        slp04.append([currentX, self.bump_function(currentY,currentZ), currentZ])
        currentX =  Lx/2.
        currentY = -Ly/2.
        currentZ = i*dz
        slp15.append([currentX, self.bump_function(currentY,currentZ), currentZ])
        currentX =  Lx/2.
        currentY =  Ly/2.
        currentZ = i*dz
        slp26.append([currentX, self.bump_function(currentY,currentZ), currentZ])
        currentX = -Lx/2.
        currentY =  Ly/2.
        currentZ = i*dz
        slp37.append([currentX, self.bump_function(currentY,currentZ), currentZ])


    fileCont = '\n'
    fileCont += 'convertToMeters 1.0; \n'
    fileCont += '\n'

    ######################################################
    # add vertices                                       #
    ######################################################
    fileCont += 'vertices\n'
    fileCont += '(\n'

    #Lowest plane
    fileCont += '    ({0:8g} {1:8g} {2:8g})\n'.format(inletP0[0],inletP0[1],inletP0[2]) # 0
    fileCont += '    ({0:8g} {1:8g} {2:8g})\n'.format(inletP1[0],inletP1[1],inletP1[2]) # 1
    fileCont += '    ({0:8g} {1:8g} {2:8g})\n'.format(inletP2[0],inletP2[1],inletP2[2]) # 2
    fileCont += '    ({0:8g} {1:8g} {2:8g})\n'.format(inletP3[0],inletP3[1],inletP3[2]) # 3
    #Highest plane
    fileCont += '    ({0:8g} {1:8g} {2:8g})\n'.format(outletP0[0],outletP0[1],outletP0[2]) # 4
    fileCont += '    ({0:8g} {1:8g} {2:8g})\n'.format(outletP1[0],outletP1[1],outletP1[2]) # 5
    fileCont += '    ({0:8g} {1:8g} {2:8g})\n'.format(outletP2[0],outletP2[1],outletP2[2]) # 6
    fileCont += '    ({0:8g} {1:8g} {2:8g})\n'.format(outletP3[0],outletP3[1],outletP3[2]) # 7

    fileCont += ');\n'

    ######################################################
    # add edges                                          #
    ######################################################
    fileCont += '\n'
    fileCont += 'edges\n'
    fileCont += '(\n'

    fileCont += '    spline 0 4 (\n'
    for spl in slp04:
      fileCont += '        ( {0:f}\t{1:f}\t{2:f} )\n'.format(spl[0], spl[1], spl[2])
    fileCont += '    )\n'

    fileCont += '    spline 1 5 (\n'
    for spl in slp15:
      fileCont += '        ( {0:f}\t{1:f}\t{2:f} )\n'.format(spl[0], spl[1], spl[2])
    fileCont += '    )\n'

    fileCont += '    spline 2 6 (\n'
    for spl in slp26:
      fileCont += '        ( {0:f}\t{1:f}\t{2:f} )\n'.format(spl[0], spl[1], spl[2])
    fileCont += '    )\n'

    fileCont += '    spline 3 7 (\n'
    for spl in slp37:
      fileCont += '        ( {0:f}\t{1:f}\t{2:f} )\n'.format(spl[0], spl[1], spl[2])
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
    y_Gs = '((0.5 0.5 {0:g}) (0.5 0.5 {1:g}))'.format(y_G, 1.0/y_G)

    bump_mid =0.5*(bump1_z+bump2_z)
    Lbump = bump_w * 8.0
    bump1_z1 = bump1_z - Lbump/2.
    bump1_z2 = bump1_z + Lbump/2.
    bump2_z1 = bump2_z - Lbump/2.
    bump2_z2 = bump2_z + Lbump/2.

    per1 = bump1_z1/Lz
    per2 = Lbump/Lz
    per3 = (bump_mid-bump1_z2)/Lz
    per4 = per3
    per5 = Lbump/Lz
    per6 = 1.0 - (per1+per2+per3+per4+per5)

    n_bumps = Lbump*5

    dz_bump = Lbump/n_bumps

    topH = Lz-bump2_z2
    Hx1 = topH/dz_bump
    res_top = 1
    sumG = 0.0
    while sumG<Hx1:
      res_top += 1
      lam = pow(z_G, 1.0/float(res_top-1))
      sumG = 0.0
      for i in range(res_top):
        sumG += pow(lam, i)

    midH = bump2_z1-bump_mid
    Hx1 = midH/dz_bump
    res_mid = 1
    sumG = 0.0
    while sumG<Hx1:
      res_mid += 1
      lam = pow(z_G, 1.0/float(res_mid-1))
      sumG = 0.0
      for i in range(res_mid):
        sumG += pow(lam, i)

    lowH = bump1_z1
    Hx1 = lowH/dz_bump
    res_low = 1
    sumG = 0.0
    while sumG<Hx1:
      res_low += 1
      lam = pow(z_G, 1.0/float(res_low-1))
      sumG = 0.0
      for i in range(res_low):
        sumG += pow(lam, i)

    total_n = 2*n_bumps + res_top + 2*res_mid + res_low

    perC1 = res_low / float(total_n)
    perC2 = n_bumps / float(total_n)
    perC3 = res_mid / float(total_n)
    perC4 = res_mid / float(total_n)
    perC5 = n_bumps / float(total_n)
    perC6 = res_top / float(total_n)

    z_Gs = '(\n' \
          '\t\t({0:g} {1:g} {2:g})\n' \
          '\t\t({3:g} {4:g} {5:g})\n' \
          '\t\t({6:g} {7:g} {8:g})\n' \
          '\t\t({9:g} {10:g} {11:g})\n' \
          '\t\t({12:g} {13:g} {14:g})\n' \
          '\t\t({15:g} {16:g} {17:g})\n' \
          '\t)'.format(
      per1, perC1, 1/z_G,
      per2, perC2, 1,
      per3, perC3, z_G,
      per4, perC4, 1/z_G,
      per5, perC5, 1,
      per6, perC6, z_G,
    )


    #####################################################################################

    fileCont += '    hex (0  1  2  3  4  5  6  7) ({0:d} {1:d} {2:d})\n' \
                '      simpleGrading (\n' \
                '        {3:g}\n' \
                '        {4:s}\n' \
                '        {5:s}\n' \
                '      )\n'. \
                format(x_res, y_res, z_res, x_G, y_Gs, z_Gs)

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
    fileCont += '        (4  5  6  7)\n'
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
    fileCont += '        ( 0   4   7   3)\n'
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
