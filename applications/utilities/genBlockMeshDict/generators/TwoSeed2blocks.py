import os, sys
import numpy as np
import genBlockMesh as gbm

__author__ = 'CmarraJ'

class TwoSeed2blocks(gbm.AbstractBaseGenerator):
  '''
  This generator creates a two seeded fracture. Cosine like seeds are at the inlet at seed_1loc and seed_2loc.
  blockMeshDict consist of two blocks. First block contains the seed and the second
  is a flat fracture.
  '''
  def __init__(self):
    self.__genName__ = 'twoseed2blocks'

    # an array of parameters for current generator
    self.parameters = []

    # the format for parameters: [name, initial_value, description
    self.parameters.append(['Lx',  100.0, 'size of the fracture in X direction'])
    self.parameters.append(['Ly',    1.0, 'size of the fracture in Y direction'])
    self.parameters.append(['Lz',  500.0, 'size of the fracture in Z direction'])
    self.parameters.append(['x_res', 200, 'number of cells in X direction'])
    self.parameters.append(['y_res',   8, 'number of cells in Y direction'])
    self.parameters.append(['z_res', 512, 'number of cells in Z direction'])
    self.parameters.append(['x_Gin',   10.0, 'grading in X direction between seeds'])
    self.parameters.append(['x_Gout',   10.0, 'grading in X direction on the sides'])
    self.parameters.append(['y_G',   1.0, 'grading in Y direction'])
    self.parameters.append(['z_G',   100.0, 'grading in Z direction'])
    self.parameters.append(['nseed',    400, 'number of points for the spline'])
    self.parameters.append(['seed_1loc',    -20.0, 'location of seed 1 in X direction'])
    self.parameters.append(['seed_2loc',    20.0, 'location of seed 2 in X direction'])
    self.parameters.append(['seed_x1',   2.0, 'size of seed 1 in X direction'])
    self.parameters.append(['seed_y1',   0.5, 'size of seed 1 in Y direction'])
    self.parameters.append(['seed_x2',   2.0, 'size of seed 2 in X direction'])
    self.parameters.append(['seed_y2',   0.5, 'size of seed 2 in Y direction'])
    self.parameters.append(['seed_z',   5.0, 'size of the seed in Z direction'])
    self.parameters.append(['x_G_conv',   0.007, 'max cell length difference from seed to side graded regions'])

  def createBlockMeshDict(self, dictFileName):
    lines = self.read_dict(dictFileName)
    empty, Lx = self.check_par('Lx', lines)
    empty, Ly = self.check_par('Ly', lines)
    empty, Lz = self.check_par('Lz', lines)
    empty, x_res = self.check_par('x_res', lines)
    empty, y_res = self.check_par('y_res', lines)
    empty, z_res = self.check_par('z_res', lines)
    empty, x_Gin = self.check_par('x_Gin', lines)
    empty, x_Gout = self.check_par('x_Gout', lines)
    empty, y_G = self.check_par('y_G', lines)
    empty, z_G = self.check_par('z_G', lines)
    empty, n = self.check_par('nseed', lines)
    empty, seed_1loc = self.check_par('seed_1loc', lines)
    empty, seed_2loc = self.check_par('seed_2loc', lines)
    empty, seed_x1 = self.check_par('seed_x1', lines)
    empty, seed_y1 = self.check_par('seed_y1', lines)
    empty, seed_x2 = self.check_par('seed_x2', lines)
    empty, seed_y2 = self.check_par('seed_y2', lines)
    empty, seed_z = self.check_par('seed_z', lines)
    empty, x_G_conv = self.check_par('x_G_conv', lines)

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
      if xi>(-seed_x1/2.+seed_1loc) and xi<(seed_x1/2.+seed_1loc):
        yi -= seed_y1/2.*(1+np.cos(2.*np.pi*(xi-seed_1loc)/seed_x1))
      if xi>(-seed_x2/2.+seed_2loc) and xi<(seed_x2/2.+seed_2loc):
        yi -= seed_y2/2.*(1+np.cos(2.*np.pi*(xi-seed_2loc)/seed_x2))
      fileCont += '        ( {0:f}\t{1:f}\t{2:f} )\n'.format(xi, yi, 0.0)
    fileCont += '    )\n'

    fileCont += '    spline 3 2 (\n'
    for xi in xx:
      yi = Ly/2.
      if xi>(-seed_x1/2.+seed_1loc) and xi<(seed_x1/2.+seed_1loc):
        yi += seed_y1/2.*(1+np.cos(2.*np.pi*(xi-seed_1loc)/seed_x1))
      if xi>(-seed_x2/2.+seed_2loc) and xi<(seed_x2/2.+seed_2loc):
        yi += seed_y2/2.*(1+np.cos(2.*np.pi*(xi-seed_2loc)/seed_x2))
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

    y_G = '((0.5 0.5 2) (0.5 0.5 0.5))'

    # Grading in X direction. We want the seeds to be more detaled
    
    Lseed1 = seed_x1*2
    G_in = 2*x_Gin
    G_out = x_Gout
    L_in1 = abs(seed_1loc) - seed_x1
    L_out1 = abs(Lx/2 - (Lseed1+L_in1))
    
    
    
    cur_in_res1=int((x_res-4)/4)
    cur_out_res1=int((x_res-4)/4)
    cur_seed_res1=x_res/2 - (cur_in_res1 + cur_out_res1)
    
    while True: 
        lam_in = pow(G_in*0.5, 1.0/float(cur_in_res1-1))
        lam_out = pow(G_out, 1.0/float(cur_out_res1-1))
        sum_lam_in = 0.0
        sum_lam_out = 0.0
        for i in range(cur_in_res1):
            sum_lam_in += pow(lam_in, i)
        for j in range(cur_out_res1):
            sum_lam_out += pow(lam_out, j)
            
        xblock_in = L_in1/sum_lam_in
        xblock_out = L_out1/sum_lam_out
        xblock_seed = Lseed1/cur_seed_res1
        
        
        if (xblock_in - xblock_seed) > 0 :
            cur_in_res1 += 1
            cur_seed_res1 -= 1
        if (xblock_out - xblock_seed) > 0 :
            cur_out_res1 += 1
            cur_seed_res1 -= 1
        if (xblock_in - xblock_seed) < 0 :
            cur_in_res1 -= 1
            cur_seed_res1 += 1
        if (xblock_out - xblock_seed) < 0 :
            cur_out_res1 -= 1
            cur_seed_res1 += 1
        if abs(xblock_in + xblock_out - 2*xblock_seed) < x_G_conv:
            cur_in_res1 -= 1
            cur_seed_res1 += 1
            break
        
    
        
    in_side_part1 = cur_in_res1/float(x_res)
    out_side_part1 = cur_out_res1/float(x_res)
    seed_part1 = cur_seed_res1/float(x_res)
    
    Lseed2 = seed_x2*2
    L_in2 = seed_2loc - seed_x2
    L_out2 = Lx/2 - (Lseed2+L_in2)
    
    
    
    cur_in_res2=int((x_res-4)/4)
    cur_out_res2=int((x_res-4)/4)
    cur_seed_res2=x_res/2 - (cur_in_res2 + cur_out_res2)
    
    while True: 
        lam_in = pow(G_in*0.5, 1.0/float(cur_in_res2-1))
        lam_out = pow(G_out, 1.0/float(cur_out_res2-1))
        sum_lam_in = 0.0
        sum_lam_out = 0.0
        for i in range(cur_in_res2):
            sum_lam_in += pow(lam_in, i)
        for j in range(cur_out_res2):
            sum_lam_out += pow(lam_out, j)
            
        xblock_in = L_in2/sum_lam_in
        xblock_out = L_out2/sum_lam_out
        xblock_seed = Lseed2/cur_seed_res2
        
        
        if (xblock_in - xblock_seed) > 0 :
            cur_in_res2 += 1
            cur_seed_res2 -= 1
        if (xblock_out - xblock_seed) > 0 :
            cur_out_res2 += 1
            cur_seed_res2 -= 1
        if (xblock_in - xblock_seed) < 0 :
            cur_in_res2 -= 1
            cur_seed_res2 += 1
        if (xblock_out - xblock_seed) < 0 :
            cur_out_res2 -= 1
            cur_seed_res2 += 1
        if abs(xblock_in + xblock_out - 2*xblock_seed) < x_G_conv:
            cur_in_res2 -= 1
            cur_seed_res2 += 1
            break
        
    
    in_side_part2 = cur_in_res2/float(x_res)
    out_side_part2 = cur_out_res2/float(x_res)
    seed_part2 = cur_seed_res2/float(x_res)

       
    x_G = '(({0:g} {1:g} {2:g}) ({3:g} {4:g} {5:g}) ({6:g} {7:g} {8:g}) ({9:g} {10:g} {11:g}) ({12:g} {13:g} {14:g}) ({15:g} {16:g} {17:g}))'.\
      format(L_out1/Lx, out_side_part1, 1/G_out,
             Lseed1/Lx, seed_part1, 1.0,
             L_in1/Lx, in_side_part1, 0.5*G_in,
             L_in2/Lx, in_side_part2, 1/(0.5*G_in),
             Lseed2/Lx, seed_part2, 1.0,
             L_out2/Lx, out_side_part2, G_out)
    
    #####################################################################################


    fileCont += '    hex (0  1  2  3  4  5  6  7) ({0:d} {1:d} {2:d})\n' \
                '      simpleGrading (\n' \
                '        {3:s}\n' \
                '        {4:s}\n' \
                '        {5:g}\n' \
                '      )\n'. \
                format(x_res, y_res, seed_res, x_G, y_G, 1.0)

    fileCont += '    hex (4  5  6  7 8 9 10 11) ({0:d} {1:d} {2:d})\n' \
                '      simpleGrading (\n' \
                '        {3:s}\n' \
                '        {4:s}\n' \
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

