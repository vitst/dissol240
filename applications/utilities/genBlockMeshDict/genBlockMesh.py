# IPython log file

import sys
from abc import ABCMeta, abstractmethod

class AbstractBaseGenerator(object):
  '''
  This is the base class for all blockMeshDict generators
  '''
  __metaclass__ = ABCMeta

  blockMeshFileName = 'blockMeshDict'
  genDictFileName = 'genDict'

  __remesh__ = False

  # Becomes True with -remesh option. Then the remeshing should be implemented.

  def __init__(self):
    pass

  def info(self):
    '''
    Prints info about the variables needed in genDict for the generator
    '''
    prntCont = '\nDescription of the generator:\n'
    prntCont += self.__doc__
    prntCont += '\ngenerator: {0:s}\n'.format(self.__genName__)
    for prm in self.parameters:
      prntCont += prm[0]+':\t'+str(prm[1])+'\t'+prm[2]+'\n'
    print(prntCont)
    return

  ############################################################################
  # read/write
  ############################################################################
  def read_dict(self, dictFileName):
    '''
    Read a file and return string lines
    '''
    lines = ''
    try:
      f = open(dictFileName, 'r')
    except IOError:
      print('\nError. Cannot open `%s`' % dictFileName)
      sys.exit(2)
    else:
      lines = f.readlines()
      f.close()
    return lines

  def writeBlockMesh(self, content):
    self.printHeader()

    file = open(self.blockMeshFileName, 'a')
    file.write( content )
    file.close()
    return
  ############################################################################

  def check_parameters(self, lines):
    '''
    Check whether parameters in generator dictionary are OK
    '''
    check = True

    ### generator name
    loc_generator = [line for line in lines if 'generator' in line]
    if len(loc_generator)==0:
      check = False
      print('\n*** Error. There is no generator name in dictionary file\n')
    elif len(loc_generator)>1:
      check = False
      print('\n*** Error. There more than one generator name in dictionary file\n')
    else:
      genName = loc_generator[0].split()[1]

      if genName!=self.__genName__:
        check = False
        print('\n*** Error. The generator name `%s`'
              'does not correspond to the name from the dictionary `%s`\n' % (self.__genName__, genName))

    for prm in self.parameters:
      #loc_check, val_par = self.check_par(prm[0], prm[1], lines)
      loc_check, val_par = self.check_par(prm[0], lines)
      check = loc_check and check

    return check, self.__genName__

  def run_check(self, dictFileName):
    lines = self.read_dict(dictFileName)
    result, name = self.check_parameters(lines)
    if result:
      print('\nThe parameters in the dictionary `%s` for the generator `%s` are OK.\n'
            % (dictFileName, name))
    else:
      print('\nThere is something wrong in the dictionary `%s` for the generator `%s`.\n' \
            'See comments above.\n' % (dictFileName, name))

    return result

  def check_par(self, parameterName, lines):
    '''
    General function to check existance and format of the parameter in gen dict
    '''

    check = True
    val_param = None

    if not parameterName in [row[0] for row in self.parameters]:
      print('\n*** Error. There is no parameter `%s`\n' % parameterName)
      sys.exit(2)
      check = False
    else:
      par_type = type(self.get_parameter(parameterName))

      loc_param = [line for line in lines if parameterName in line]
      if len(loc_param)==0:
        check = False
        print('\n*** Error. There is no `%s` in the dictionary file\n' % parameterName)
      elif len(loc_param)>1:
        check = False
        print('\n*** Error. There more than one `%s` in the dictionary file\n' % parameterName)
      else:
        try:
          val_param = par_type(loc_param[0].split()[1])
          self.set_parameter(parameterName, val_param)
        except IndexError:
          check = False
          print('\n*** Error. `%s` does not have a value.' % parameterName)
        except ValueError:
          check = False
          print('\n*** Error. `%s` has a wrong number format.' % parameterName)

    return check, val_param

  def get_parameter(self, parameterName):
    ind = [row[0] for row in self.parameters].index(parameterName)
    return self.parameters[ind][1]
  def set_parameter(self, parameterName, value):
    ind = [row[0] for row in self.parameters].index(parameterName)
    self.parameters[ind][1] = value
    return

  def generate(self):
    '''
    Create a generator dictionary
    '''
    file = open(self.genDictFileName, 'w')
    fileCont = 'generator: {0:s}\n'.format(self.__genName__)
    for prm in self.parameters:
      fileCont += prm[0]+':\t'+str(prm[1])+'\t'+prm[2]+'\n'
    file.write( fileCont )
    file.close()

    return

  def printHeader(self):
    file = open(self.blockMeshFileName, 'a')

    h = "/*------------------------------*- C++ -*--------------------------------*\\ \n"
    h = h + "| =========                 |                                             | \n"
    h = h + "| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox       | \n"
    h = h + "|  \\\\    /   O peration     | Version:  2.4.0                             | \n"
    h = h + "|   \\\\  /    A nd           | Web:      www.openfoam.com                  | \n"
    h = h + "|    \\\\/     M anipulation  |                                             | \n"
    h = h + "\\*-----------------------------------------------------------------------*/ \n"

    h = h + "FoamFile \n"
    h = h + "{ \n"
    h = h + "    version     2.0; \n"
    h = h + "    format      ascii; \n"
    h = h + "    class       dictionary; \n"
    h = h + "    object      blockMeshDict; \n"
    h = h + "} \n"
    h = h + "// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * // \n"

    file.write( h )
    file.close()
    return

  @abstractmethod
  def createBlockMeshDict(self):
    '''
    Creates a blockMeshDict. It should return True or False depending on the operation,
    whether it was successful of not.
    '''
    pass






















'''

#Cylindric/Elliptical*******************************************************************

def points5blocks(n, R=1.0, h=0.25, H=500.0, ES=0.5, alpha=1.5):
  '' '
  fracture with round boundary
  R - y direction ellipse radius
  h - distance between tips
  H - height of the channel
  n - number of points
  ES - width added horizontally to circle to create an ellipse
  alpha - x direction ellipse radius
  '' '

  #xy Resolution
  p = 16
  
  #Resolution in the z direction
  z = 1024

  #Outward Grading
  q = 10
  
  #Tangential Grading
  w = 1
  
  #Z grading
  e = 10  
  
  printHeader()
  file = open(fileName, 'a')
  fileCont = '\n'
  fileCont += 'convertToMeters 1.0; \n'
  fileCont += '\n'

  # add vertices
  fileCont += 'vertices\n'
  fileCont += '(\n'

  phi = np.pi / 4
  dphi = np. pi / 2 / n

  fileCont += '    ( %f\t%f\t%f )\n' % (alpha * np.cos(5 * phi), R * np.sin(5 * phi), 0.0)
  fileCont += '    ( %f\t%f\t%f )\n' % (alpha * np.cos(7 * phi), R * np.sin(7 * phi), 0.0)
  fileCont += '    ( %f\t%f\t%f )\n' % (                 -h-ES,                  -h, 0.0)
  fileCont += '    ( %f\t%f\t%f )\n' % (                  h+ES,                  -h, 0.0)
  fileCont += '    ( %f\t%f\t%f )\n' % (                 -h-ES,                   h, 0.0)
  fileCont += '    ( %f\t%f\t%f )\n' % (                  h+ES,                   h, 0.0)
  fileCont += '    ( %f\t%f\t%f )\n' % (alpha * np.cos(3 * phi), R * np.sin(3 * phi), 0.0)
  fileCont += '    ( %f\t%f\t%f )\n' % (    alpha * np.cos(phi),     R * np.sin(phi), 0.0)

  fileCont += '    ( %f\t%f\t%f )\n' % (alpha * np.cos(5 * phi), R * np.sin(5 * phi),   H)
  fileCont += '    ( %f\t%f\t%f )\n' % (alpha * np.cos(7 * phi), R * np.sin(7 * phi),   H)
  fileCont += '    ( %f\t%f\t%f )\n' % (                 -h-ES,                  -h,   H)
  fileCont += '    ( %f\t%f\t%f )\n' % (                  h+ES,                  -h,   H)
  fileCont += '    ( %f\t%f\t%f )\n' % (                 -h-ES,                   h,   H)
  fileCont += '    ( %f\t%f\t%f )\n' % (                  h+ES,                   h,   H)
  fileCont += '    ( %f\t%f\t%f )\n' % (alpha * np.cos(3 * phi), R * np.sin(3 * phi),   H)
  fileCont += '    ( %f\t%f\t%f )\n' % (    alpha * np.cos(phi),     R * np.sin(phi),   H)

  fileCont += ');\n'

  # add edges
  fileCont += '\n'
  fileCont += 'edges\n'
  fileCont += '(\n'

  #*********************************************************************************************************
  
  #inlet

  fileCont += '    spline 0 1 (\n'
  for i in range(0, n+1, 1):
    fileCont += '        ( %f\t%f\t%f )\n' % (alpha * np.cos(5 * phi + i * dphi), R * np.sin(5 * phi + i * dphi), 0.0)
  fileCont += '    )\n'

  fileCont += '    spline 6 7 (\n'
  for i in range(n, -1, -1):
    fileCont += '        ( %f\t%f\t%f )\n' % (alpha * np.cos(phi + i * dphi), R * np.sin(phi + i * dphi), 0.0)
  fileCont += '    )\n'

  fileCont += '    spline 0 6 (\n'
  for i in range(n, -1, -1):
    fileCont += '        ( %f\t%f\t%f )\n' % (alpha * np.cos(3 * phi + i * dphi), R * np.sin(3 * phi + i * dphi), 0.0)
  fileCont += '    )\n'
  
  fileCont += '    spline 1 7 (\n'
  for i in range(0, n+1, 1):
    fileCont += '        ( %f\t%f\t%f )\n' % (alpha * np.cos(7 * phi + i * dphi), R * np.sin(7 * phi + i * dphi), 0.0)
  fileCont += '    )\n'  
  
  #outlet  
  
  fileCont += '    spline 8 9 (\n'
  for i in range(0, n+1, 1):
    fileCont += '        ( %f\t%f\t%f )\n' % (alpha * np.cos(5 * phi + i * dphi), R * np.sin(5 * phi + i * dphi), H)
  fileCont += '    )\n'

  fileCont += '    spline 14 15 (\n'
  for i in range(n, -1, -1):
    fileCont += '        ( %f\t%f\t%f )\n' % (alpha * np.cos(phi + i * dphi), R * np.sin(phi + i * dphi), H)
  fileCont += '    )\n'

  fileCont += '    spline 8 14 (\n'
  for i in range(n, -1, -1):
    fileCont += '        ( %f\t%f\t%f )\n' % (alpha * np.cos(3 * phi + i * dphi), R * np.sin(3 * phi + i * dphi), H)
  fileCont += '    )\n'
  
  fileCont += '    spline 9 15 (\n'
  for i in range(0, n+1, 1):
    fileCont += '        ( %f\t%f\t%f )\n' % (alpha * np.cos(7 * phi + i * dphi), R * np.sin(7 * phi + i * dphi), H)
  fileCont += '    )\n' 

  fileCont += ');\n'

  # add blocks
  fileCont += '\n'
  fileCont += 'blocks\n'
  fileCont += '(\n'

  fileCont += '    hex (0  1  3  2   8   9  11  10) (%d %d %d) simpleGrading (%f %f %f)\n' % (p, p, z, w, q, e)
  fileCont += '    hex (7  6  4  5  15  14  12  13) (%d %d %d) simpleGrading (%f %f %f)\n' % (p, p, z, w, q, e)

  fileCont += '    hex (6  0  2  4  14   8  10  12) (%d %d %d) simpleGrading (%f %f %f)\n' % (p, p, z, w, q, e)
  fileCont += '    hex (1  7  5  3   9  15  13  11) (%d %d %d) simpleGrading (%f %f %f)\n' % (p, p, z, w, q, e)

  fileCont += '    hex (2  3  5  4  10  11  13  12) (%d %d %d) simpleGrading (1 1 %f)\n' % (p, p, z, e)

  fileCont += ');\n'
  # *****************************************************************************

  # add boundaries
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

  # add boundaries
  fileCont += '\n'
  fileCont += 'mergePatchPairs\n'
  fileCont += '(\n'
  fileCont += ');\n\n'

  fileCont += '// ********************************************************************* //\n'


  file.write( fileCont )
  file.close()

  return

#***************************************************************************************

#Curved*********************************************************************************
def curvedMeshParameters(n, R=0.25, H=500.0, a=0, b=0, m=5, f=1):
  '' '
  fracture with curved boundary
  R - amplitude of curved boundary
  n - number of points on the boundary
  H - height of the channel [z dir]
  a - half of Y direction length (average amplitude of the fracture)
  b - half of X direction length
  m - number of waves on curved boundary
  f - -1 for asymmetry. +1 for symmetry about Z.
 
  f = -1 ********* f = 1
  
   |__|  ********* |__|
   (__(  *********(____)
    )__) ********* )__(  ******^
   (__(  *********(____) ******|
   |__|  ********* |__|  ******z *** y->
  '' '
  
 
  #x Resolution  
  p = 4
  #y Resolution
  o = 8  
  #Resolution in the z direction
  z = 4096
  
  smallZ = z / 250
  largeZ = z - 2 * smallZ

  #Y Grading
  q = 1
  
  #X Grading
  w = 1
  
  #Z grading
  e = 1

  
  
  printHeader()
  file = open(fileName, 'a')
  fileCont = '\n'
  fileCont += 'convertToMeters 1.0; \n'
  fileCont += '\n'

  # add vertices
  fileCont += 'vertices\n'
  fileCont += '(\n'

  # ome = np.pi / 100
  L = H / 250
  curveH = H - 2 * L
  ome = np.pi * 2 * m / curveH 
  phi = np.pi / 2
  dphi = curveH / float(n)
  print "curveH= %f ome= %f dphi= %f" % (curveH, ome, dphi)
  
  
  y1 =  0.5*(1 - a*np.cos(2*m*np.pi*0/float(n)))
  print "y1 = %f" % (y1)
  
  y2 = -0.5*(1 - f*a*np.cos(2*m*np.pi*0/float(n)))
  print "y2 = %f" % (y2)
  
  #Lowest plane
  fileCont += '    ( %f\t%f\t%f )\n' % ( -b, y2, 0.0)
  fileCont += '    ( %f\t%f\t%f )\n' % (  b, y2, 0.0)
  fileCont += '    ( %f\t%f\t%f )\n' % (  b, y1, 0.0)
  fileCont += '    ( %f\t%f\t%f )\n' % ( -b, y1, 0.0)
  #Start of Curve
  fileCont += '    ( %f\t%f\t%f )\n' % ( -b, y2, L)
  fileCont += '    ( %f\t%f\t%f )\n' % (  b, y2, L)
  fileCont += '    ( %f\t%f\t%f )\n' % (  b, y1, L)
  fileCont += '    ( %f\t%f\t%f )\n' % ( -b, y1, L)
  #End of Curve
  fileCont += '    ( %f\t%f\t%f )\n' % ( -b, y2, H - L)
  fileCont += '    ( %f\t%f\t%f )\n' % (  b, y2, H - L)
  fileCont += '    ( %f\t%f\t%f )\n' % (  b, y1, H - L)
  fileCont += '    ( %f\t%f\t%f )\n' % ( -b, y1, H - L)
  #Highest plane
  fileCont += '    ( %f\t%f\t%f )\n' % ( -b, y2, H)
  fileCont += '    ( %f\t%f\t%f )\n' % (  b, y2, H)
  fileCont += '    ( %f\t%f\t%f )\n' % (  b, y1, H)
  fileCont += '    ( %f\t%f\t%f )\n' % ( -b, y1, H)

  fileCont += ');\n'

  # add edges
  fileCont += '\n'
  fileCont += 'edges\n'
  fileCont += '(\n'

  #*********************************************************************************************************
  
  fileCont += '    spline 4 8 (\n'
  for i in range(0, n+1, 1):
    fileCont += '        ( %f\t%f\t%f )\n' % (-b, -0.5*(1 - f*a*np.cos(2*m*np.pi*i/float(n))), dphi * i + L)
  fileCont += '    )\n'
  
  fileCont += '    spline 5 9 (\n'
  for i in range(0, n+1, 1):
    fileCont += '        ( %f\t%f\t%f )\n' % ( b, -0.5*(1 - f*a*np.cos(2*m*np.pi*i/float(n))), dphi * i + L)
  fileCont += '    )\n'
  
  fileCont += '    spline 6 10 (\n'
  for i in range(0, n+1, 1):
    fileCont += '        ( %f\t%f\t%f )\n' % ( b, 0.5*(1 - a*np.cos(2*m*np.pi*i/float(n))), dphi * i + L)
  fileCont += '    )\n'
  
  fileCont += '    spline 7 11 (\n'
  for i in range(0, n+1, 1):
    fileCont += '        ( %f\t%f\t%f )\n' % (-b, 0.5*(1 - a*np.cos(2*m*np.pi*i/float(n))), dphi * i + L)
  fileCont += '    )\n'
  
  
  fileCont += ');\n'

  # add blocks
  fileCont += '\n'
  fileCont += 'blocks\n'
  fileCont += '(\n'
  
  fileCont += '    hex (0  1  2  3  4  5  6  7) (%d %d %d) simpleGrading (%f %f %f)\n' % (p, o, smallZ, w, q, e)
  fileCont += '    hex (4  5  6  7  8  9 10 11) (%d %d %d) simpleGrading (%f %f %f)\n' % (p, o, largeZ, w, q, e)
  fileCont += '    hex (8  9 10 11 12 13 14 15) (%d %d %d) simpleGrading (%f %f %f)\n' % (p, o, smallZ, w, q, e)

  fileCont += ');\n'
  # *****************************************************************************

  # add boundaries
  fileCont += '\n'
  fileCont += 'boundary\n'
  fileCont += '(\n'

  fileCont += '    outlet\n'
  fileCont += '    {\n'
  fileCont += '      type patch;\n'
  fileCont += '      faces\n'
  fileCont += '      (\n'
  fileCont += '        (12  13  14  15)\n'
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
  fileCont += '        ( 4   5   9   8)\n'  
  fileCont += '        ( 8   9  13  12)\n'
    
  fileCont += '        ( 2   3   7   6)\n'
  fileCont += '        ( 6   7  11  10)\n'
  fileCont += '        (10  11  15  14)\n'
  fileCont += '      );\n'
  fileCont += '    }\n'
  fileCont += '\n'
  
  fileCont += '    periodicx1\n'
  fileCont += '    {\n'
  fileCont += '      type cyclic;\n'
  fileCont += '      neighbourPatch periodicx2;\n'
  fileCont += '      faces\n'
  fileCont += '      (\n'
  #fileCont += '        ( 4   7   3   0)\n'
  #fileCont += '        ( 8  11   7   4)\n'
  #fileCont += '        (12  15  11   8)\n'
  fileCont += '        ( 0   4   7   3)\n'
  fileCont += '        ( 4   8  11   7)\n'
  fileCont += '        ( 8  12  15  11)\n'
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
  fileCont += '        ( 9  10  14  13)\n'
  fileCont += '      );\n'
  fileCont += '    }\n'
  fileCont += '\n'

  fileCont += ');\n'
  # *****************************************************************************

  # add boundaries
  fileCont += '\n'
  fileCont += 'mergePatchPairs\n'
  fileCont += '(\n'
  fileCont += ');\n\n'

  fileCont += '// ********************************************************************* //\n'


  file.write( fileCont )
  file.close()

  return

#***************************************************************************************

#Seeded*********************************************************************************
def seeded2blocks(n, A=1.0, H=5000.0, xL=500.0):
  # fracture with round boundary
  # n - number of points
  # A - amplitude of seed (A = 0.25 will result in a seed that goes from y=+-0.5 to y=+-1.0)
  # H - height of the channel
  # xL - 1/2 of total x directional length (length in each direction from origin)
    
  #Directional Lengths
  yL = 0.5
  h = H / 100
  x_seed = 5
  x_frac = xL - x_seed
  
  #Directional Res
  xR = 100
  yR = 16
  zR_1 = 40
  zR_2 = 984
  
  #Directional Grad
  xG_1 = 10.0
  xG_2 = 1/xG_1
  yG_1 = 2.0
  yG_2 = 1/yG_1
  zG_1 = 1.0
  zG_2 = 10.0
  
  #spline eqn
  n_frac = float(n) * ( xL -h ) / xL
  n_seed = float(n) - n_frac  
  W = 10
  x = np.linspace(-0.5*W,0.5*W,n_seed+1)
  y = 0.5 + A*(1+np.cos(2*np.pi*x/W))
  
  printHeader()
  file = open(fileName, 'a')
  fileCont = '\n'
  fileCont += 'convertToMeters 1.0; \n'
  fileCont += '\n'

  # add vertices
  fileCont += 'vertices\n'
  fileCont += '(\n'

  phi = n_frac/2
  dphi = xL/phi
  ome = x_seed/dphi

  fileCont += '    ( %f\t%f\t%f )\n' % (-xL, -yL,   0.0)
  fileCont += '    ( %f\t%f\t%f )\n' % ( xL, -yL,   0.0)
  fileCont += '    ( %f\t%f\t%f )\n' % (-xL,  yL,   0.0)
  fileCont += '    ( %f\t%f\t%f )\n' % ( xL,  yL,   0.0)
  fileCont += '    ( %f\t%f\t%f )\n' % (-xL, -yL,   h)
  fileCont += '    ( %f\t%f\t%f )\n' % ( xL, -yL,   h)
  fileCont += '    ( %f\t%f\t%f )\n' % (-xL,  yL,   h)
  fileCont += '    ( %f\t%f\t%f )\n' % ( xL,  yL,   h)
  fileCont += '    ( %f\t%f\t%f )\n' % (-xL, -yL,   H)
  fileCont += '    ( %f\t%f\t%f )\n' % ( xL, -yL,   H)
  fileCont += '    ( %f\t%f\t%f )\n' % (-xL,  yL,   H)
  fileCont += '    ( %f\t%f\t%f )\n' % ( xL,  yL,   H)

  fileCont += ');\n'

  # add edges
  fileCont += '\n'
  fileCont += 'edges\n'
  fileCont += '(\n'

  #*********************************************************************************************************
  
  #inlet

  fileCont += '    spline 0 1 (\n'
  for i in range(-phi, -ome+1, 1):
    fileCont += '        ( %f\t%f\t%f )\n' % (dphi * i, -yL, 0.0)
  for i in range(0, n_seed+1, 1):
    fileCont += '        ( %f\t%f\t%f )\n' % (x[i], -y[i], 0.0)
  for i in range(ome, phi+1, 1):
    fileCont += '        ( %f\t%f\t%f )\n' % (dphi * i, -yL, 0.0)
  fileCont += '    )\n'

  fileCont += '    spline 2 3 (\n'
  for i in range(-phi, -ome+1, 1):
    fileCont += '        ( %f\t%f\t%f )\n' % (dphi * i, yL, 0.0)
  for i in range(0, n_seed+1, 1):
    fileCont += '        ( %f\t%f\t%f )\n' % (x[i], y[i], 0.0)
  for i in range(ome, phi+1, 1):
    fileCont += '        ( %f\t%f\t%f )\n' % (dphi * i, yL, 0.0)
  fileCont += '    )\n'

  fileCont += ');\n'

  # add blocks
  fileCont += '\n'
  fileCont += 'blocks\n'
  fileCont += '(\n'


  fileCont += '   hex (0  1  3  2  4  5  7  6) (%d %d %d)\n' \
              '   simpleGrading (\n' % (xR, yR, zR_1)
  fileCont += '       (\n' \
              '             (0.5 0.5 %f)\n' \
              '             (0.5 0.5 %f)\n' \
              '       )\n' % (xG_2, xG_1)
  fileCont += '       (\n' \
              '             (0.5 0.5 %f)\n' \
              '             (0.5 0.5 %f)\n' \
              '       )\n' % (yG_1, yG_2)
  fileCont += '       %d \n' \
              '   )\n' % (zG_1)
  
  fileCont += '   hex (4  5  7  6  8  9 11 10) (%d %d %d)\n' \
              '   simpleGrading (\n' % (xR, yR, zR_2)
  fileCont += '       (\n' \
              '             (0.5 0.5 %f)\n' \
              '             (0.5 0.5 %f)\n' \
              '       )\n' % (xG_2, xG_1)
  fileCont += '       (\n' \
              '             (0.5 0.5 %f)\n' \
              '             (0.5 0.5 %f)\n' \
              '       )\n' % (yG_1, yG_2)
  fileCont += '       %d \n' \
              '   )\n' % (zG_2)
              
  fileCont += ');\n'
  # *****************************************************************************

  # add boundaries
  fileCont += '\n'
  fileCont += 'boundary\n'
  fileCont += '(\n'

  fileCont += '    outlet\n'
  fileCont += '    {\n'
  fileCont += '      type patch;\n'
  fileCont += '      faces\n'
  fileCont += '      (\n'
  fileCont += '        ( 8  10  11  9)\n'
  fileCont += '      );\n'
  fileCont += '    }\n'
  fileCont += '\n'

  fileCont += '    inlet\n'
  fileCont += '    {\n'
  fileCont += '      type patch;\n'
  fileCont += '      faces\n'
  fileCont += '      (\n'
  fileCont += '        ( 0   1   3   2)\n'
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
  fileCont += '        ( 0   2   6   4)\n'
  fileCont += '        ( 4   6  10   8)\n'
  fileCont += '      );\n'
  fileCont += '    }\n'
  fileCont += '\n'
  
  fileCont += '    periodicx2\n'
  fileCont += '    {\n'
  fileCont += '      type cyclic;\n'
  fileCont += '      neighbourPatch periodicx1;\n'
  fileCont += '      faces\n'
  fileCont += '      (\n'
  fileCont += '        ( 1   3   7   5)\n'
  fileCont += '        ( 5   7  11   9)\n'
  fileCont += '      );\n'
  fileCont += '    }\n'
  fileCont += '\n'

  fileCont += ');\n'
  # *****************************************************************************

  # add boundaries
  fileCont += '\n'
  fileCont += 'mergePatchPairs\n'
  fileCont += '(\n'
  fileCont += ');\n\n'

  fileCont += '// ********************************************************************* //\n'


  file.write( fileCont )
  file.close()

  return
  
#***************************************************************************************

def curve_function(z, y, A, wl):
  return y - A*(1-np.cos(2*np.pi*z/wl))



def curveMeshParameters(n, H, R):
  # fracture with curved boundary
  # All Values set by runCurvedGen.py [these are just dummy var]
  # R - amplitude of curved boundary
  # n - number of points on the boundary [for the spline]
  # H - height of the channel [z dir]
  #

  dz = H/float(n-1)

  # x size
  x_size = 2.0
  # y size
  y_size = 1.0
  # z size is H

  # x Resolution
  x_resol = 4
  # y Resolution
  y_resol = 16
  # Resolution in the z direction
  z_resol = 1024 #8192

  # grading in x direction
  x_G = 1.0
  # grading qin y direction
  y_G = 1.0
  # grading in z direction
  z_G = 1.0


  # the position of bumps
  z_1 = 100
  z_2 = 200
  wave_length = 25
  amplitude = R/2.


  printHeader()
  file = open(fileName, 'a')
  fileCont = '\n'
  fileCont += 'convertToMeters 1.0; \n'
  fileCont += '\n'

  # add vertices
  fileCont += 'vertices\n'
  fileCont += '(\n'



  # ome = np.pi / 100
  '' '
  L = H / 250
  curveH = H - 2 * L
  ome = np.pi * 2 * m / curveH
  phi = np.pi / 2
  dphi = curveH / float(n)
  print "curveH= %f ome= %f dphi= %f" % (curveH, ome, dphi)
  y1 =  0.5*(1 - a*np.cos(2*m*np.pi*0/float(n)))
  print "y1 = %f" % (y1)
  y2 = -0.5*(1 - f*a*np.cos(2*m*np.pi*0/float(n)))
  print "y2 = %f" % (y2)
  '' '


  #Lowest plane
  fileCont += '    ( %f\t%f\t%f )\n' % ( -x_size/2., -y_size/2., 0.0)   # 0
  fileCont += '    ( %f\t%f\t%f )\n' % (  x_size/2., -y_size/2., 0.0)   # 1
  fileCont += '    ( %f\t%f\t%f )\n' % (  x_size/2.,  y_size/2., 0.0)   # 2
  fileCont += '    ( %f\t%f\t%f )\n' % ( -x_size/2.,  y_size/2., 0.0)   # 3
  #Highest plane
  fileCont += '    ( %f\t%f\t%f )\n' % ( -x_size/2., -y_size/2., H)     # 4
  fileCont += '    ( %f\t%f\t%f )\n' % (  x_size/2., -y_size/2., H)     # 5
  fileCont += '    ( %f\t%f\t%f )\n' % (  x_size/2.,  y_size/2., H)     # 6
  fileCont += '    ( %f\t%f\t%f )\n' % ( -x_size/2.,  y_size/2., H)     # 7

  fileCont += ');\n'

  # add edges
  fileCont += '\n'
  fileCont += 'edges\n'
  fileCont += '(\n'

  #*********************************************************************************************************

  fileCont += '    spline 0 4 (\n'
  for i in range(0, n):
    currentX = -x_size/2.
    currentY = -y_size/2.
    currentZ = i*dz
    if( currentZ > z_1-wave_length/2. and currentZ < z_1+wave_length/2.):
      rel_z = currentZ - (z_1-wave_length/2.)
      currentY = curve_function(rel_z, currentY, amplitude, wave_length)
    elif(currentZ > z_2-wave_length/2. and currentZ < z_2+wave_length/2.):
      rel_z = currentZ - (z_2-wave_length/2.)
      currentY = curve_function(rel_z, currentY, amplitude, wave_length)
    else:
      pass

    fileCont += '        ( %f\t%f\t%f )\n' % (currentX, currentY, currentZ)

  fileCont += '    )\n'

  fileCont += '    spline 1 5 (\n'
  for i in range(0, n):
    currentX = x_size/2.
    currentY = -y_size/2.
    currentZ = i*dz
    if( currentZ > z_1-wave_length/2. and currentZ < z_1+wave_length/2.):
      rel_z = currentZ - (z_1-wave_length/2.)
      currentY = curve_function(rel_z, currentY, amplitude, wave_length)
    elif(currentZ > z_2-wave_length/2. and currentZ < z_2+wave_length/2.):
      rel_z = currentZ - (z_2-wave_length/2.)
      currentY = curve_function(rel_z, currentY, amplitude, wave_length)
    else:
      pass

    fileCont += '        ( %f\t%f\t%f )\n' % (currentX, currentY, currentZ)

  fileCont += '    )\n'

  fileCont += '    spline 2 6 (\n'
  for i in range(0, n):
    currentX = x_size/2.
    currentY = y_size/2.
    currentZ = i*dz

    if( currentZ > z_1-wave_length/2. and currentZ < z_1+wave_length/2.):
      rel_z = currentZ - (z_1-wave_length/2.)
      currentY = curve_function(rel_z, currentY, amplitude, wave_length)
    elif(currentZ > z_2-wave_length/2. and currentZ < z_2+wave_length/2.):
      rel_z = currentZ - (z_2-wave_length/2.)
      currentY = curve_function(rel_z, currentY, amplitude, wave_length)
    else:
      pass

    fileCont += '        ( %f\t%f\t%f )\n' % (currentX, currentY, currentZ)

  fileCont += '    )\n'

  fileCont += '    spline 3 7 (\n'
  for i in range(0, n):
    currentX = -x_size/2.
    currentY = y_size/2.
    currentZ = i*dz

    if( currentZ > z_1-wave_length/2. and currentZ < z_1+wave_length/2.):
      rel_z = currentZ - (z_1-wave_length/2.)
      currentY = curve_function(rel_z, currentY, amplitude, wave_length)
    elif(currentZ > z_2-wave_length/2. and currentZ < z_2+wave_length/2.):
      rel_z = currentZ - (z_2-wave_length/2.)
      currentY = curve_function(rel_z, currentY, amplitude, wave_length)
    else:
      pass


    fileCont += '        ( %f\t%f\t%f )\n' % (currentX, currentY, currentZ)

  fileCont += '    )\n'

  # 0.5*(1 - a*np.cos(2*m*np.pi*i/float(n)))


  fileCont += ');\n'

  # add blocks
  fileCont += '\n'
  fileCont += 'blocks\n'
  fileCont += '(\n'

  fileCont += '   hex (0  1  2  3  4  5  6  7) (%d %d %d)\n' \
              '   simpleGrading (\n' \
              '       %f\n' % (x_resol, y_resol, z_resol, x_G)

  fileCont += '       (\n' \
              '             (0.5 1.0 10.0)\n' \
              '             (0.5 1.0  0.1)\n' \
              '       )\n'

  fileCont += '       %f\n' \
              '   )\n' % (z_G)

  fileCont += ');\n'
  # *****************************************************************************

  # add boundaries
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
  # *****************************************************************************

  # add boundaries
  fileCont += '\n'
  fileCont += 'mergePatchPairs\n'
  fileCont += '(\n'
  fileCont += ');\n\n'

  fileCont += '// ********************************************************************* //\n'


  file.write( fileCont )
  file.close()


 # plt.plot( x, y, 'ro')
 # plt.show()

  return

'''


