#!/usr/bin/env python3.4

import os, sys, getopt

listOfAllGenerators = {}

package_dir = 'generators'
package_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), package_dir)

# import generator classes from generators directory
for filename in os.listdir(package_path):
  modulename, ext = os.path.splitext(filename)

  if modulename[0] != '_' and ext == '.py':
    subpackage = '{0}.{1}'.format(package_dir, modulename)
    obj = getattr(__import__(subpackage, globals(), locals(), [modulename]), modulename)
    listOfAllGenerators.update({obj().__genName__ : obj})

short_description="\nThis tool helps to create a `blockMeshDict` file for the OpenFoam's blockMesh."

description='''
  {0:s}\n
  This is a collection of `blockMeshDict` generators which are used to create sophisticated
  geometries for the simple OpenFOAM tool `blockMesh`. A dictionary file is needed.
  Usual structure of a dictionary file:

    name of the generator(generator: *name*)
    list of parameters (*name*: *value*)
    .....

  To see an example of a dictionary file for simple flat fracture generator run:
    `runGen.py -i simple`

  To list all available generators run:
    `runGen.py -l` of `runGen.py --list`

  !Important! Since grading is a complicated parameter usually only simpleGrading is implemented.
  Advanced grading can be implemented within a generator like it is done in seeded2blocks.
  Check `Seeded2blocks.py`.
'''.format(short_description)

package_name="blockMesh dictionary generator"
__doc__ = description
__version__ = 0.1


def help(short):
  if short:
    print(short_description)
  else:
    print(description)
  print('options:\n'
        '\t-h\t\t\tshort help\n'
        '\t--help\t\t\textended help\n'
        '\t-l [or --list]\t\t\tlist of all possible generators\n'
        '\t-i [or --info]  <generatorName>\tprint description of all parameters for the generator\n'
        '\t-c [or --check] <fileName>\tcheck parameters in generator dictionary file\n'
        '\t-d [or --dfile] <fileName>\tgenerate `blockMeshDict` using generator dictionary file\n'
        '\t-g <generatorName>\tgenerate `genDict` using generatorName\n'
        'examples:\n'
        '\trunGen.py -d <name_of_the_dictionary_file>\n'
        '\trunGen.py -g <name_of_the_generator>\n')


def generate(generatorName):
  if os.path.isfile('genDict'):
    choice = ''
    while choice!='y' and choice!='n':
      # need to know what is the version of the interpreter
      if sys.version_info[0]==2:
        choice=raw_input('\nThe file `genDict` exists already. Do you want to delete the file first? (y/n): ')
      else:
        choice=input('\nThe file `genDict` exists already. Do you want to delete the file first? (y/n): ')
    if choice=='y':
      os.remove('genDict')
    else:
      return

  try:
    genClass = listOfAllGenerators[generatorName]
  except KeyError:
    print('\n*** Error. There is no generator `{0:s}`.\n'
          'Run `runGen.py -l` or `runGen.py --list` to see the list of '
          'all available generators.\n'.format(generatorName))
    sys.exit(2)

  generator = genClass()
  generator.generate()

  print('\nA template dictionary file `genDict` for the generator {0:s} is created\n'.format(generatorName))
  return

def check_parameters(dictName):
  try:
    f = open(dictName, 'r')
  except IOError:
    print('\n*** Error. Cannot open `%s`\n' % dictName)
  else:
    generatorLine = f.readline()
    f.close()
    if 'generator:' in generatorLine:
      generatorName = generatorLine.split()[1]
      try:
        genClass = listOfAllGenerators[generatorName]
      except KeyError:
        print('\n*** Error. There is no generator `%s`.\n'
              'Change generator name in your dictionary file.\n'
              'Run `runGen.py -l` or `runGen.py --list` to see the list of '
              'all available generators.\n' % generatorName)
        sys.exit(2)

      generator = genClass()
      if not generator.run_check(dictName):
        sys.exit(2)
    else:
      print('The name of the generator in `%s` was not found.\n'
            'You can always generate dictionary file by running:\n'
            ' runGen.py -g *nameOfAGenerator*' % dictName)
  return

def main(argv):
  try:
    opts, args = getopt.getopt(argv,"hvld:g:i:c:",["help", "version", "dfile=", "list", "info=", "check=", "remesh"])
  except getopt.GetoptError:
    print('\n*** Error. Options are not correct.')
    help(True)
    sys.exit(2)
  
  if len(opts)==0:
    print('for help run -h option')

  for opt, arg in opts:
    if opt == "-h":
      help(True)
      sys.exit()
    elif opt =="--help":
      help(False)
      sys.exit()
    elif opt in ("-v", "--version"):
      print("{0:s}. Version {1:g}".format(package_name, __version__))
      sys.exit()
    elif opt in ("-l", "--list"):
      print('\nList of the `blockMeshDict` generatorors: ')
      #for key, value in listOfAllGenerators.items():
      for key in sorted(listOfAllGenerators.keys()):
        print('  %s' % key)
      print('')
      sys.exit()
    elif opt == "-g":
      generate(arg)
      sys.exit()
    elif opt in ("-i", "--info"):
      try:
        genClass = listOfAllGenerators[arg]
      except KeyError:
        print('\n*** Error. There is no generator `%s`.\n'
              'Run `runGen.py -l` or `runGen.py --list` to see the list of '
              'all available generators.\n' % arg)
        sys.exit(2)
      generator = genClass()
      generator.info()

    elif opt in ("-c", "--check"):
      check_parameters(arg)
    elif opt in ("-d", "--dfile"):
      check_parameters(arg)
      try:
        f = open(arg, 'r')
      except IOError:
        print('\n*** Error. Can not open `%s`\n' % arg)
      else:
        generatorLine = f.readline()
        f.close()
        if 'generator:' in generatorLine:
          generatorName = generatorLine.split()[1]
          try:
            genClass = listOfAllGenerators[generatorName]
          except KeyError:
            print('There is no generator %s. Choose one from the list...' % generatorName)
            sys.exit(2)

          # create a generator
          generator = genClass()
          # check whether a file exists
          if os.path.isfile(generator.blockMeshFileName):
            choice = ''
            while choice!='y' and choice!='n':
              if sys.version_info[0]==2:
                choice=raw_input('The file `{0:s}` exists already. Do you want to delete the file first? (y/n): '.format(generator.blockMeshFileName) )
              else:
                choice=input('The file `{0:s}` exists already. Do you want to delete the file first? (y/n): '.format(generator.blockMeshFileName) )
            if choice=='y':
              os.remove(generator.blockMeshFileName)
            else:
              return
          for opt1, arg1 in opts:
            if opt1 == "--remesh":
              generator.__remesh__ = True

          generated = generator.createBlockMeshDict(arg)

          if not generated:
            print('\nSomething went wrong during the generation :(\n')
            sys.exit(2)

          print('\n`{0:s}` successfully generated!\n'.format(generator.blockMeshFileName))
        else:
          print('The name of the generator in %s was not found.\n'
                'You can always generate dictionary file by running:\n'
                ' runGen.py -g *nameOfAGenerator*' % arg)

if __name__ == "__main__":
  main(sys.argv[1:])

