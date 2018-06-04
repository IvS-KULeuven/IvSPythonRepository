#!/usr/bin/env python
"""
Generate epydoc documentation of IvS repository and insert images.

The documentation will be stored in the subdirectory "doc/html" of the IvS root
directory (i.e. where this file is located). If previously generated
documentation exists, it will be overwriten.

The main page of the documentation can be found in "doc/html/index.html".

Updating the documentation is a simple as rerunning the script, you do not need
to delete the results from the previous run.

For more command-line options, see

    $:> python makedoc.py -h
"""
#-- import necessary modules
import os
import subprocess
import shutil
import glob
import webbrowser
import argparse
from ivs.aux import termtools


#-- skip documentation generation of the following third party modules:
skip = ['lmfit']

#-- do you want to immediately show the contents in your default webbrowser each
#   time the documentation is generated (False, 'current tab' or 'new tab')?
parser = argparse.ArgumentParser(description='Build IvS Python repository documentation')
parser.add_argument('-o','--open',action='store_true',help='Open documentation in webbrowser after creation')
parser.add_argument('-v','--verbose',action='store_true',help='Increase output verbosity in Epydoc')
parser.add_argument('-w','--width',type=int,default=75,help='Width of images (percentage)')
args = parser.parse_args()

#-- remember the current working directory; you can generate the documentation
#   in any directory if you really want to.
this_dir = os.path.dirname(os.path.abspath(__file__))

#-- collect all files and directories for which to generate documentation
output = subprocess.check_output('git ls-files',shell=True)
alldirs = output.strip().split('\n')

#-- we do a quick selection here since we don't want to generate documentation
#   of third-party modules that are installed alongside the repository
alldirs = [ff for ff in alldirs if os.path.splitext(ff)[1]=='.py' and not os.path.basename(ff)[:4]=='test']
for iskip in skip:
    alldirs = [ff for ff in alldirs if not iskip in ff]

#-- build documentation
cmd = 'epydoc --html '+" ".join(alldirs)+\
            ' -o doc/html --parse-only --graph all {}'.format('-v' if args.verbose else '')
print("Building documentation using the command:\n {}".format(cmd))
flag = subprocess.call(cmd,shell=True)

#-- check if all went well
if flag:
    print("Could not execute command, do you have epydoc installed?")
    raise SystemExit

#-- Epydoc cannot insert images in the HTML code itself, we do this manually.
#   for this, we look into the HTML code and replace all ]include figure]
#   occurrences with the HTML code for image insertion. We copy the html file
#   to a temporary copy, change the html code, save it and replace the old file.
if not os.path.isdir('doc/images'):
    print("You don't have an image directory; images are not inserted")
    raise SystemExit

#-- we can't do everything in a generic way just yet, we first treat the
#   exception to the general rule.
'''shutil.move('doc/html/ivs.timeseries.freqanalyse-module.html','doc/html/ivs.timeseries.freqanalyse-module_.html')
ff = open('doc/html/ivs.timeseries.freqanalyse-module_.html','r')
oo = open('doc/html/ivs.timeseries.freqanalyse-module.html','w')
ESC = chr(27)
for line in ff.readlines():
    if r'p = pl.ylim(0,0.018)' in line:
        oo.write(line+'\n\n')
        oo.write(r"<center><img src='../images/timeseries_freqanalyse_01.png' alt='[image example]' /></center>"+'\n\n')
        #-- save cursor, remove line, print and reset cursor
        message = 'Added image to ivs.timeseries.freqanalyse'
        termtools.overwrite_line(message)
    else:
        oo.write(line)
ff.close()
oo.close()
os.remove('doc/html/ivs.timeseries.freqanalyse-module_.html')'''

#-- now run through all the other ones. We need to explicitly treat line breaks.
files = sorted(glob.glob('doc/html/*module.html'))
files+= sorted(glob.glob('doc/html/*class.html'))
for myfile in files:
    shutil.move(myfile,myfile+'_')
    ff = open(myfile+'_','r')
    oo = open(myfile,'w')

    line_break = False
    for line in ff.readlines():

        if ']include figure]' in line or line_break:
            filename = line.split(']')[-2].strip()
            oo.write(r"<center><img src='../images/{}' alt='[image example]' width='{}%'/></center>".format(filename,args.width)+'\n\n')
            #-- save cursor, remove line, print and reset cursor
            message = 'Added image {} to {}\r'.format(filename,myfile)
            termtools.overwrite_line(message)
            oo.write('\n\n')
            line_break = False

        elif ']include' in line:
            line_break = True

        else:
            oo.write(line)

    ff.close()
    oo.close()
    os.remove(myfile+'_')
#-- that's it! Open the documentation if requested
if os.path.isfile(os.path.join(this_dir,'doc/html/index.html')):
    print("HTML documentation is now available in doc/html")
else:
    print("Something went wrong during documentation generation")
    raise SystemExit
if args.open:
    print("Documentation will be opened in your default webbrowser")
    webbrowser.open_new_tab('doc/html/index.html')
