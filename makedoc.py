#!/usr/bin/env python
"""
Generate epydoc documentation of IVS repository and insert images.
"""
# -*- coding: utf-8 -*-
import os
import time
import sys
import commands
import shutil
import glob

this_dir = os.getcwd()

#-- collect all files and directories for which to generate documentation

stat,output = commands.getstatusoutput('git ls-files')
alldirs = output.split('\n')
os.chdir(this_dir)
alldirs = [ff for ff in alldirs if os.path.splitext(ff)[1]=='.py' and not os.path.basename(ff)[:4]=='test']
alldirs = [ff for ff in alldirs if not 'uncertainties' in ff]
alldirs = [ff for ff in alldirs if not 'lmfit' in ff]

#-- build documentation
print "Building documentation:"
cmd = 'epydoc --html '+" ".join(alldirs)+\
            ' -o doc --parse-only --graph all'# -v'
print cmd
os.system(cmd)


#-- add images
if not os.path.isdir('doc/images'):
    print "You don't have an image directory; images are not inserted"
    raise SystemExit


shutil.move('doc/ivs.timeseries.freqanalyse-module.html','doc/ivs.timeseries.freqanalyse-module_.html')
ff = open('doc/ivs.timeseries.freqanalyse-module_.html','r')
oo = open('doc/ivs.timeseries.freqanalyse-module.html','w')
for line in ff.readlines():
    if r'p = pl.ylim(0,0.018)' in line:
        oo.write(line+'\n\n')
        oo.write(r"<center><img src='images/timeseries_freqanalyse_01.png' alt='[image example]' /></center>"+'\n\n')
        print 'Added image to ivs.timeseries.freqanalyse'
    else:
        oo.write(line)
ff.close()
oo.close()
os.remove('doc/ivs.timeseries.freqanalyse-module_.html')


files = sorted(glob.glob('doc/*module.html'))
files+= sorted(glob.glob('doc/*class.html'))
for myfile in files:
    shutil.move(myfile,myfile+'_')
    ff = open(myfile+'_','r')
    oo = open(myfile,'w')
    
    line_break = False
    for line in ff.readlines():
        
        if ']include figure]' in line or line_break:
            filename = line.split(']')[-2].strip()
            oo.write(r"<center><img src='images/%s' alt='[image example]' width='75%%'/></center>"%(filename)+'\n\n')
            print 'Added image %s to %s'%(filename,myfile)
            oo.write('\n\n')
            line_break = False
        
            
        elif ']include' in line:
            line_break = True
            
            
        else:
            oo.write(line)
        
    
    ff.close()
    oo.close()
    os.remove(myfile+'_')
    
