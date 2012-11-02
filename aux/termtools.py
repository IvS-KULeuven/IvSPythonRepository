"""
Tools for cursor and color control in the terminal

Example usage (this does not work in the documentation):

>>> print green('this text is green')
this text is green
>>> print blink('this text blinks')+' this text does not blink'
this text blinks this text does not blink

You can combine all the possibilities, as the functions are generated on the fly:

>>> print blink_green_bgred_bold('blinking green bold text on red background')
blinking green bold text on red background
    
"""
from __future__ import print_function
import functools
import inspect
import sys
import types
import subprocess
import time
RED="\[\033[0;35m\]"
YELLOW="\[\033[0;33m\]"
instructs = {'black':"\033[30m",
          'red':"\033[31m",
          'green':"\033[32m",
          'yellow':"\033[33m",
          'blue':"\033[34m",
          'magenta':"\033[35m",
          'cyan':"\033[36m",
          'white':"\033[37m",
          'bgblack':"\033[40m",
          'bgred':"\033[41m",
          'bggreen':"\033[42m",
          'bgyellow':"\033[43m",
          'bgblue':"\033[44m",
          'bgmagenta':"\033[45m",
          'bgcyan':"\033[46m",
          'bgwhite':"\033[47m",
          'blink':'\033[5m',
          'underline':'\033[4m',
          'bold':'\033[1m',
          'reset':'\033[m'}

def overwrite_line(message):
    """
    Save cursor at current position, clear current line, print the message and reset the cursor.
    
    @param message: message to print to the screen
    @type message: str
    """
    ESC=chr(27)
    print('{ESC}[s{ESC}[2K{message}{ESC}[u'.format(ESC=ESC,message=message),end='')
    

def line_at_a_time(fileobj):
    """
    Return one line at a time from a file-like object.
    
    >>> #p1 = subprocess.Popen('ls',shell=True,stdout=subprocess.PIPE)
    >>> #for line in line_at_a_time(p1.stdout):
    ... #    print 'next line:',line.strip()
    >>> #retcode = p1.wait()
    
    use C{os.kill(p1.pid,SIGKILL)} to with C{SIGKILL} from C{signal} standard
    module to kill the process.
    
    return model_number
    Works around the iter behavior of pipe files in
    Python 2.x, e.g., instead of "for line in file" you can
    write "for line in line_at_a_time(file)"
    """
    while True:
        line = fileobj.readline()
        if not line:
            return
        yield line    


def subprocess_timeout(command, time_out):
    """
    
    Kill a running subprocess after a certain amount of time.
    
    Command represents the command for the process you would give in a terminal e.g. 'ls -l' in a terminal becomes ["ls", "-l"] or 'firefox' becomes ["firefox"]'; time_out is expressed in seconds. If the process did not complete before time_out, the process is killed.
    
    @param command: command to run
    @type command: str
    """

    # launching the command
    c = subprocess.Popen(command)

    # now waiting for the command to complete
    t = 0
    while t < time_out and c.poll() is None:
        time.sleep(1)  # (comment 1)
        t += 1

    # there are two possibilities for the while to have stopped:
    if c.poll() is None:
        # in the case the process did not complete, we kill it
        c.terminate()
        # and fill the return code with some error value
        returncode = -1  # (comment 2)

    else:                 
        # in the case the process completed normally
        returncode = c.poll()

    return returncode

class CallInstruct:
    """
    Generate a callable function on-the-fly.
    
    The function takes text as an input and will first print all the terminal
    instructions, then the text and then reset the terminal settings to the
    normal value.
    """
    def __init__(self,instr):
        """
        Remember which instructions to print.
        """
        self.instr = instr
    def __call__(self,text):
        """
        Print the instructions, the text and reset the terminal.
        """
        return "".join([instructs[i] for i in self.instr.split('_')]) + text + '\033[m'



class MemoryMonitor(object):

    def __init__(self, username=None):
        """Create new MemoryMonitor instance."""
        self.username = username
        self.time = []
        self.free = []
        self.used = []
    
    def usage(self):
        """Return int containing memory used by user's processes."""
        #self.process = subprocess.Popen("ps -u %s -o rss | awk '{sum+=$1} END {print sum}'" % self.username,
        #                                shell=True,
        #                               stdout=subprocess.PIPE,
        #                                )
        process = subprocess.Popen("free -mto",shell=True,stdout=subprocess.PIPE)
        self.time.append(time.time())
        used,free = process.communicate()[0].split('\n')[1].split()[2:4]
        self.free.append(free)
        self.used.append(used)
        



class Wrapper(object):
  """
  Wrap the module so that the C{getattr} function can be redefined.
  """
  def __init__(self, wrapped):
    self.wrapped = wrapped
  def __getattr__(self, name):
    # Perform custom logic here
    try:
      return getattr(self.wrapped, name)
    except AttributeError:
      return CallInstruct(name) # Some sensible default

#-- Wrap the module so that the C{getattr} function can be redefined.
sys.modules[__name__] = Wrapper(sys.modules[__name__])