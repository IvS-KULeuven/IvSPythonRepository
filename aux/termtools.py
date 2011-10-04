"""
Tools for cursor and color control in the terminal
"""
import functools
import inspect

def add_text(fctn):
    """
    Add possibility of text to each function.
    """
    @functools.wraps(fctn)
    def with_text(*args):
        #-- how many arguments does this function need?
        arguments = inspect.getargspec(fctn)[0]
        #-- if no extra arguments are given, just pass them
        if len(args)==len(arguments):
            return fctn(*args)
        #-- if an extra argument is given, assume it is text
        elif len(args)==len(arguments)+1:
            return fctn(*args[:-1])+ args[-1] + reset()
        else:
            raise ValueError('too many arguments')
    return with_text
        
colors = {'black':"\033[30m",
          'red':"\033[31m",
          'green':"\033[32m",
          'yellow':"\033[33m",
          'blue':"\033[34m",
          'magenta':"\033[35m",
          'cyan':"\033[36m",
          'white':"\033[37m"}
#BGGREEN=`echo "\033[42m"` # background to green 
#BGCYAN=`echo "\033[45m"` # background to cyan 
#BGWHITE=`echo "\033[47m"` # background to white 
#BGYELLOW=`echo "\033[43m"` # background to yellow 
#BGBLACK=`echo "\033[40m"` # background to black 
#BGBLUE=`echo "\033[44m"` # background to blue 
#BGRED=`echo "\033[41m"` # background to red 
#BGCYAN=`echo "\033[45m"` # background to cyan 
#BGMAGENTA=`echo "\033[45m"` # background to magenta 

#BGWHITEFGBLUE=`echo "\033[91m"` # background to white & foreground to blue 
##BGWHITEFGBLUE=`echo "\033[94m"` # background to white & foreground to blue 
#BGWHITEFGMAGENTA=`echo "\033[95m"` # background to white & fg to magenta 
#REVERSE=`echo "\033[7m"` #terminal to reverse video 


def reset():
    """
    Reset terminal to normal color.
    """
    return '\033[m'

@add_text
def blink():
    """
    Set terminal to blink
    """
    return '\033[5m'

@add_text
def underline():
    """
    Set terminal to underline
    """
    return '\033[4m'

@add_text
def nounderline():
    """
    Set terminal not to underline
    """
    return '\033[24m'

@add_text
def bold():
    """
    Set terminal text to bold.
    """
    return '\033[1m'

@add_text
def color(clr):
    return colors[clr]


if __name__=="__main__":
    print color('green')
    print 'this should be green'
    print reset()
    print 'this should be normal'
    print color('red','this should be red')
    print 'this should be normal again'
    for color_ in colors:
        print color(color_,'this should be %s'%(color_))
    print blink('this should be blinking')+' this should not be blinking'