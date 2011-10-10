"""
Read or download files from the internet.
"""
import urllib
from ivs.aux import decorators

@decorators.retry_http(3)
def download(link,filename=None):
    """
    Download a file from a link.
    
    If you want to download the contents to a file, supply C{filename}. The
    function will return that filename as a check.
    
    If you want to read the contents immediately from the url, just give the link,
    and a fileobject B{and} the url object will be returned. Remember
    to close the url after finishing reading!
    
    @parameter link: the url of the file
    @type link: string
    @parameter filename: the name of the file to write to (optional)
    @type filename: str
    @return: output filename(, url object)
    @rtype: string(, FancyURLopener)
    """
    url = urllib.FancyURLopener()
    myfile,msg = url.retrieve(link,filename=filename)
    if filename is not None:
        url.close()
        return myfile
    else:
        return myfile,url