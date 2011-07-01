# -*- coding: utf-8 -*-
"""
Convert a XML page to a dictionary
"""

import xml.parsers.expat
import urllib

class XMLParser:
    """
    Convert a XML page to a dictionary
    """
    def __init__(self,text,logger=None):
        self.text = text
        self.content = {}
        self.elem = []
        self.logger = logger
        p = xml.parsers.expat.ParserCreate()
        p.StartElementHandler = self.start_element
        p.EndElementHandler = self.end_element
        p.CharacterDataHandler = self.char_data
        if logger:
            logger.debug("XMLParser initialised")
        p.Parse(text)
        
    def start_element(self,name,attrs):
        """
        If a new element is found, at it to the dictionary
        """
        if len(self.elem)==0:
            self.content[name] = {}
            self.elem.append(name)
        else:
            curdict = self.content
            for ie in self.elem:
                curdict = curdict[ie]
            if not name in curdict.keys():
                curdict[name] = {}
            if not self.elem[-1]==name:
                self.elem.append(name)
            if attrs:
                key = attrs[attrs.keys()[0]]
                #print "START2",self.elem
                #print "START3",name,attrs,key,curdict[name].keys()
                if key not in curdict[name].keys():
                    curdict[name][key] = {}
                    self.elem.append(key)
    
    def end_element(self,name):
        """
        Remove the element from the queu (and everything after it)
        
        @parameter name: designation of the element
        @type name: string
        """
        #-- search for last occurence of 'name', and delete everything after it
        if len(self.elem)>1:
            index = len(self.elem) - 1 - self.elem[::-1].index(name)
            self.elem = self.elem[:index]
    
    def char_data(self,data):
        """
        Add the value of an element to the dictionary with its designation.
        
        @parameter data: value of the element
        @type data: string
        """
        #-- go to the final dictionary (most nested)
        curdict = self.content
        for ie in self.elem[:-1]:
            curdict = curdict[ie]
        #-- try to make it a float, if possible
        try:
            data = float(data)
        except:
            pass
        #-- if there is already a value for this element, make it a list and
        #   and the value, or just add the value.
        if curdict[self.elem[-1]] and not isinstance(curdict[self.elem[-1]],list):
            curdict[self.elem[-1]] = [curdict[self.elem[-1]]]
            curdict[self.elem[-1]].append(data)
        elif curdict[self.elem[-1]]:
            curdict[self.elem[-1]].append(data)
        else:
            curdict[self.elem[-1]] = data
        if self.logger:
            self.logger.debug("... %s: %s"%(self.elem,data))
        