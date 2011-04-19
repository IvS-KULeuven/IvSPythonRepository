# -*- coding: utf-8 -*-
"""
A simple interface to work with a database saved on the hard disk. 

Author: Robin Lombaert

"""

import os
import cPickle
import time


class Database(dict):
    
    '''
    A database class.
    
    The class creates and manages a dictionary saved to the hard disk. 
    
    It functions as a python dictionary with the extra option of synchronizing
    the database instance with the dictionary saved on the hard disk. 
    
    No changes will be made to the hard disk copy, unless Database.sync() is 
    called.
    
    Note that changes made on a deeper level than the (key,value) pairs of the 
    Database (for instance in the case where value is a dict() type itself) 
    will not be automatically taken into account when calling the sync() 
    method. The key for which the value has been changed on a deeper level has 
    to be added to the Database.__changed list by calling addChangedKey(key)
    manually.
    
    Running the Database.sync() method will not read the database from the hard
    disk if no changes were made or if changes were made on a deeper level 
    only. In order to get the most recent version of the Database, without 
    having made any changes, use the .read() method. Note that if changes were 
    made on a deeper level, they will be lost.
    
    Example:
    
    >>> import os
    >>> from ivs.io import database
    >>> filename = 'mytest.db'
    >>> db = database.Database(filename)
    No database present at mytest.db. Creating a new one.
    >>> db['test'] = 1
    >>> db['test2'] = 'robin'
    >>> db.sync()
    >>> db2 = database.Database(filename)
    >>> print db2['test']
    1
    >>> print db2['test2'] 
    robin
    >>> db2['test'] = 2
    >>> db2.sync()
    >>> db.sync()
    >>> print db['test']
    1
    >>> db.read()
    >>> print db['test']
    2
    >>> del db2['test2']
    >>> db2.sync()
    >>> print db['test2']
    robin
    >>> db.read()
    >>> print db['test2']
    Traceback (most recent call last):
      File "<stdin>", line 1, in <module>
    KeyError: 'test2'
    >>> test_dict = dict()
    >>> db['test'] = test_dict
    >>> db.sync()
    >>> db2.read()
    >>> print db2['test']
    {}
    >>> db['test']['test'] = 1
    >>> db.sync()
    >>> db2.read()
    >>> print db2['test']
    {}
    >>> db.addChangedKey('test')
    >>> db.sync()
    >>> db2.read()
    >>> print db2['test']
    {'test': 1}
    >>> db.setdefault('test','defkey')
    {'test': 1}
    >>> db.setdefault('test3','defval')
    'defval'
    >>> db.sync()
    >>> db2.read()
    >>> print db2['test3']
    defval
    >>> os.system('rm %s'%filename)
    0
    '''
    
    
    def __init__(self,db_path):
        
        '''
        Initializing a Database class.
        
        Upon initialization, the class will read the dictionary saved at the 
        db_path given as a dictionary.
        
        Note that cPickle is used to write and read these dictionaries.
        
        If no database exists at db_path, a new dictionary will be created.
        
        @param db_path: The path to the database on the hard disk.
        @type db_path: string
  
        '''
        
        super(Database, self).__init__()
        self.db_path = db_path
        self.read()
        self.__changed = []
        self.__deleted = []
      
      
      
    def __delitem__(self,key):
        
        '''
        Delete a key from the database.
        
        This deletion is also done in the hard disk version of the database 
        when the sync() method is called. 
        
        This method can be called by using syntax:
        del db[key]
        
        @param key: a dict key that will be deleted from the Database in memory
        @type key: a type valid for a dict key
        
        '''
        
        self.__deleted.append(key)
        return super(Database,self).__delitem__(key)
        
    
    
    
    def __setitem__(self,key,value):
        
        '''
        Set a dict key with value. 
        
        This change is only added to the database saved on the hard disk when
        the sync() method is called. 
        
        The key is added to the Database.__changed list.
        
        This method can be called by using syntax:
        db[key] = value
        
        @param key: a dict key that will be added to the Database in memory
        @type key: a type valid for a dict key
        @param key: value of the key to be added
        @type value: any
        
        '''
        
        self.__changed.append(key)
        return super(Database,self).__setitem__(key,value)
        
        
        
    def setdefault(self,key,*args):
        
        '''
        Return key's value, if present. Otherwise add key with value default
        and return. 
        
        Database.__changed is updated with the key if it is not present yet.
        
        @param key: the key to be returned and/or added.
        @type key: any valid dict() key
        @param args: A default value added to the dict() if the key is not 
        present. If not specified, default defaults to None.
        @type args: any type
        @return: key's value or default
                
        '''
        
        if not self.has_key(key):
            self.__changed.append(key)
        return super(Database,self).setdefault(key,*args)
        
        
            
    def pop(self,key,*args):
        
        '''
        If database has key, remove it from the database and return it, else 
        return default. 
    
        If both default is not given and key is not in the database, a KeyError
        is raised. 
        
        If deletion is successful, this change is only added to the database 
        saved on the hard disk when the sync() method is called. 
        
        The key is added to the Database.__deleted list, if present originally.
                
        @param key: a dict key that will be removed from the Database in memory
        @type key: a type valid for a dict key
        @param args: value of the key to be returned if key not in Database
        @type args: any
        @return: value for key, or default
        
        '''
        
        if self.has_key(key):
            self.__deleted.append(key)
        return super(Database,self).pop(key,*args)
        
        
    def popitem(self):
        
        '''
        Remove and return an arbitrary (key, value) pair from the database.
        
        A KeyError is raised if the database has an empty dictionary.
        
        If removal is successful, this change is only added to the database 
        saved on the hard disk when the sync() method is called. 
        
        The removed key is added to the Database.__deleted list.
                
        @return: (key, value) pair from Database
        
        '''
        
        (key,value) = super(Database,self).popitem()
        self.__deleted.append(key)
        return (key,value)
            
            
        
    def update(self,*args,**kwargs):
        
        '''
        
        Update the database with new entries, as with a dictionary. 
        
        This update is not synched to the hard disk! Instead Database.__changed
        includes the changed keys so that the next sync will save these changes
        to the hard disk.
        
        @param args: A dictionary type object to update the Database.
        @type args: dict()
        @keyword kwargs: Any extra keywords are added as keys with their values.
        @type kwargs: any type that is allowed as a dict key type.
        
        '''
        
        self.__changed.extend(kwargs.keys())
        self.__changed.extend(args[0].keys())
        return super(Database,self).update(*args,**kwargs)
               
               
        
    def read(self):
        
        '''
        Read the database from the hard disk.
        
        Whenever called, the database in memory is updated with the version 
        saved on the hard disk.
        
        Any changes made outside the session of this Database() instance will
        be applied to the database in memory! 
        
        Any changes made to existing keys in current memory before calling 
        read() will be undone! Use sync() instead of read if you want to keep
        current changes inside the session. 
        
        If no database is present at the path given to Database() upon 
        initialisation, a new Database is made by saving an empty dict() at the
        requested location.        
        
        Reading and saving of the database is done by cPickle-ing the dict(). 
        
        '''
        
        try:
            dbfile = open(self.db_path,'r')
            while True:
                try:
                    db = cPickle.load(dbfile)
                    break
                except ValueError:
                    print 'Loading database failed: ValueError ~ insecure '+\
                          'string pickle. Waiting 10 seconds and trying again.' 
                    time.sleep(10)
            dbfile.close()
            self.clear()
            super(Database,self).update(db)
        except IOError:
            print 'No database present at %s. Creating a new one.'%self.db_path
            self.__save()
                
                
                
    def sync(self):
        
        ''' 
        
        Update the database on the harddisk and in the memory.
         
        The database is read anew, ie updated with the hard disk version to 
        account for any changes made by a different program. Next, the changes
        made to the database in memory are applied, before saving the database
        to the hard disk again.
        
        Any items deleted from the database in memory will also be deleted from
        the version saved on the hard disk!
        
        The keys that are changed explicitly are all listed in self.__changed,
        to which entries can be added manually using the addChangedKey method, 
        or automatically by calling .update(), .__setitem__() or .setdefault().
        
        '''
        
        if self.__changed or self.__deleted:
            current_db = dict([(k,v) 
                               for k,v in self.items() 
                               if k in set(self.__changed)])
            self.read()
            self.__deleted = list(set(self.__deleted))
            while self.__deleted:
                try:
                    super(Database,self).__delitem__(self.__deleted.pop())
                except KeyError:
                    pass
            super(Database,self).update(current_db)
            self.__save()
            self.__changed = []
    
    
    
    def __save(self):
        
        '''
        
        Save a database. 
        
        Only called by Database() internally. Use sync() to save the Database
        to the hard disk.
        
        Reading and saving of the database is done by cPickle-ing the dict(). 
        
        '''
        
        dbfile = open(self.db_path,'w')
        cPickle.dump(self,dbfile)
        dbfile.close()
    
    
    
    def addChangedKey(self,key):
        
        '''
        Add a key to the list of changed keys in the database.
        
        This is useful if a change was made to an entry on a deeper level, 
        meaning that the __set__() method of Database() is not called directly.
        
        If the key is not added to this list manually, it will not make it into
        the database on the hard disk when calling the sync() method.
        
        @param key: the key you want to include in the next sync() call.
        @type key: string
        
        '''
        
        self.__changed.append(key)
    
    
    
    def getDeletedKeys(self):
        
        '''
        Return a list of all keys that have been deleted from the database in
        memory.
        
        @return: list of keys
        @rtype: list
        '''
        
        return self.__deleted
    
    
    
    def getChangedKeys(self):
        
        '''
        Return a list of all keys that have been changed in the database in
        memory.
        
        @return: list of keys
        @rtype: list
        '''
        
        return self.__changed
        
if __name__ == "__main__":
    import doctest
    doctest.testmod()
    