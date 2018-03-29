Documentation
-------------

*   You can read online documentation at https://ivs-kuleuven.github.io/IvSPythonRepository/index.html

*   (Recommended) The repository also comes ready to produce local documentation, instructions can be found below.

Installation
------------

*   Clone the git repository to create a local copy as a folder called "ivs". Make sure the directory is re-named as "ivs".

    .. code-block:: bash

        $ cd python
        $ git clone https://github.com/IvS-KULeuven/IvSPythonRepository.git ivs

    This will clone all repository files into the ~/python/ivs folder. Be aware, however, that only the python scripts and the documentation are being cloned, not the (numerous and sometimes huge) datafiles that come along with it, containing, for example, limbdarkening coefficients. The directories of these files can be found, and specified, in ~ivs/config.py.

*   Updating your own clone of the IvS python repository to the most recent version can be done with:

    .. code-block:: bash

        $ cd ivs
        $ git pull

*   The IvS repository now comes with an Anaconda python environment file to help avoid user problems with the python modules. To create the repository anaconda environment run the following commands:

    .. code-block:: bash

        $ cd ivs
        $ conda env create -f IvS_repo.yml

    If you do not use anaconda you can find the repository dependencies in this file. Make sure to activate your environment each time before using the repository, either by sourcing it in your .bash_profile or running:

    .. code-block:: bash

        $ source activate ivs_repo

*   Make sure that your python path points to the directory above the ivs folder:

    If your python repository is in, e.g., ~/python/ivs/ , you can put in your .bash_profile:

    .. code-block:: bash

        export PYTHONPATH=/home/YOURNAME/python:$PYTHONPATH

    Warning: don't put ~/python/ivs in your Python path, but just ~/python.

    Once added to the PYTHONPATH you can import ivs modules using, for example:

    .. code-block:: python

        from ivs.statistics import linearregression

*   The IvS Python repository contains mostly python routines. Some of the time-critical functions, however, are written in fortran. To compile them you can run

    .. code-block:: bash

        $ python config.py compile

    If you want to specify your own fortran compiler you can do so with

    .. code-block:: bash

        $ python config.py compile f77

    Note: This requires the F2PY module which has previously caused problems... Sometimes the compilation process may fail. If so, try to compile spectra/pyrotin4.f manually, and then retry the automatic compilation (repeat for all troublesome files):

    .. code-block:: bash

        $ cd spectra/
        $ f2py --fcompiler=gfortran -c pyrotin4.f -m pyrotin4
        $ cd ../
        $ python config.py compile

*   In the config file you may also change the paths where the data catalogs (variable: data_dir) can be found, if you are not using the default locations (i.e. you are outside the institute).


*   To generate the documentation, simply run the script

    .. code-block:: bash

        $ python makedoc.py

    in the repository's root folder. This assumes that 'epydoc' is available which is
    already installed on all IvS computers. On your own laptop, you can get it from
    http://epydoc.sourceforge.net.

    Open "/doc/html/index.html" in your favorite browser and start browsing!
    Whenever you change something yourself in your local branch or you pull changes
    from someone else, you can re-run the makedoc.py script.


* Happy computing!

Encountered errors and their solutions:
---------------------------------------

:Q: When I run "python config.py compile", I get the following error:

    .. code-block:: bash

        numpy.distutils.fcompiler.CompilerNotFound: gnu95: f90 nor f77

:A: Install gfortran
