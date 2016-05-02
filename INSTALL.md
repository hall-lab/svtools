## Installing `svtools`

### Table of Contents
1. [Preparing your Python environment](#python-env)
2. [Installing from a conda package](#conda-install)
3. [Installing from PyPI using pip](#pip-install)
4. [Installing directly from the git repo](#git-install)
5. [Installing from a github release tarball](#tarball-install)

### <a name="python-env"></a> Preparing your Python environment
`svtools` requires Python 2.7.  Before proceeding you need to prepare your Python environment.  We highly recommend that you manage your installation of `svtools` with the [pip][2] package manager as shown below.  pip is installed by pyenv by default.  Using pip will allow you to uninstall `svtools` easily.  You might want to use [pyenv virtualenv][4] to create a virtual environment. There are instructions for setting up a virtualenv on the [pyenv github site](https://github.com/yyuu/pyenv/blob/master/README.md).

The creation of the pyenv virtual environment and activation looks like this:

    pyenv virtualenv 2.7.9 svtools-2.7.9
    pyenv activate svtools-2.7.9

### <a name="conda-install"></a> Installing using a [conda][1] package
We will eventually have a conda package available for install, but this version doesn't currently have one. 

### <a name="pip-install"></a> Installing from [PyPI][3] using [pip][2]
Once you have your Python environment set up you should be able to install the `svtools` package from PyPI:

    pip install svtools

_note:_ On older systems you may encounter an [error during pysam installation](https://github.com/pysam-developers/pysam/issues/262). This can be solved by specifying a version of [pysam][10] greater than 0.8.1 and less than 0.9.0

    pip install 'pysam>=0.8.1,<0.9.0'

You can spot check your `svtools` install by running:

    svtools --version

### <a name="git-install"></a> Installing directly from the git repo
Once you have your Python environment set up you will want to clone `svtools` from the [hall-lab github repository][5].  

    git clone https://github.com/hall-lab/svtools.git svtools_test
    cd svtools_test

Or you can discover the release tags on the repo and checkout the latest version

    git tag -l

when you discover the version tag you wish to install you can switch to that version using

    git checkout tags/v0.2.0b1

_note:_ you can ignore the warning about "You are in 'detached HEAD' state."

OR, you can just proceed to install from master.

Now install the dependencies suggested in the requirements file

    pip install statsmodels

Installing statsmodels can take a few minutes, but it satisfies the requirement for numpy, pandas, and scipy.
_note:_ On older systems you may encounter an [error during pysam installation](https://github.com/pysam-developers/pysam/issues/262). This can be solved by specifying a version of [pysam][10] greater than 0.8.1 and less than 0.9.0 

    pip install 'pysam>=0.8.1,<0.9.0'

Now we can use pip to install `svtools` from within the repo. If you are not already in the directory

    cd svtools_test

or just

    pip install .

Finally we can spot check our `svtools` installation and observe the version number.

    svtools --version

### <a name="tarball-install"></a> Installing from a github release tarball

Once you have your python environment set up, visit the [svtools releases github page][6].  Select the latest release and use the `Source code (tar.gz)` link to download a tarball of the source.

Navigate to the download location on your filesystem and use the tar command:

    tar -xvzf svtools-0.2.0b1.tar.gz
    
to expand the archive.  [While you wait, enjoy this cartoon from xkcd][7]. 

Now enter the directory that has been created.

    cd svtools-0.2.0b1

Now install the dependencies suggested in the requiremnts files:

    pip install statsmodels

Installing statsmodel can take a few minutes, but it satisfies the requirement for numpy, pandas, and scipy.

In our environment we need to specify [pysam][10] versions greater than 0.8.1 and less than 0.9.0:

    pip install 'pysam>=0.8.1,<0.9.0'
    pip install .

Finally we can spot check our `svtools` installation and observe the version number with the following command:

    svtools --version

[1]: http://conda.pydata.org/docs/
[2]: https://pypi.python.org/pypi/pip/
[3]: https://pypi.python.org/pypi
[4]: https://github.com/yyuu/pyenv-virtualenv
[5]: https://github.com/hall-lab/svtools
[6]: https://github.com/hall-lab/svtools/releases
[7]: https://xkcd.com/1168/
[8]: https://github.com/warner/python-versioneer
[9]: https://github.com/hall-lab/svtools/releases
[10]: https://github.com/pysam-developers/pysam

