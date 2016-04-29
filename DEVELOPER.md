This document is intended to provide guidelines for software developers either maintaining or 
contributing to `svtools` development.

## Installing `svtools`

1. From a [conda][1] package <sup>[1](#conda-install)</sup>
2. Using the [pip][2] package from [pypi][3] <sup>[2](#pip-install)</sup>
3. From the git repo <sup>[3](#git-install)</sup>
4. Downloading a github tarball <sup>[4](#tarball-install)</sup>

### <a name="conda-install"></a> From a conda package TODO 

### <a name="pip-install"></a> Using pip package from pypi

`svtools` is based on `python 2.7`.  First you will need to prepare your python environment.  You might want to use [pyenv virtualenv][4].
The creation of the pyenv virtual environment and activation looks like this:

    pyenv virtualenv 2.7.9 svtools_install_instructions-2.7.9
    pyenv activate svtools_install_instructions-2.7.9

Now you will need to satisfy the [`pysam`][10] dependency:

    pip install pysam>=0.8.1,<0.9.0

Then you should be able to install the `svtools` package from pypi:

    pip install svtools

You can spot check your `svtools` install by running:

    svtools --version

### <a name="git-install"></a> From the git repo

`svtools` is based on `python 2.7`.  First you will need to prepare your python environment.  You might want to use [pyenv virtualenv][4].
The creation of the pyenv virtual environment and activation looks like this:

    pyenv virtualenv 2.7.9 svtools_install_from_repo-2.7.9
    pyenv activate svtools_install_from_repo-2.7.9

Once you have your python environment set up you will want to check out `svtools` from the [hall-lab github repository][5].  The master branch contains the bleeding edge version.

    git clone https://github.com/hall-lab/svtools.git svtools_test
    cd svtools_test

Or you can discover the release tags on the repo and checkout the latest version

    git tag -l

when you discover the version tag you wish to install you can switch to that version using

    git checkout tags/v0.2.0b1

_note:_ you can ignore the warning about "You are in 'detached HEAD' state."

OR, you can just proceed to install from master.

Now install the dependencies suggested in the requiremnts files

    pip install nose
    pip install coverage
    pip install statsmodels

Installing `statsmodels` can take a few minutes, but it satisfies the requirement for `numpy`, `pandas`, and `scipy`.

In our environment we need to specify [`pysam`][10] versions greater than 0.8.1 and less than 0.9.0

    pip install 'pysam>=0.8.1,<0.9.0'

Now we can use `pip` to install `svtools` from within the repo. If you are not already in the directory

    cd svtools_test

or just

    pip install .

Finally we can spot check our `svtools` installation and observe the version number.

    svtools --version

### <a name="tarball-install"></a> Downloading a github tarball

`svtools` is based on `python 2.7`.  First you will need to prepare your python environment.  You might want to use [pyenv virtualenv][4].
The creation of the pyenv virtual environment and activation looks like this:

    pyenv virtualenv 2.7.9 svtools_install_instructions-2.7.9
    pyenv activate svtools_install_instructions-2.7.9

Visit the [svtools releases github page][6].  Select the latest release and use the `Source code (tar.gz)` link to download a tarball of the source.

Navigate to the download location on your filesystem and use the `tar` command:

    tar -xvzf svtools-0.2.0b1.tar.gz
    
to expand the archive.  [While you wait enjoy this cartoon from xkcd][7]. 

Now enter the directory that has been created.

    cd svtools-0.2.0b1

Now install the dependencies suggested in the requiremnts files:

    pip install nose
    pip install coverage
    pip install statsmodels

Installing `statsmodel` can take a few minutes, but it satisfies the requirement for `numpy`, `pandas`, and `scipy`.

In our environment we need to specify [`pysam`][10] versions greater than 0.8.1 and less than 0.9.0:

    pip install 'pysam>=0.8.1,<0.9.0'
    pip install .

Finally we can spot check our `svtools` installation and observe the version number with the following command:

    svtools --version

## Releasing a new version

`svtools` manages its versions using [python-versioneer][8].  New versions are derived and generated from the names of tags in git. To release a new version, all that is needed is to tag the correct commit with an annotated tag. We do not prepend versions with a 'v' character.

For example, to release version 0.0.1 from the current commit:

    git tag -a 0.0.1 -m 'v0.0.1'
    git push --tags

Next navigate to the github [Releases page][9], draft a release using the tag you just generated and add release information.

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

