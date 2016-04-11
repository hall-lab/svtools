# Developer Documentation

This document is intended to provide guidelines for software developers either maintaining or 
contributing to svtools development.

## Installing svtools

1. from a repo
1. from a conda thing
1. using pip package on pypi
1. downloading a github tarball

# From the git repo
prepare your python environment
  you might want to use pyenv virtual-env, you can learn about that here https://github.com/yyuu/pyenv-virtualenv 
  we require python 2.7.X 
  this process might look like this
  <pre><code>pyenv virtualenv 2.7.9 svtools_install_instructions-2.7.9
  pyenv activate</pre></code> 
check out svtools from the hall-lab github repository
 <pre><code>git clone https://github.com/hall-lab/svtools.git svtools_test</code></pre>
install dependencies as suggested in the requiremnts files
<pre><code>pip install nose
pip install coverage
pip install statsmodels</pre></code>
installing statsmodel can take a while, but in our testing it satisfies the requirement for numpy,pandas, and scipy
in our local environment we need to specify pysam versions less than 0.9.0, YMMV
<pre><code>pip install 'pysam>=0.8.1,<0.9.0'</pre></code>

use pip to install svtools from within the repo
<pre><code>cd svtools
pip install .</pre></code>
test your svtools installation
<pre><code>svtools --version</pre></code>

# From a conda package TODO
# Using pip package on pypi TODO
# Downloading a github tarball TODO

## Releasing a new version

svtools manages its versions using [python-versioneer](https://github.com/warner/python-versioneer). 
New versions are derived and generated from the names of tags in git. To release a new version, all 
that is needed is to tag the correct commit with an annotated tag. We do not prepend versions with a 
'v' character.

For example, to release version 0.0.1 from the current commit:
```
git tag -a 0.0.1 -m 'v0.0.1'
git push --tags
```

Next navigate to the github [Releases page](https://github.com/hall-lab/svtools/releases), draft a 
release using the tag you just generated and add release information.
