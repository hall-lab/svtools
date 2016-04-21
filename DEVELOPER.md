# Developer Documentation

This document is intended to provide guidelines for software developers either maintaining or 
contributing to svtools development.

## Installing svtools

1. Using pip package from pypi
1. From the git repo
1. From a conda package
1. Downloading a github tarball

### Using pip package from pypi
First you will need to prepare your python environment.
You might want to use pyenv virtualenv.
You can learn about that here <https://github.com/yyuu/pyenv-virtualenv>
We require python 2.7.X
The creation of the pyenv and activation looks like this in our environment
<pre><code>pyenv virtualenv 2.7.9 svtools_install_instructions-2.7.9
pyenv activate svtools_install_instructions-2.7.9</pre></code>
Now you will need to satisfy the pysam dependency
<pre><code>pip install pysam>=0.8.1,<0.9.0</pre></code>
Then you should be able to install the svtools package from pypi
<pre><code>pip install svtools</pre></code>
You can spot check your svtools install by running
<pre><code>svtools --version</pre></code>

<h3> From the git repo</h3>
First you will need to prepare your python environment.
You might want to use pyenv virtualenv. 
You can learn about that here <a href="https://github.com/yyuu/pyenv-virtualenv">https://github.com/yyuu/pyenv-virtualenv</a> 
We require python 2.7.X 
The creation of the pyenv and activation looks like this in our environment
<pre><code>pyenv virtualenv 2.7.9 svtools_install_from_repo-2.7.9
pyenv activate svtools_install_from_repo-2.7.9</pre></code> 
Once you have your python environment set up you will want to check out svtools from the hall-lab github repository.
To be on the bleeding edge you can install from master
<pre><code>git clone https://github.com/hall-lab/svtools.git svtools_test</code></pre>
<pre><code>cd svtools_test</pre></code>
Or you can discover the release tags on the repo and checkout the latest version
<pre><code>git tag -l</code></pre>
when you discover the version tag you wish to install you can switch to that version using
<pre><code>git checkout tags/v0.2.0b1</pre></code>
note: you can ignore the warning about "You are in 'detached HEAD' state."
OR you can just proceed to install from master.
Now install the dependencies suggested in the requiremnts files
<pre><code>pip install nose
pip install coverage
pip install statsmodels</pre></code>
Installing statsmodel can take a few minutes, but it satisfies the requirement for numpy, pandas, and scipy.
In our environment we need to specify pysam versions greater than 0.8.1 and less than 0.9.0
<pre><code>pip install 'pysam>=0.8.1,<0.9.0'</pre></code>
Now we can use pip to install svtools from within the repo
If you are not already in the directory
<pre><code>cd svtools_test</pre></code>
or just
<pre><code>pip install .</pre></code>
Finally we can spot check our svtools installation and observe the version number.
<pre><code>svtools --version</pre></code>

### From a conda package TODO
### Downloading a github tarball TODO
First you will need to prepare your python environment.
You might want to use pyenv virtualenv.
You can learn about that here <a href="https://github.com/yyuu/pyenv-virtualenv">https://github.com/yyuu/pyenv-virtualenv</a>
We require python 2.7.X
The creation of the pyenv and activation looks like this in our environment
<pre><code>pyenv virtualenv 2.7.9 svtools_install_instructions-2.7.9
pyenv activate svtools_install_instructions-2.7.9</pre></code>
Visit the svtools github page.
<https://github.com/hall-lab/svtools/>
Use the <pre>Download Zip</pre> button about half way down the page on the right hand side.  
Navigate to the download location on your filesystem and use unzip to expand the archive.
<pre><code>unzip svtools-master.zip</pre></code>
Now enter the directory that has been created.
<pre><code>cd svtools-master</pre></code>
To satisfy the pysam dependency of svtools you will want to install a version greater than 0.8.1 and less than 0.9.0
<pre><code>pip install 'pysam>=0.8.1,<0.9.0'</pre></code>
  
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
