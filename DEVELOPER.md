# Developer Documentation

This document is intended to provide guidelines for software developers either maintaining or 
contributing to svtools development.

## Installing svtools

1. from a repo
1. from a conda thing
1. using pip package on pypi
1. downloading a github tarball
 
check out svtools
 <pre><code>git clone https://github.com/hall-lab/svtools.git svtools</code></pre>
pip install in the repo
 cd svtools
 pip install .
profit

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
