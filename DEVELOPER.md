# Developer Documentation

This document is intended to provide guidelines for software developers either maintaining or 
contributing to `svtools` development.

## Releasing a new version

### Tagging a release on github
`svtools` manages its versions using [python-versioneer](https://github.com/warner/python-versioneer). 
New versions are derived and generated from the names of tags in git. To release a new version, all 
that is needed is to tag the correct commit with an annotated tag. Always prepend versions with a 
'v' character.

For example, to release version 0.0.1 from the current commit:
```
git tag -a v0.0.1 -m 'v0.0.1'
git push --tags
```

Next navigate to the github [Releases page](https://github.com/hall-lab/svtools/releases), draft a 
release using the tag you just generated and add release information (a description of changes made since the last release). This will create an entry on the Github Releases page and upload the release to Zenodo.

### Build a pip package and upload to PyPI
Now that you have a new release, upload the package to [PyPI](https://pypi.python.org/pypi). To do so, you'll need a PyPI account, configuration for both the standard and test PyPI servers and necessary permissions on the `svtools` package. 

These instructions assume you have committed no additional changes after tagging the new release.

1. Build the new release and test by uploading to the PyPI test server.
  ```
  python setup.py sdist bdist_wheel upload -r pypitest
  ```
2. Verify that the package appears and information looks correct at https://testpypi.python.org/pypi
3. Build and upload the package to PyPI itself.
  ```
  python setup.py sdist bdist_wheel upload
  ```
4. In a fresh virtual environment, verify that the new package installs.
  ```
  pyenv virtualenv 2.7.9 test_new_package
  pyenv activate test_new_package
  pip install svtools
  ```
