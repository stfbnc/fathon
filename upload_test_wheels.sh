#!/bin/bash

echo "Uploading on Test Pypi"
python3 -m twine upload --verbose --skip-existing wheelhouse/* -r testpypi -u "$TEST_PYPI_USR" -p "$TEST_PYPI_PWD"

exit 0
