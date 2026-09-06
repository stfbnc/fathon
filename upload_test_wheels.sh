#!/bin/bash

echo "Uploading on Test Pypi"
python3 -m twine upload --verbose --skip-existing wheelhouse/* -r testpypi -u __token__ -p "$TEST_PYPI_TOKEN"

exit 0
