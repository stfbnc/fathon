#!/bin/bash

upload_ret=1
if [ "$PYPI_UPLOAD" = "test" ]; then
	echo "Uploading on testpypi"
	upload_ret=$(python3 -m twine upload --skip-existing wheelhouse/* -r testpypi -u "$TEST_PYPI_USR" -p "$TEST_PYPI_PWD")
elif [ "$PYPI_UPLOAD" = "release" ]; then
	echo "Uploading on pypi"
	upload_ret=$(python3 -m twine upload --skip-existing wheelhouse/* -u "$PYPI_USR" -p "$PYPI_PWD")
fi

exit "$upload_ret"

