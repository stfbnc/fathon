#!/bin/bash

upload_ret=1
if [ "$PYPI_UPLOAD" = "test" ]; then
	echo "Uploading on test pypi"
	upload_ret=$(python3 -m twine upload --skip-existing wheelhouse/* -r testpypi -u "$TEST_PYPI_USR" -p "$TEST_PYPI_PWD")
elif [ "$PYPI_UPLOAD" = "release" ]; then
	echo "Uploading on pypi"
	upload_ret=$(python3 -m twine upload --skip-existing wheelhouse/* -u "$PYPI_USR" -p "$PYPI_PWD")
fi

if [ "$upload_ret" = 0 ]; then
	echo "Failed to upload wheels"
else
	echo "Successfully uploaded wheels"
fi

exit 0

