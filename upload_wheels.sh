#!/bin/bash

upload_ret=$(twine upload --skip-existing wheelhouse/* -r testpypi -u "$TEST_PYPI_USR" -p "$TEST_PYPI_PWD")

if [ "$upload_ret" = 0 ]; then
	echo "Succesfully uploaded wheels"
else
	echo "Cannot upload wheels, already existing"
fi

exit 0

