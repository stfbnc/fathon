#!/bin/bash

upload_ret=1
if [ "$PYPI_UPLOAD" = "test" ]; then
	#upload_ret=$(
	python3 -m twine upload --skip-existing wheelhouse/* -r testpypi -u "$TEST_PYPI_USR" -p "$TEST_PYPI_PWD" #)
elif [ "$PYPI_UPLOAD" = "release" ]; then
	#upload_ret=$(
	python3 -m twine upload --skip-existing wheelhouse/* -u "$PYPI_USR" -p "$PYPI_PWD" #)
fi

#if [ "$upload_ret" = 0 ]; then
#	echo "Succesfully uploaded wheels"
#else
#	echo "Cannot upload wheels, already existing"
#fi

#exit 0

