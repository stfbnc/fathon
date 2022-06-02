Write-Host "Uploading on Test Pypi"
python -m twine upload --verbose --skip-existing wheelhouse\* --repository-url https://test.pypi.org/legacy/ -u $env:TEST_PYPI_USR -p $env:TEST_PYPI_PWD

Exit 0
