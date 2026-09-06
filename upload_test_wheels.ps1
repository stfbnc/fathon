Write-Host "Uploading on Test Pypi"
python -m twine upload --verbose --skip-existing wheelhouse\* --repository-url https://test.pypi.org/legacy/ -u __token__ -p $env:TEST_PYPI_TOKEN

Exit 0
