if ($env:PYPI_UPLOAD -eq "test") {
    Write-Host "Uploading on Test Pypi"
    python -m twine upload --verbose --skip-existing wheelhouse\* --repository-url https://test.pypi.org/legacy/ -u $env:TEST_PYPI_USR -p $env:TEST_PYPI_PWD
} elseif ($env:PYPI_UPLOAD -eq "release") {
    Write-Host "Uploading on Pypi"
    python -m twine upload --verbose --skip-existing wheelhouse\* -u $env:PYPI_USR -p $env:PYPI_PWD
}

Exit 0
