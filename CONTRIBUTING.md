# Contributing to fathon

If you have found bugs, or you want to propose a new feature:

1. describe the problem and specify the details of your setup (OS, Python version, environment, etc.);
2. provide a minimal working example that can be easily reproduced and tested;
3. open an issue with the previous information in the GitHub repository with tag <code>bug</code> or <code>proposed feature</code>.

PRs are also welcomed.

For any other question, open an issue with tag <code>question</code>.



## Documentation

Documentation is written in [reStructuredText](http://docutils.sourceforge.net/rst.html), built with <code>sphinx</code> and placed in the <code>docs</code> folder.

The script  <code>docs_gen.sh</code> generates `.py` files with only docstrings, since `sphinx` does not support Cython.

The variable `SPHINXBUILD` in the `Makefile` depends on the particular user's system and should be probably changed.

For any problem with the documentation, open an issue with tag <code>documentation</code>.