# Contributing to fathon



## Bug fixes and new features

If you have found or fixed bugs, or you want to propose or add a new feature, please follow these instructions.

#### To report a bug:

1. Describe the problem accurately and specify the details of your setup (OS, Python version, environment, etc.);
2. Provide a minimal working example that can be easily reproduced and tested;
3. Open an issue with the previous information in the Github repository with tag <code>bug</code>.

#### To fix a bug:

1. Fork the repository;
2. If the bug involves Cython code and you know Cython:
   - Fix the bug;
   - Describe the bug, provide a minimal working example that can be easily reproduced and tested to reproduce the bug, and explain how you fixed it;
   - Create a pull request including what described in the previous point.
3. If the bug involves Cython code and you do not know Cython:
   - Simply report the bug as described in "To report a bug".
4. If the bug involves Python code:
   - Follow the same instructions given in point 2.

#### To propose a new feature:

1. Open an issue in the Github repository with tag  <code>proposed feature</code>;
2. Please try to give a detailed description of the feature, providing references if available.

#### To add a new feature:

1. Fork the repository;

2. If the new feature needs to be written in Cython and you know Cython:

   - Develop the new code;

   - Describe the new feature and provide references if available, provide a minimal working example that can be easily reproduced and tested, and list all the added or modified repository's files;

   - Create a pull request including what described in the previous point.

3. If the new feature needs to be written in Cython and you do not know Cython:

   - Develop the new code in Python;
   - Describe the new feature and provide references if available, provide the new code with information on where to include it, and provide a minimal working example that can be easily reproduced and tested;
   - Open an issue with the previous information in the Github repository with tag  <code>new feature</code>;
   - I will review it and translate the code to Cython.

4. If the new feature can be written in Python:

   - Follow the same instructions given in point 2.

#### For any other question regarding the package:

1. Open an issue as described in "To report a bug", but with tag <code>question</code>, or write to fathon.package@gmail.com.



## Code formatting

The current code is formatted using tabs.

In case you are going to modify `fathon`, please use the same docstring format already present in the source files.



## Documentation

Documentation is written in [reStructuredText](http://docutils.sourceforge.net/rst.html), built with <code>sphinx</code> and placed in the <code>docs</code> folder. If you have contributed, please update also the documentation.

- Add or modify the function/class in the corresponding <code>.rst</code> file in the folder <code>fun_class</code>, and follow the same naming convention if a new file is needed;
- In case of a new file, add it at the end of <code>index.rst</code>;
- Install <code>sphinx</code> and <code>numpydoc</code>;
  - Run <code>docs_gen.sh</code> inside the <code>docs</code> folder (the variable `SPHINXBUILD` in the `Makefile` should be probably changed), and check the result file <code>_build/html/index.html</code>.

For any problem with the documentation, open an issue with tag <code>documentation</code>.