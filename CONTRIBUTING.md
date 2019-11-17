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

##

## Code formatting

The current code is formatted using tabs.

In case you are going to modify *fathon*, please use the same docstring format already present in the source files.

## 

## Documentation

Contributing to SHARPy's documentation benefits everyone. As a developer, writing documentation helps you better understand what you have done and whether your functions etc make logical sense. As a user, any documentation is better than digging through the code. The more we have documented, the easier the code is to use and the more users we can have.

If you want to contribute by documenting code, you have come to the right place.

SHARPy is documented using Sphinx and it extracts the documentation directly from the source code. It is then sorted into directories automatically and a human readable website generated. The amount of work you need to do is minimal. That said, the recipe for a successfully documented class, function, module is the following:

1. Your documentation has to be written in ReStructuredText (rst). I know, another language... hence I will leave a few tips:

   - Inline code is written using two backticks ````

   - Inline math is written as `:math:`1+\exp^{i\pi} = 0` `. Don't forget the backticks!

   - Math in a single or multiple lines is simple:

     ```
         .. math:: 1 + \exp{i\pi} = 0
     ```

   - Lists in ReStructuredText are tricky, I must admit. Therefore, I will link to some [examples](http://docutils.sourceforge.net/docs/user/rst/quickref.html#enumerated-lists). The key resides in not forgetting the spaces, in particular when you go onto a second line!

   - The definite example list can be found [here](http://docutils.sourceforge.net/docs/user/rst/quickref.html).

2. Titles and docstrings, the bare minimum:

   - Start docstrings with `r` such that they are interpreted raw:

     ```
     r"""
     My docstring
     """
     ```

   - All functions, modules and classes should be given a title that goes in the first line of the docstring

   - If you are writing a whole package with an `__init__.py` file, even if it's empty, give it a human readable docstring. This will then be imported into the documentation

   - For modules with several functions, the module docstring has to be at the very top of the file, prior to the `import` statements.

3. We use the Google documentation style. See [description](https://google.github.io/styleguide/pyguide.html#38-comments-and-docstrings).

4. Function arguments and returns:

   - Function arguments are simple to describe: 

     ```
     def func(arg1, arg2):
     """Summary line.
     
     Extended description of function.
     
     Args:
       arg1 (int): Description of arg1
       arg2 (str): Description of arg2
     
     Returns:
       bool: Description of return value
     
     """
     return True
     ```

5. Solver settings:

   - If your code has a settings dictionary, with defaults and types then make sure that:

     - They are defined as class variables and not instance attributes.

     - Define a `settings_types`, `settings_default` and `settings_description` dictionaries.

     - After all your settings, update the docstring with the automatically generated settings table. You will need to import the `sharpy.utils.settings` module

       ```
       settings_types = dict()
       settings_default = dict()
       settings_description = dict()
       
       # keep adding settings
       
       settings_table = sharpy.utils.settings.SettingsTable()
       __doc__ += settings_table.generate(settings_types, settings_default ,settings_description)
       ```

6. See how your docs looks like!

   - Once you are done, run the following `SHARPy` command:

   ```
   sharpy any_string -d
   ```

   - If you are making minor updates to docstrings (i.e. you are not documenting a previously undocumented function/class/module) you can simply change directory to  `sharpy/docs` and run

   ```
   make html
   ```

   - Your documentation will compile and warnings will appear etc. You can check the result by opening

   ```
   docs/build/index.html
   ```

   and navigating to your recently created page.

   - Make sure that **before committing** any changes in the documentation you update the entire `docs` directory by running

   ```
   sharpy any_string -d
   ```

Thank you for reading through this and contributing to make SHARPy a better documented, more user friendly code!