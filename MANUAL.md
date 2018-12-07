# How to install Python 2.7 packages (including GDAL and PCraster) <img src="./graphs/header_logo.png" align="right" />
Author: Joanna Suliga

If you are here it means that you are likely working with GIS software too. Commonly used software ArcGIS installs its own Python packages under the default path like for example: C:\Python27\ArcGIS10.1. It's very important to keep ArcGIS Python packages intact or some of ArcGIS functions can get corrupted! If you want to add your own packages or update existing ones, then you may prefer to work with another directory like: C:\Python27\Lib\site-packages. 

This manual is designed to help you with learning how to install Python packages (including GDAL and PCraster) without corrupting any other software or Windows environments. There are of course other ways to
solve this problem but please kindly find below my guidelines. If you know how to improve the manual, please don't hesitate to contact me! Let's make this manual better together!

Some remarks:
Unfortunately, manual installation of Python packages requires knowing the exact dependencies between all desired (and possible additional) packages.
Some platforms like Anaconda are very useful in previewing/installing/updating Python packages and getting Spyder Python interpreter but if not managed carefully, will likely corrupt ArcGIS packages. Especially older versions of ArcGIS are very vulnerable because they can only work with the older and not supported versions of Python packages. In addition, ArcGIS tends to overwrite the original Python interpreter & shell IDLE with its own version (also called IDLE…) that ignores windows environments specified by User. 

 Step | Explanation
------|-------------
 1-4  | How to import Python packages
 5    | How to install GDAL 
 6    | How to install PCraster
 7    | How to install Visual Studio Code

# Step 1. Install Python and ArcGIS
---

1.	Install Python 2.7
2.	Install ArcGIS, QGIS… 
3.	Install [Microsoft Visual C++ Compiler for Python 2.7](https://www.microsoft.com/EN-US/DOWNLOAD/confirmation.aspx?id=44266) from MS  
4.	Adjust Windows environments: File explorer -> This PC (Right click) -> Properties -> Advanced system settings -> Environment Variables… -> System variables -> Path (click and add new directories as follow:
a.	C:\Python27\
b.	C:\Python27\Scripts
c.	C:\Python27\Lib\site-packages

# Step 2. Set the original IDLE as default
---
ArcGIS installs its own Python interpreter which has the exact same name (IDLE) as the original Python interpreter. Even more, ArcGIS-IDLE sets itself as the default interpreter for all .py files. ArcGIS-IDLE is only searching for packages inside C:\Python27\ArcGIS10.1 and ignores Windows environments settings. Copy-pasting new packages to original ArcGIS folder may corrupt ArcPy therefore it’s highly advised to set original IDLE as default interpreter again and run scripts using packages located in C:\Python27\Lib\site-packages.
For Windows 10

1.	Right- click any .py file and select -> Open with…
2.	Tick the box ‘Always use this app to open .py files’
3.	Select -> More apps
4.	Scroll down and select -> Look for another app on this PC
5.	Navigate to C:\Python27\Lib\idlelib and select ‘idle’ 
Advanced alternative solution: Start -> type "run" > type "regedit"

# Step 3. Create a test .py file and open in the original python interpreter
---

1.	Create a text file and save as .py. 
2.	Open with IDLE
3.	Check python paths by selecting File -> Path Browser. Paths fixed in Step 1 should be there! 
4.	Test importing packages by adding lines ex. ‘import os’ ‘import numpy’ ‘import pandas’ ‘import gdal’ to the script and press F5. At this step, importing may fail.

# Step 4. Importing pure Python packages with PIP
---

1.	Open cmd (Windows command prompt) as administrator
2.	We will use inbuilt Python package PIP to preview/install/update. Assume we want to install numpy: Type in cmd: pip install numpy
Test importing your packages with test.py
Complex packages may require additional packages to be-preinstalled. Command ‘pip install’ should do this for you.
3.	Errors? Let’s see what could have gone wrong…
a.	Run CMD as administrator
b.	Microsoft Visual C++ 9.0 is missing
c.	pip install wheel (see Step 5.2)
Test importing your package with test.py
4.	Here you can find more commands to be typed in cmd:
python --version
pip --version pip is the inbuilt python package
pip list returns a list of installed using PIP packages
pip show ‘nameOfThePackage’ returns the version of the package ‘nameOfThePackage’ and paths
pip search ‘nameOfThePackage’ returns available packages
pip install ‘nameOfThePackage’ pulls and installs the newest package ‘nameOfThePackage’
pip uninstall ‘nameOfThePackage’
pip install --upgrade ‘nameOfThePackage’
pip freeze returns a list of packages installed with pip

# Step 5. Install GDAL
---
GDAL isn’t a pure Python code therefore it requires getting extra binary library for Windows. To pull installation file from the server, firstly download a binary library called [wheel](https://www.lfd.uci.edu/~gohlke/pythonlibs/#gdal). Select a version compatible with your Windows and python bit version. To install GDAL follow those steps:

1.	Open cmd (Windows command prompt) as administrator
2.	For GDAL 2.2.4 for 64bit Windows and Python 2.7 type in cmd:
pip install C:\Users\NameOfTheUser\Downloads\GDAL-2.2.4-cp27-cp27m-win_amd64.whl
pip install gdal
3.	Errors? Let’s see what could have gone wrong…
a.	Run cmd as administrator
b.	Microsoft Visual C++ 9.0 is missing
c.	pip install wheel (see Step 4.2)

# Step 6. Install PCRaster
---
PCRaster allows opening and previewing maps in .asc and .tif format without opening QGIS/ArcGIS.

1.	Download and unzip [PCraster](http://pcraster.geo.uu.nl/downloads/latest-release/) in C:\Program Files
2.	Adjust Windows environments: File explorer -> This PC (Right click) -> Properties -> Advanced system settings -> Environment Variables… -> System variables 
-> Path (click and add C:\Program Files\pcraster-4.1.0_x86-64\bin)
-> PYTHONPATH (click and add C:\Program Files\pcraster-4.1.0_x86-64\python)
3.	To auto-open .asc .tif files with Aguila from PCraster repeat points 1-4 from the Step 4. Navigate to C:\Program Files\pcraster-4.1.0_x86-64\bin and select aguila.

# Step 7. Install Visual Studio Code
---
An alternative for Spyder interpreter might be Visual Studio Code.

1.	Download and install [Visual Studio Code](https://code.visualstudio.com/Download)
2.	Make sure that you have the following paths in your ‘Path’ Windows environmental variables
a.	C:\Python27\
b.	C:\Python27\Scripts
c.	C:\Python27\Lib\site-packages
d.	C:\Program Files\Microsoft VS Code\bin

(Special thanks for Celray!)

