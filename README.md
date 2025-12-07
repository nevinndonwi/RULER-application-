# RULER-application-

RULER - A component measuring approach to the sperm measuring problem (make measuring faster and easier than the currently manual used methods today). The final computing for biomedicine project of Nevin Ndonwi.
* I created a GUI that helps users in measuring sperm cells in a given sample image. 
* The application reports the length of the sperm that is specified by the user.
  
<b>Application setup and Initialization: </b>

* There are a couple of ways you can run this application. The already built executables are present at this link in Google drive: https://drive.google.com/drive/folders/16GqBbQO1myYrcalqOrmq4C-Ca_wnQcrR?usp=sharing 

* If you are able to use Linux then download the Linux version, it has the most up to date UI dependencies and it is able to render the zoom in and zoom out, move up, down, left, and right,  functionalities for more ease of use and user interaction. In addition, it will also enable you to be able to save the actual images of the measurments and sperm overlay windows in order to save your work and see the results if you ever want to look at them again. 

* If you have to use the windows exe file, then it will still work for being able to do the measuring functionalities, Due to limitations in executable file generation for windows exes and the fact I used the cv2 dependency for my UI, it will work but it will not be as robust as the linux version. 

* If you prefer to clone the repository, then you can look at requirements.txt to see all of the dependencies required for this application and attempt to build and run them all yourself as well. The only one that may be annoying to figure out is cv2 but everything else should be runnable unless you run into additional issues. You can try
* <code> pip install -r dependencies.txt </code>

and then run the script with:

*  <code>python3 FINAL_BIO_PROJECT.py</code>

<b>Application Usage: </b>

* Once you start the application, the gui containing your file system will appear and you can select an image of the sperm that you are trying to analyze.
*  From there, there are two windows that appear: One that contains the instructions and keys for different operations within the application and another containing the image you are analyzing. Following those instructions for the operations you want to do will help in your measurements. 

<b>Current Features: </b>

- You can load an image (jpeg, png, etc, anything supported by the application)

<b>While in auto-measure mode:<b>
- measure the different components of the sperm.
- invert the image to test if that improves the algorithm's ability to measure or to make hand tracing easier
- select and adjust which found components are included in the measurement.
- draw lines or pixels in order to fill gaps or holes or connect sperm parts that the algorithm did not catch.
- move around within the image (up, down, left, right) and zoom in/out
- Save images of the measurement results and/or the sperm highlighting.
- Return the length of the sperm in pixels and in millimeters. 

<b>While in manual-threshold adjustment/ component identification mode:<b>
- Select different components and designate them as different colors in order to get the measurement of different components to count as the same sperm or different sperm. 
- measure the different components of the sperm.
- invert the image to test if that improves component identificatiom
- select and adjust which found components are included in the measurement.
- draw lines or pixels in order to fill gaps or holes or connect sperm parts that the algorithm did not catch.
- move around within the image (up, down, left, right) and zoom in/out
- Save images of the measurement results and/or the sperm highlighting.
- Return the length of the sperm in millimeters.
- Save the summary to csv or save the images for future consideration.


Have fun in your project exploration!




