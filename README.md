# RULER-application-
RULER - A component measuring approach to the sperm measuring problem (make measuring faster and easier than the currently used methods today). The final computing for biomedicine project of Nevin Ndonwi.

<b>Application setup and Initialization: </b>

* There are a couple of ways you can run this application. The already built executables are present at this link in Google drive: https://drive.google.com/drive/folders/16GqBbQO1myYrcalqOrmq4C-Ca_wnQcrR?usp=sharing 

* If you are able to use Linux then download the Linux version, it has the most up to date UI dependencies and it is able to render the zoom in and zoom out, move up, down, left, and right,  functionalities for more ease of use and user interaction. In addition, it will also enable you to be able to save the actual images of the measurments and sperm overlay windows in order to save your work and see the results if you ever want to look at them again. 

* If you have to use the windows exe file, then it will still work for being able to do the measuring functionalities, Due to limitations in executable file generation for windows exes and the fact I used the cv2 dependency for my UI, it will work but it will not be as robust as the linux version. 

* If you prefer to clone the repository, then you can look at requirements.txt to see all of the dependencies required for this application and attempt to build and run them all yourself as well. The only one that may be annoying to figure out is cv2 but everything else should be runnable unless you run into additional issues. You can try pip install -r dependencies.txt and then run the script with:  <code>python3 FINAL_BIO_PROJECT.py</code>

<b>Application Usage: </b>
Once you start the application, the gui containing your file system will appear and you can select an image of the sperm that you are trying to analyze. From there, there are two windows that appear one that contains the instructions and keys for different operations within the application, and another containing the image you are analyzing. following those instructions for the operations you want to do will help in your measurements. 



