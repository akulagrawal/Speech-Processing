Simple Yes/No Detector
----------------------------------------------------------------------------------------------------------
## Steps to implement

- Read CoolEdit output file
- Correct DC Shift
- Normalize Data
- Crop out useful data
- Divide signals into frames and compute E and ZCR
- Decide Yes/No based on E and ZCR
----------------------------------------------------------------------------------------------------------
## Steps for Code Execution

- Open the project in Visual Studio 2010.
- In the main function, the initial parameters {keywords[], nSamples, dir} have to be set correctly.
- Save all the CoolEdit recorded text files in a single directory. Now, set the "dir" variable to the absolute path of this directory ending with "/".
  Currently, "Assignment1/Assignment1/Data" is used. If another directory has to be used, "dir" variable has to be updated.
- Currently, the code is written for testing total 20 files, 10 files named "yes_1.txt", "yes_2.txt", ..., "yes_10.txt", and other 10 files named "no_1.txt", "no_2.txt", ..., "no_10.txt".
- Hence, if the keywords and nSamples is not changed, these 20 files must be present for correct compilation of code.
- Else, keywords[] and nSamples can be changed accordingly. The files are expected to be in a format as described below:
	for keyword in keywords:
		there are <nSamples> files named "keyword_1.txt", "keyword_2.txt", ..., "keyword_<nSamples>.txt"
  Thus, if there is a single file, it can be renamed "file_1.txt". In this case, keywords = {"file"} and nSamples = 1.
- After setting the 3 variables, the code can be compiled and run from Visual Studio 2010 and he output will be obtained on the standard output console window.
- Segmented text files will be saved by appending "_seg" to the original file name, in the same directory.