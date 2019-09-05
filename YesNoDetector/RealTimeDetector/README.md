Real-Time Yes/No Detector
-----------------------------------------------------------------------------------------------
## Steps to implement

- Instead of taking CoolEdit output data, use Recording_Module.exe to record voice and store in text file.
-----------------------------------------------------------------------------------------------
## Steps for Code Execution

- Open the project in Visual Studio 2010.
- Copy the Recording_Module folder in a location whose address doesn't contain spaces.
- In the main function, the initial parameter exeCommand has to be set accordingly as described below.
- Usage of exeCommand:
	system("PATH TIME input_file.wav input_file.txt");
  PATH: Complete path to "Recording_Module.exe" file.
  TIME: Time(positive integer) in seconds
  input_file: Save recording as this name (.wav and .txt files)
  This will record for TIME seconds and amplitude values of samples will be stored in input_file.txt(no header).