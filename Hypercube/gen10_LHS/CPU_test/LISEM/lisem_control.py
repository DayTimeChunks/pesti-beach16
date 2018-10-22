import os
import shlex
import subprocess

print(os.getcwd())

path_run = r"..\\Alteck_test.run"
command_line = r"lisem -b -r " + str(path_run)
args = shlex.split(command_line)

print(args)

runs = 0
while runs < 2:
    runs += 1
    print("Running...")
    code = subprocess.call(args)  # Calling LISEM!

    if code == -1073741819:
        print("Code was good: ", str(code))
    else:
        print("Different code: ", str(code))
        print("Exiting test loop")
        break

print("Finished all runs!! Total runs: " + str(runs))




"""
Notes to self...

From LISEM documentation:

Get LISEM to run, with the '.run' as input file

cmd:
lisem -b -r Alteck_test.run

Notes:
lisem -b -ni -no -r runfilename.run -c options

-r REQUIRED Provide the runfile name (including full path and extension)
-b Optional Use LISEM in batch mode with interface
-ni Optional Use LISEM without interface
-no Optional Use LISEM with only error output
-c Optional Provide options, these will overwrite those from the run file. 

Options should have the following format [<option name>=<Value>;<option2 name>=Value2>;...] 
The option name should be an exact copy of the name provided within this document, which is 
identical to the name in the run file. For example: -c [SwitchChannel=0;SwitchNoErosion=1] 
This would turn off the usage of a channel and the usage of erosion.

All options, settings and map names can be set with the '-c' argument.

"""
