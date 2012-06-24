Reading .nex files
------------------

Use readNexFile function to read the contents of *.nex file in Matlab.
For example:

nexFileData = readNexFile(filePath);
nexFileData1 = readNexFile(); % if file path is not specified, Open File dialog will be used

nexFileData is a structure that contains all the data from .nex file.
See comments in the readNexFile.m file for more info on the contents of nexFileData sctucture.



Writing .nex files
------------------

Use writeNexFile function to write the contents of nexFileData structure to .nex file
For example:

writeResult = writeNexFile(nexFileData, filePath);
% modify nexFileData here...
writeResult1 = writeNexFile(nexFileData); % if file path is not specified, Save File dialog will be used

nexFileData is a structure that contains all the data from .nex file.
See comments in the writeNexFile.m file for more info on the contents of nexFileData sctucture.

The easieast way to save data in .nex format is to load nexFileData from a .nex file via readNexFile(),
then modify nexFileData and save it using writeNexFile.


