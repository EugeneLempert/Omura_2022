Greetings curious one. This repository is for some of the work done in the following publication:
Omura M, Takabatake Y, Lempert E, Benjamin-Hong S, D'Hulst C, Feinstein P. A genetic platform for functionally profiling odorant receptors in olfactory cilia ex vivo. Sci Signal. 2022 Aug 9;15(746):eabm6112. doi: 10.1126/scisignal.abm6112. Epub 2022 Aug 9. PMID: 35944068.

The specific work is the RT-qPCR results in Supplemental Figures S1 and S5B.

The Packages_Functions.R file contains the necessary libraries and functions to run the analysis scripts. Run it first.

The Figure_S1.R file will import the raw qPCR and melt curve files from a relative location, which requires the files themselves and the modification of the file location. This file expects the raw data to be located in an "Input" folder located at the same level as the working directory for the current R session. Similarly, the AA10 file is expected at the same level and should be imported before the rest of the code.

The Figure_S5B.R file has the same data source and thus the import step can be skipped if the Input_list and AA10 files already exist.

Over time, I might clean up some of the code, but perhaps it is best left in its current state as a testiment to my progress.

On 11.30.22, I will not be providing any of the raw data or the AA10 file without direct communication. I do not know the licensing agreements of the publisher as of this time, but I will update as I learn more.  

The information surrounding these analysis files, how they reached their current state, and my rationale up to that point- I can make available through direct communication.

Feel free to contact me at genielamp16@gmail.com to comment, question, or criticize my work.
