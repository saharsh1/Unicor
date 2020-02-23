UNICOR
=========================



1. ###### Download or clone the UNICOR package from GitHub. 

   Currently only the MATLAB version of the coded is available, and located in the `./matlab` subdirectory.  The contents of this folder are: 

   - engine (folder)
   - Example (folder)
   - Unicor2execution.m (matlab m file)
   - UnicorConfigDefault.ini (a default UNICOR initialization file)

   

   The `51Peg_Example` folder contains an example case of 51Peg:
   	51Peg.mat - a mat file containing the observed spectra 
   	UnicorConfig.ini - an ini file that defines the run configuration of UNICOR for the 51Peg example.

   The `HD95363_Example` folder contains a sample case of HD95363:
   	HD95363.mat - a mat file containing the observed spectra (measured by the eShel spectrograph)
   	UnicorConfig.ini - an ini file that defines the run configuration of UNICOR for the HD95363 example.

   

2. ###### Basic usage instructions:

   A file named `UnicorConfig.ini` should be always present in the data directory. This file contains the necessary parameters for UNICOR RV extraction for the spectra given in the mat file. If no config file is found, the default config file will be used.  Execute the code by typing the command:

   ```matlab
   [time,vels,evels,CCF] = Unicor2execution(data_path, mat file name>)
   ```

   

The code will generate the following variables:

- time - observation time, as were saved in the mat file.
   - vels - radial velocities, in km/s.
   - evels - errors on the radial velocities.
   - CCF - an array of structures containing the CCF  data for each of the observations.
   
   
   
   For example, in order to analyze the example spectra of HD95363, type:
   
   ```matlab
   [jd,vels,evels,CCF] = Unicor2execution('HD95363_Example', 'HD95363')
   ```
   
   **Note**: *If you run Unicor2execution() without any parameters, a file dialog will prompt you for the mat file.*



3. ###### Input data structure:

   The mat file contains a cell array. Each cell contains a structure that holds that data for a single observation, in the following fields : 

   - ​	name - a string with the stars name.

   - ​	bcv - baricentric velocity ( inkm/s) of the spectrograph relative to the solar system barycenter.

   - ​	bjd - barycentric Julian date of the observation. 

   - ​	wv - a matrix of wavelengths at each pixel of each oreder (number of spectrum points x norders). 

   - ​	sp - a matrix with the same dimensions of wv with the flux values.

     **Note:** *wv and sp must have the same number of rows (np) in all the orders and observations.*
     *This can be achieved either by trimming the  wv and sp vectors to the same common length, or by padding the wv and sp with NaN's to a common number of points.*
     

------

