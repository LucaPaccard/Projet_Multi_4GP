# Projet_Multi_4GP
This documentation presents the Graphic User Interface (GUI) developed during the 4th year multidisciplinary project at the Physics Department of INSA Toulouse.  
Its aim is to study quantum dots directed assembly. It combines a quartz crystal microbalance, an injection system and a camera module.  
The collaborating laboratory of this project is the Laboratory of Physics and Chemistry of Nano-Objects (LPCNO).

## Name
openQCM Next GUI Software — Injection and Camera Implementations

## Programming language
Python

## Date 
2022-01-10

## Author
openQCM Team - Marco  
INSA Team - Luca PACCARD, Arthur LEMAIRE

## Description
This software application is developed in Python, which is an open source, object - oriented and suited for scientific application programming language. 
Python makes the software program easy to modify and develop for custom application.

It combines different parts: 

The openQCM NEXT part exploits all the main functionalities of the openQCM NEXT device. Real time monitoring of frequency and dissipation on the fundamental mode and overtone harmonics. It is possible to acquire almost simultaneously 5 sweep signals and elaborate the frequency and dissipation measurement in roughly 700 msec. In addition, the application allows to control and monitor the sensor module temperature in real time.

An injection module has been added with a modular pressure-based flow controller from Fluigent (https://www.fluigent.com/). The Flow EZ™ (https://www.fluigent.com/research/instruments/pressure-flow-controllers/lineup-series/flow-ez/) is the most advanced system available for pressure-based flow control. The compact device stands near the microfluidic device, allowing the user to minimize bench space use without the need of a PC. One can be operational and generate data rapidly. The Flow EZ™ supports reservoir sizes from 2 mL to one liter laboratory bottles. One can use large reservoirs and maintain continuous, pulseless flow for days without refilling.

A camera modude is present to monitor in real time the processes happening in the microbalance. it is capable of taking screenshots and capturing images for a given duration and framerate.

The application is based mainly on multiprocessing package (https://docs.python.org/3/library/multiprocessing.html).

## Intended Audience
Science/Research/Engineering

## Requirements
Requirements:
- Python 3.7 (verified compatibility with Python 3.6) (https://www.python.org/).
- Anaconda3-5.3.0 
     External Packages:
     - PyQt5 5.9.2 (https://pypi.org/project/PyQt5/).
     - PySerial 3.4 (https://pypi.org/project/pyserial/).
     - PyQtGraph 0.10.0 (http://www.pyqtgraph.org/).

Other internal packages used:
- multiprocessing, numpy, scipy, setuptools, io, platform, sys, enum, argparse, cvs, time, datetime, logging, progressbar, pandas etc.

## Installation instructions/guide:
### Source code
Download source code here:
https://github.com/LucaPaccard/Projet_Multi_4GP

Windows, macOS
  1.  Download and install Anaconda3 for Python 3.7 version Anaconda3-5.3.0  https://www.anaconda.com/download/
  2.  Open Anaconda3 prompt (Windows) or terminal (macOS) and type (install/upgrade Python packages) : 
        conda install pyqtgraph pyserial 

Linux
  1.  Type the command below by replacing username with that of your pc change permission of Anaconda3    sudo chown -R username:username /home/username/anaconda3
  2.  Open Anaconda3 terminal  and type (install/upgrade Python packages) : 
        conda install pyqtgraph pyserial
  3.  Set permission on serial port 
        sudo usermod -a -G uucp username
        sudo usermod -a -G dialout username
  4.  Logout and Login

### Qt Designer 
Qt Designer is a tool used to easily create graphical interfaces. To apply changes to the GUI, download and install it at : https://build-system.fman.io/qt-designer-download
    
The Qt Designer file for the software can be found here:  
    `...\OpenQCM\openQCM_Next_py_0.1.2_source\OPENQCM\openQCM\res\mainWindow_new.ui`
The file is saved under `.ui` extension.

### Link between Qt Designer interface and software functionnalities
To link the GUI to custom functionalities, the file extension has to be change. 
  1. Open Anaconda3 terminal (Windows) or terminal (macOS) and go to the following directory:
        `...\OpenQCM\openQCM_Next_py_0.1.2_source\OPENQCM\openQCM\res\`
  2. Type the command:
        `pyuic5 mainWindow_new.ui -o mainWindow_new_ui.py`
  3. In the same directory, type the command:
        `cp mainWindow_new_ui.py ../ui/`  
    This will copy the new file with the `.py` extension in the `..\openQCM\ui` directory.

### Adding new functionnalities
The file that has to be edited in order to add new functionnalities is the following:
    `...\OpenQCM\openQCM_Next_py_0.1.2_source\OPENQCM\openQCM\ui\mainWindow.py`
* The link between the GUI created on Qt Designer and `mainWindow.py` is made with the function `_configure_signals(self)`.
* To add a function linked to a graphical element on the GUI such as a button follow this example:
     * On Qt Designer, each element is given a name in the 'objectName' box in the 'Property Editor' panel. For the button used to start the injection of fluid, this button is called `pButton_Start_Injection`.
     * We will then associate a function to this button when it is clicked in `_configure_signals(self)`:
     `self.ui.pButton_Start_Injection.clicked.connect(self.Start_Injection)`  
     * Important: On Qt Designer, the text displayed on this button is 'Start Injection'. It is not the name used for the software.
* You can now edit your function. It will be launched when the button is clicked.
       
## INSA implemented functions
Si tu veux ajouter du code c'est comme ça
```python
  def updateViews1():
  	self._plt0.clear()
  	if self._get_source() != SourceType.multiscan:
  		self._plt1.clear()
  ```

## Usage
Start the application from Anaconda3 prompt (Windows) or terminal (macOS)
1.  Launch Anaconda3 prompt or terminal
2.  Browse to the openQCM Python software main directory:
          `...\openQCM_Next_py_0.1.2_source\OPENQCM`
3.  launch the python application main GUI by typing the command:
          `python app.py`
        


## Links
- [website] https://openqcm.com/
- [github]  https://github.com/openQCM 

- [website] http://lpcno.insa-toulouse.fr/

## Contact
openQCM
- [mail] info@openqcm.com

INSA
- [mail] paccard@insa-toulouse.fr
- [mail] a_lemair@insa-toulouse.fr
