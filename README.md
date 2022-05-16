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
  1.  Type the command below by replacing username with that of your pc change permission of Anaconda3    `sudo chown -R username:username /home/username/anaconda3`
  2.  Open Anaconda3 terminal  and type (install/upgrade Python packages) : 
        `conda install pyqtgraph pyserial`
  3.  Set permission on serial port 
        `sudo usermod -a -G uucp username`
        `sudo usermod -a -G dialout username`
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

### Adding new functionalities
The file that has to be edited in order to add new functionalities is the following:
    `...\OpenQCM\openQCM_Next_py_0.1.2_source\OPENQCM\openQCM\ui\mainWindow.py`
* The link between the GUI created on Qt Designer and `mainWindow.py` is made with the function `_configure_signals(self)`.
* To add a function linked to a graphical element on the GUI such as a button follow this example:
     * On Qt Designer, each element is given a name in the 'objectName' box in the 'Property Editor' panel. For the button used to start the injection of fluid, this button is called `pButton_Start_Injection`.
     * We will then associate a function to this button  in `_configure_signals(self)` when it is clicked:
     `self.ui.pButton_Start_Injection.clicked.connect(self.Start_Injection)`  
     * Important: On Qt Designer, the text displayed on this button is 'Start Injection'. It is not the name used for the software.
* You can now edit your function. It will be launched when the button is clicked.

## Usage
Start the application from Anaconda3 prompt (Windows) or terminal (macOS)
1.  Launch Anaconda3 prompt or terminal
2.  Browse to the openQCM Python software main directory:
          `...\openQCM_Next_py_0.1.2_source\OPENQCM`
3.  launch the python application main GUI by typing the command:
          `python app.py`
       
## INSA implemented functions
### Functions related to the usage of the injection equipment 

This function allow the program to get the status of the instrument, to know if they are connected and ready to be used, it initialises the connection. 

```Python
def InstrumentsInfo(self):
        print('')
        print('INSTRUMENTS INFORMATIONS:')
        # Detect all controllers
        SNs, types = fgt_detect()
        controllerCount = len(SNs)
        print('Number of controllers detected: {}'.format(controllerCount))

        # List all found controllers' serial number and type
        for i, sn in enumerate(SNs):
            print('Detected instrument at index: {}, ControllerSN: {}, type: {}' \
                  .format(i, sn, str(types[i])))

        print('')

        ## Initialize specific instruments
        # Initialize only specific instrument controllers here If you do not want
        # a controller in the list or if you want a specific order (e.g. LineUP
        # before MFCS instruments), rearrange parsed SN table
        fgt_init(SNs)

        ## Get the number of channels of each type

        # Get total number of initialized pressure channels
        print('Total number of pressure channels: {}'.format(fgt_get_pressureChannelCount()))

        # Get total number of initialized pressure channels
        print('Total number of sensor channels: {}'.format(fgt_get_sensorChannelCount()))

        # Get total number of initialized TTL channels
        print('Total number of TTL channels: {}'.format(fgt_get_TtlChannelCount()))

        # Get total number of initialized valve channels
        print('Total number of valve channels: {}'.format(fgt_get_valveChannelCount()))

        print('')

        ## Get detailed information about all controllers

        controllerInfoArray = fgt_get_controllersInfo()
        for i, controllerInfo in enumerate(controllerInfoArray):
            print('Controller info at index: {}'.format(i))
            print(controllerInfo)
            print('')

        ## Get detailed information about all pressure channels

        pressureInfoArray = fgt_get_pressureChannelsInfo()
        for i, pressureInfo in enumerate(pressureInfoArray):
            print('Pressure channel info at index: {}'.format(i))
            print(pressureInfo)
            print('')

        ## Get detailed information about all sensor channels

        sensorInfoArray, sensorTypeArray = fgt_get_sensorChannelsInfo()
        for i, sensorInfo in enumerate(sensorInfoArray):
            print('Sensor channel info at index: {}'.format(i))
            print(sensorInfo)
            print("Sensor type: {}".format(sensorTypeArray[i]))
            print('')

        ## Get detailed information about all TTL channels

        ttlInfoArray = fgt_get_TtlChannelsInfo()
        for i, ttlInfo in enumerate(ttlInfoArray):
            print('TTL channel info at index: {}'.format(i))
            print(ttlInfo)
            print('')

        valveInfoArray, valveTypeArray = fgt_get_valveChannelsInfo()
        for i, valveInfo in enumerate(valveInfoArray):
            print('Valve channel info at index: {}'.format(i))
            print(valveInfo)
            print("Valve type: {}".format(valveTypeArray[i]))
            print('')

        ## Close the session
        fgt_close()
```

The following functions are designed to start and stop the injection process, by setting the pressure either to the desired value during the desired time, or by setting it to 0. The start function initialises the injection. 

```Python
def Start_Injection(self):
        pressure = self.ui.spinBox_Pressure.value()
        global injection_time
        injection_time = self.ui.spinBox_Injection_Time.value()

        ## Initialize the session
        # This step is optional, if not called session will be automatically created
        fgt_init()

        # Set pressure to the given pressure on first pressure channel of the list
        # mbar is the default unit at initialization
        fgt_set_pressure(0, pressure)

        y = True
        t_end = time.time() + injection_time/1000
        while time.time() < t_end:
            if Stop_Injection_Bool and y:
                fgt_set_pressure(0, 0)
                y = False
        fgt_set_pressure(0, 0)
```

```Python
def Stop_Injection(self):
        global Stop_Injection_Bool
        Stop_Injection_Bool = True
        fgt_set_pressure(0,0)
```

### Functions related to the usage of the camera :

The following functions enable to turn the camera on, and to start to display the video on the GUI. 

This function turn the camera on. 

```Python
def start_camera(self):
        camera_port = int(self.ui.cBox_Camera_Port.currentText())
        if self.capture is None:
            self.capture = cv2.VideoCapture(camera_port) # argument may change with camera
            self.capture.set(cv2.CAP_PROP_FRAME_HEIGHT, 480)    # height of saved image   [480]
            self.capture.set(cv2.CAP_PROP_FRAME_WIDTH, 640)     # width of saved image    [640]
        self.timer.start()
```

This function displays then the video on the GUI directly. 

```Python
def displayImage(self, img, window=True):
        qformat = QtGui.QImage.Format_Indexed8
        if len(img.shape) == 3:
            if img.shape[2] == 4:
                qformat = QtGui.QImage.Format_RGBA8888
            else:
                qformat = QtGui.QImage.Format_RGB888
        outImage = QtGui.QImage(img, img.shape[1], img.shape[0], img.strides[0], qformat)
        outImage = outImage.rgbSwapped()
        outImage = outImage.mirrored(1,0)

        if window:
            # Set the size of the camera widget on the UI (200x200px)
            self.ui.imgLabel.setPixmap(QtGui.QPixmap.fromImage(outImage.scaled(200, 200)))
```

And finally this function updates the frame of the camera in the GUI. 

```Python
def update_frame(self):
        ret, image = self.capture.read()
        self.displayImage(image, True)
  ```

The following functions are designed to save on the computer pictures from the camera. The first one enables to take a single picture, while the second one takes several pictures according to the capture parameters selected.

```Python
def screenshot(self):
        flag, frame = self.capture.read()
        path = 'screenshots'
        if flag:
            name = "Screenshot_{}.png".format(datetime.now().strftime("%Y-%m-%d_%H-%M-%S-%f"))
            cv2.imwrite(os.path.join(path, name), frame)
            
```
  
```Python
  def capture_image(self):
        capture_rate = self.ui.spinBox_Capture_Rate.value()
        capture_time = self.ui.spinBox_Capture_Time.value()/1000  # capture time in ms on the ui but in sec in functions
        dt = 1/capture_rate

        dirname = "Recorded_images_{}".format(datetime.now().strftime("%Y-%m-%d_%H-%M-%S-%f"))   # creates new directory for every recording
        os.makedirs('recorded_images/{}'.format(dirname))
```
  
``` Python
    def record_image():
            flag, frame = self.capture.read()

            path = 'recorded_images/{}'.format(dirname)

            if flag:
                name = "Photoluminescence_{}.png".format(datetime.now().strftime("%Y-%m-%d_%H-%M-%S-%f"))
                cv2.imwrite(os.path.join(path, name), frame)

        thread1 = perpetualTimer(dt, record_image)
        thread1.start()
        T.sleep(capture_time)
        thread1.cancel()
```
### Function related to the simultaneous injection and image capture process

The GUI allows the operator to automate the recording of pictures during the injection process, using the functions below. 

The following function sets the recording parameters according to the ones chosen by the operator.

```Python
def Injection_Capture(self):
        capture_rate = self.ui.spinBox_Capture_Rate.value()
        injection_time = self.ui.spinBox_Injection_Time.value() / 1000                  # injection time
        time_before_injection = self.ui.spinBox_Time_Before_Injection.value() / 1000    # image capture time before injection
        time_after_injection = self.ui.spinBox_Time_After_Injection.value() / 1000      # image capture time after injection

        dt = 1 / capture_rate

        dirname = "Recorded_images_{}".format(datetime.now().strftime("%Y-%m-%d_%H-%M-%S-%f"))  # creates new directory for every recording
        os.makedirs('recorded_images/{}'.format(dirname))

        def record_image():
               flag, frame = self.capture.read()

               path = 'recorded_images/{}'.format(dirname)
               if flag:
                   name = "Photoluminescence_{}.png".format(datetime.now().strftime("%Y-%m-%d_%H-%M-%S-%f"))
                   cv2.imwrite(os.path.join(path, name), frame)

        def injection():
            pressure = self.ui.spinBox_Pressure.value()
            fgt_init()
            fgt_set_pressure(0, pressure)

        thread1 = perpetualTimer(dt, record_image)
        
        thread1.start()
        T.sleep(time_before_injection)
        injection()
        T.sleep(injection_time)
        fgt_set_pressure(0, 0)
        T.sleep(time_after_injection)
        thread1.cancel()
```


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
