# INSA add-ons created by Luca PACCARD and Arthur LEMAIRE

# from openQCM.ui.mainWindow_ui import Ui_Controls, Ui_Info, Ui_Plots

from __future__ import print_function  # For Fluigent

from openQCM.ui.mainWindow_new_ui import Ui_MainWindow

# from openQCM.ui.ui_controls import Ui_Controls
# from openQCM.ui.ui_info import Ui_Info
# from openQCM.ui.ui_plots import Ui_Plots

from pyqtgraph import AxisItem
import pyqtgraph as pg

from PyQt5 import QtCore, QtGui
from PyQt5.QtCore import *

# VER 0.1.2
# importing from PyQt5 or PySide2
'''
try:
    from PyQt5 import QtCore, QtGui
except:
    from PySide2 import  QtCore, QtGui
'''

from openQCM.core.worker import Worker
from openQCM.processors.Serial import SerialProcess
from openQCM.core.constants import Constants, SourceType, DateAxis, NonScientificAxis
from openQCM.ui.popUp import PopUp
from openQCM.common.logger import Logger as Log

import numpy as np
import sys
import serial

import time
from numpy import loadtxt
from openQCM.core.ringBuffer import RingBuffer

from time import sleep

###########################################################################
# INSA imports
###########################################################################
# Injection
from Fluigent.SDK import fgt_detect, fgt_init, fgt_close
from Fluigent.SDK import fgt_get_controllersInfo
from Fluigent.SDK import fgt_get_pressureChannelCount, fgt_get_pressureChannelsInfo
from Fluigent.SDK import fgt_get_sensorChannelCount, fgt_get_sensorChannelsInfo
from Fluigent.SDK import fgt_get_TtlChannelCount, fgt_get_TtlChannelsInfo
from Fluigent.SDK import fgt_get_valveChannelCount, fgt_get_valveChannelsInfo
from Fluigent.SDK import fgt_set_pressure, fgt_get_pressure, fgt_get_pressureRange

# Camera
import os
import cv2
from PyQt5 import QtCore, QtGui, QtWidgets, uic
from datetime import datetime
import time as T
from threading import Timer


########################################################################
# INSA : Global variables for injection functions
########################################################################
Injection_Bool = True
Stop_Injection_Bool = False


########################################################################
# INSA : perpetualTimer class for threaded events (simultaneous injection + image capture)
########################################################################

class perpetualTimer():

   def __init__(self,t,hFunction):
      self.t=t
      self.hFunction = hFunction
      self.thread = Timer(self.t,self.handle_function)

   def handle_function(self):
      self.hFunction()
      self.thread = Timer(self.t,self.handle_function)
      self.thread.start()

   def start(self):
      self.thread.start()

   def cancel(self):
      self.thread.cancel()


##########################################################################################
# Package that handles the UIs elements and connects to worker service to execute processes
##########################################################################################

TAG = ""  # "[MainWindow]"

class MainWindow(QtGui.QMainWindow):

    ###########################################################################
    # Initializes methods, values and sets the UI
    ###########################################################################
    def __init__(self, samples=Constants.argument_default_samples):

        #:param samples: Default samples shown in the plot :type samples: int.
        # to be always placed at the beginning, initializes some important methods
        QtGui.QMainWindow.__init__(self)

        self.ui = Ui_MainWindow()
        self.ui.setupUi(self)

        # Shared variables, initial values
        self._plt0 = None
        self._plt1 = None
        self._plt2 = None
        # TODO 2m self._plt3 DISSIPATION
        # TODO delete
        # self._plt3 = None
        self._plt4 = None

        self._pltD = None

        self._timer_plot = None
        self._readFREQ = None
        self._QCS_installed = None
        self._ser_control = 0
        self._ser_error1 = 0
        self._ser_error2 = 0
        self._ser_err_usb = 0

        # internet connection variable
        self._internet_connected = False

        # Reference variables
        self._reference_flag = False
        self._vector_reference_frequency = None
        self._vector_reference_dissipation = None
        self._vector_1 = None
        self._vector_2 = None

        # Instantiates a Worker class
        self.worker = Worker()

        # Populates comboBox for sources
        self.ui.cBox_Source.addItems(Constants.app_sources)

        # Init combo box for PID setting
        self.ui.cBox_PID.addItems(Constants.PID_default_settings)
        # set default value to #1 openqcm setting
        self.ui.cBox_PID.setCurrentIndex(Constants.PID_Setting_default_index)

        # Configures specific elements of the PyQtGraph plots
        self._configure_plot()

        # Configures specific elements of the QTimers
        self._configure_timers()

        # Configures the connections between signals and UI elements
        self._configure_signals()

        # Populates combo box for serial ports
        self._source_changed()
        # set calibration to default measurement mode
        self.ui.cBox_Source.setCurrentIndex(SourceType.calibration.value)

        # TODO delete sbox number of samples
        # self.ui.sBox_Samples.setValue(samples)  #samples

        # set style temperature indicator
        self.ui.indicator_temperature.setStyleSheet(
            "background-color:white; padding: 2 px; border-style: inset; border-color: gray; border-width: 1px;")
        # set style freqwuency and dissapation indicator
        self.ui.F0.setStyleSheet(
            "background-color:white; padding: 2 px; border-style: inset; border-color: gray; border-width: 1px;")
        self.ui.D0.setStyleSheet(
            "background-color:white; padding: 2 px; border-style: inset; border-color: gray; border-width: 1px;")
        self.ui.F3.setStyleSheet(
            "background-color:white; padding: 2 px; border-style: inset; border-color: gray; border-width: 1px;")
        self.ui.D3.setStyleSheet(
            "background-color:white; padding: 2 px; border-style: inset; border-color: gray; border-width: 1px;")
        self.ui.F5.setStyleSheet(
            "background-color:white; padding: 2 px; border-style: inset; border-color: gray; border-width: 1px;")
        self.ui.D5.setStyleSheet(
            "background-color:white; padding: 2 px; border-style: inset; border-color: gray; border-width: 1px;")
        self.ui.F7.setStyleSheet(
            "background-color:white; padding: 2 px; border-style: inset; border-color: gray; border-width: 1px;")
        self.ui.D7.setStyleSheet(
            "background-color:white; padding: 2 px; border-style: inset; border-color: gray; border-width: 1px;")
        self.ui.F9.setStyleSheet(
            "background-color:white; padding: 2 px; border-style: inset; border-color: gray; border-width: 1px;")
        self.ui.D9.setStyleSheet(
            "background-color:white; padding: 2 px; border-style: inset; border-color: gray; border-width: 1px;")

        # enable ui
        self._enable_ui(True)

        # disable temperature control at startup
        self._Temperature_PID_Setting_isEnabled(False)

        ###############################################################################################################
        self.get_web_info()
        # Gets the QCS installed on the device (not used now)
        # self._QCS_installed = PopUp.question_QCM(self, Constants.app_title, "Please choose the Quartz Crystal Resonator installed on the openQCM-1 Device (default 5MHz if exit)")

        # TODO my serial for temeprature communication
        self._my_serial = serial.Serial()

        # TODO multi frequency buffer definition # FREQUENCY
        self._frequency_buffer = RingBuffer(Constants.ring_buffer_samples)
        self._frequency_buffer_1 = RingBuffer(Constants.ring_buffer_samples)
        self._frequency_buffer_2 = RingBuffer(Constants.ring_buffer_samples)

        # TODO multi dissipation buffer definition  # DISSIPATION
        self._dissipation_buffer = RingBuffer(Constants.ring_buffer_samples)
        self._dissipation_buffer_1 = RingBuffer(Constants.ring_buffer_samples)
        self._dissipation_buffer_2 = RingBuffer(Constants.ring_buffer_samples)

        self.ui.label_Temperature_state.setStyleSheet(
            "background-color: yellow; border: 1px solid gray; border-radius:2px; padding: 2 px;")

        # multiscan array selector
        self.scan_selector = [0, 0, 0, 0, 0]

        # set T and PID to default values at init
        self._set_PID_T_default()

        # VER 0.1.2
        # multiscan y-range limit lists
        self._y_freq_max = [0, 0, 0, 0, 0]
        self._y_freq_min = [0, 0, 0, 0, 0]
        self._y_diss_max = [0, 0, 0, 0, 0]
        self._y_diss_min = [0, 0, 0, 0, 0]

    # https://stackoverflow.com/questions/63182608/colcount-not-working-for-legenditem-in-pyqtgraph-with-pyqt5-library
    # TODO legend horizontal layout
    def setColumnCount(self, legend, columnCount):
        '''
        change the orientation of all items of the legend
        '''

        def _addItemToLayout(legend, sample, label):
            col = legend.layout.columnCount()
            row = legend.layout.rowCount()
            if row:
                row -= 1
            nCol = legend.columnCount * 2
            # FIRST ROW FULL
            if col == nCol:
                for col in range(0, nCol, 2):
                    # FIND RIGHT COLUMN
                    if not legend.layout.itemAt(row, col):
                        break
                if col + 2 == nCol:
                    # MAKE NEW ROW
                    col = 0
                    row += 1
            legend.layout.addItem(sample, row, col)
            legend.layout.addItem(label, row, col + 1)

        legend.columnCount = columnCount
        legend.rowCount = int(len(legend.items) / columnCount)
        for i in range(legend.layout.count() - 1, -1, -1):
            legend.layout.removeAt(i)  # clear layout
        for sample, label in legend.items:
            _addItemToLayout(legend, sample, label)
        legend.updateSize()

    ###########################################################################
    # Starts the acquisition of the selected serial port
    ###########################################################################
    def start(self):

        import os
        os.system('cls' if os.name == 'nt' else 'clear')

        # This function is connected to the clicked signal of the Start button.
        # print("")
        print(TAG, 'Clicked START')
        Log.i(TAG, "Clicked START")

        # INSA INJECTION : boolean reset to true each time the "START" button is pushed
        global Injection_Bool
        Injection_Bool = True

        # TODO Disable temperature control button
        # DEV
        # =============================================================================
        #         self.ui.pButton_Tswitch_OFF.setEnabled(False)
        #         self.ui.pButton_Tswitch_ON.setEnabled(False)
        # =============================================================================

        # Instantiates process
        self.worker = Worker(QCS_on=self._QCS_installed,
                             port=self.ui.cBox_Port.currentText(),
                             speed=self.ui.cBox_Speed.currentText(),
                             samples=Constants.argument_default_samples,
                             source=self._get_source(),
                             export_enabled=False)

        # SINGLE
        # ---------------------------------------------------------------------
        if self.worker.start():
            # Gets frequency range
            # self._readFREQ = self.worker.get_frequency_range()

            # Duplicate frequencies
            self._reference_flag = False
            # self._vector_reference_frequency = list(self._readFREQ)
            self._reference_value_frequency = 0
            self._reference_value_dissipation = 0

            # init frequency reference array values
            self._reference_value_frequency_array = [0, 0, 0, 0, 0]
            # init dissipation reference array values
            self._reference_value_dissipation_array = [0, 0, 0, 0, 0]

            self._labelref1 = "not set"
            self._labelref2 = "not set"
            # progressbar variables
            self._completed = 0
            self._ser_control = 0
            # error variables
            self._ser_error1 = 0
            self._ser_error2 = 0
            self._ser_err_usb = 0
            ##### other useful location #########
            # self.get_web_info()
            #####

            # SINGLE
            # -----------------------------------------------------------------
            if self._get_source() == SourceType.serial:
                # TODO DELETE
                # overtones_number = len(self.worker.get_source_speeds(SourceType.serial))

                # for single scan Gets frequency range
                self._readFREQ = self.worker.get_frequency_range()

                self._vector_reference_frequency = list(self._readFREQ)
                self._overtones_number_all = len(self.worker.get_source_speeds(SourceType.serial))

                # TODO set the quartz sensor

                if (float(self.worker.get_source_speeds(SourceType.serial)[
                              self._overtones_number_all - 1]) > 4e+06 and float(
                        self.worker.get_source_speeds(SourceType.serial)[self._overtones_number_all - 1]) < 6e+06):
                    label_quartz = "@5MHz_QCM"
                elif (float(self.worker.get_source_speeds(SourceType.serial)[
                                self._overtones_number_all - 1]) > 9e+06 and float(
                        self.worker.get_source_speeds(SourceType.serial)[self._overtones_number_all - 1]) < 11e+06):
                    label_quartz = "@10MHz_QCM"

                #  TODO set the legend in single mode
                overtone_selected = self._overtones_number_all - self.ui.cBox_Speed.currentIndex() - 1
                # TODO PRINT THE LEGEND
                # frequency
                self._plt2.plot(
                    pen=pg.mkPen(color=Constants.plot_color_multi[overtone_selected], width=Constants.plot_line_width),
                    name=Constants.name_legend[overtone_selected])
                # dissipation
                self._pltD.plot(
                    pen=pg.mkPen(color=Constants.plot_color_multi[overtone_selected], width=Constants.plot_line_width),
                    name=Constants.name_legend[overtone_selected])

                # =============================================================================
                #                 # VER 0.1.2
                #                 # add phase plot additional axis
                #                 self._plt0.scene().addItem(self._plt1)
                # =============================================================================
                # clear plot
                self.clear()
                self._plt1.clear()

            # CALIBRATION
            # -----------------------------------------------------------------
            elif self._get_source() == SourceType.calibration:

                label_quartz = self.ui.cBox_Speed.currentText()

                # VER 0.1.2
                # add phase plot additional axis
                self._plt0.scene().addItem(self._plt1)
                # clear plot
                self.clear()
                self._plt1.clear()



            # MULTISCAN
            # -----------------------------------------------------------------
            elif self._get_source() == SourceType.multiscan:

                # TODO get number of overtones and do nothing apparently
                self._overtones_number_all = len(self.worker.get_source_speeds(SourceType.serial))

                # TODO get the quartz crystal fundamental frequency
                if (float(self.worker.get_source_speeds(SourceType.multiscan)[
                              self._overtones_number_all - 1]) > 4e+06 and float(
                        self.worker.get_source_speeds(SourceType.multiscan)[self._overtones_number_all - 1]) < 6e+06):
                    label_quartz = "@5MHz_QCM"
                elif (float(self.worker.get_source_speeds(SourceType.multiscan)[
                                self._overtones_number_all - 1]) > 9e+06 and float(
                        self.worker.get_source_speeds(SourceType.multiscan)[self._overtones_number_all - 1]) < 11e+06):
                    label_quartz = "@10MHz_QCM"

                # TODO redefine the array here
                self._arr = np.zeros((self._overtones_number_all, Constants.ring_buffer_samples))

                # legend
                for idx in range(self._overtones_number_all):
                    # TODO PRINT THE LEGEND
                    # frequency
                    self._plt2.plot(
                        pen=pg.mkPen(color=Constants.plot_color_multi[idx], width=Constants.plot_line_width),
                        name=Constants.name_legend[idx])
                    # dissipation
                    self._pltD.plot(
                        pen=pg.mkPen(color=Constants.plot_color_multi[idx], width=Constants.plot_line_width),
                        name=Constants.name_legend[idx])

                # init radio button
                self.ui.radioBtn_F0.setChecked(True)
                self.ui.radioBtn_F3.setChecked(True)
                self.ui.radioBtn_F5.setChecked(True)
                self.ui.radioBtn_F7.setChecked(True)
                self.ui.radioBtn_F9.setChecked(True)

                self._update_scan_selector()

                # =============================================================================
                #                 # VER 0.1.2
                #                 # add phase plot additional axis
                #                 self._plt0.scene().removeItem(self._plt1)
                # =============================================================================

                # clear plot
                self.clear()

            # SET TIMER UPDATE
            # -----------------------------------------------------------------
            self._timer_plot.start(Constants.plot_update_ms)

            # CALL UPDATE PLOT
            # -----------------------------------------------------------------
            self._timer_plot.timeout.connect(self._update_plot)  # moved from _configure_timers mothod
            self._enable_ui(False)

            # VER 0.1.2
            # auto scale frequency and dissipation plt
            self._plt2.enableAutoRange(axis='y', enable=True)
            self._pltD.enableAutoRange(axis='y', enable=True)

            if self._get_source() == SourceType.calibration:
                self.ui.pButton_Clear.setEnabled(False)  # insert
                self.ui.pButton_Reference.setEnabled(False)  # insert
                self.ui.pButton_Reference_Not.setEnabled(False)

        else:
            print(TAG, "Warning: port is not available!")
            Log.i(TAG, "Warning: port is not available")
            PopUp.warning(self, Constants.app_title,
                          "Warning: Selected Port [{}] is not available!".format(self.ui.cBox_Port.currentText()))

    ###########################################################################
    # Stops the acquisition of the selected serial port
    ###########################################################################
    def stop(self):

        # This function is connected to the clicked signal of the Stop button.
        self.ui.infostatus.setStyleSheet('background: white; padding: 1px; border: 1px solid #cccccc')
        self.ui.infostatus.setText("<font color=#000000 > Program Status Stanby</font>")
        self.ui.infobar.setText("<font color=#0000ff > Infobar </font>")

        # TODO Enable temperature control button
        self.ui.pButton_Tswitch_OFF.setEnabled(True)
        self.ui.pButton_Tswitch_ON.setEnabled(True)

        # remove legend item
        for idx in range(self._overtones_number_all):
            self._legend_f.removeItem(Constants.name_legend[idx])
            self._legend_D.removeItem(Constants.name_legend[idx])

        # set temperature to default value
        # self.ui.doubleSpinBox_Temperature.setValue( Constants.Temperature_Set_Value )

        self._timer_plot.stop()
        self._enable_ui(True)
        self.worker.stop()

        # add a delay to prevent the error caused by the serial com port open
        time.sleep(1)

        # turn off the peltier
        self.Temperature_Control_OFF()

        # add a little delay
        time.sleep(0.5)
        # reset temperature and PID to default values
        self._set_PID_T_default()

        # set pid setting combo box to default factory
        self.ui.cBox_PID.setCurrentIndex(Constants.PID_Setting_default_index)

    ###########################################################################
    # SET TEMPERATURE
    ###########################################################################
    def temperatureSet(self):
        '''
        print ("SET TEMPERATURE")
        var = self.ui.doubleSpinBox_Temperature.value()
        print ('T' + str(int(var)) + '\n')

        self._my_serial.port = self.ui.cBox_Port.currentText()
        self._my_serial.baudrate = Constants.serial_default_speed #115200
        self._my_serial.stopbits = serial.STOPBITS_ONE
        self._my_serial.bytesize = serial.EIGHTBITS
        self._my_serial.timeout = Constants.serial_timeout_ms
        self._my_serial.writetimeout = Constants.serial_writetimeout_ms

        # Gets the state of the serial port
#        if not self._my_serial.isOpen():
#            # OPENS the serial port
#            self._my_serial.open()
#            cmd = 'T' + str(int(var)) + '\n'
#            self._my_serial.write(cmd.encode())
##            print(self._my_serial.readall)
##            buffer = self._my_serial.read(self._my_serial.in_waiting).decode(Constants.app_encoding)
##            print(buffer)
        # OPENS the serial port
        self._my_serial.open()
        cmd = 'T' + str(int(var)) + '\n'
        self._my_serial.write(cmd.encode())
        self._my_serial.close()
        '''

        var = self._get_temperature()
        cmd = 'T' + str(int(var)) + '\n'

        print("Set Temperature =  ", var / 1000)

        # serial port parameter
        self._my_serial.port = self.ui.cBox_Port.currentText()
        self._my_serial.baudrate = Constants.serial_default_speed  # 115200
        self._my_serial.stopbits = serial.STOPBITS_ONE
        self._my_serial.bytesize = serial.EIGHTBITS
        self._my_serial.timeout = Constants.serial_timeout_ms
        self._my_serial.writetimeout = Constants.serial_writetimeout_ms

        # check if a process is NOT running
        if (self.worker.is_running() == False):
            # open the serial port
            self._my_serial.open()
            # write set temperature command
            self._my_serial.write(cmd.encode())
            # close serial
            self._my_serial.close()

    def _get_temperature(self):
        _var = self.ui.doubleSpinBox_Temperature.value() * 1000
        _path = Constants.manual_frequencies_path
        # np.savetxt( _path,  np.row_stack( [_var, 1] ) )
        # np.savetxt( _path,  [_var])
        # return self.ui.doubleSpinBox_Temperature.value()

        _var_cycling_time = self.ui.spinBox_Cycling_Time.value()
        _var_P_share = self.ui.spinBox_P_Share.value()
        _var_I_Share = self.ui.spinBox_I_Share.value()
        _var_D_Share = self.ui.spinBox_D_Share.value()
        _var_bool = 1

        # VER 0.1.2
        # get temperature control boolean
        param = loadtxt(Constants.manual_frequencies_path)
        _ctrl_bool = param[6]

        np.savetxt(_path, np.row_stack(
            [_var, _var_cycling_time, _var_P_share, _var_I_Share, _var_D_Share, _var_bool, _ctrl_bool]), fmt='%d')

        return _var

    # TEMPERATURE CONTROL FUNCTION
    ###########################################################################

    def Temperature_Control_ON(self):

        # change the indicator color
        self.ui.label_Temperature_state.setStyleSheet(
            "background-color: rgba(0, 142, 192, 0.4); border: 1px solid gray; border-radius: 2px;  padding: 2 px;")

        print("Temperature Control ON")
        # set the temperature control ONLY in NOT measuring mode
        self._my_serial.port = self.ui.cBox_Port.currentText()
        self._my_serial.baudrate = Constants.serial_default_speed  # 115200
        self._my_serial.stopbits = serial.STOPBITS_ONE
        self._my_serial.bytesize = serial.EIGHTBITS
        self._my_serial.timeout = Constants.serial_timeout_ms
        self._my_serial.writetimeout = Constants.serial_writetimeout_ms

        # Gets the state of the serial port
        # if not self._my_serial.isOpen():
        # VER 0.1.2
        # check if a process is NOT running to verify the serial is available
        if (self.worker.is_running() == False):
            # OPENS the serial port
            self._my_serial.open()
            var = 1
            cmd = 'X' + str(int(var)) + '\n'
            self._my_serial.write(cmd.encode())
            self._my_serial.close()

        # =============================================================================
        #         # VER 0.1.2 TODO
        #         elif  ( self.worker.is_running() == True ):
        #             # DEV
        #             print ("the worker is running ")
        # =============================================================================

        # VER 0.1.2
        # save temperature control boolean to config file
        _path = Constants.manual_frequencies_path
        param = loadtxt(Constants.manual_frequencies_path)
        _ctrl_bool = 1
        np.savetxt(_path, np.row_stack([param[0], param[1], param[2], param[3], param[4], param[5], _ctrl_bool]),
                   fmt='%d')

        # enable - disable UI control
        self._Temperature_PID_Setting_isEnabled(True)

    def Temperature_Control_OFF(self):

        # change the led color
        self.ui.label_Temperature_state.setStyleSheet(
            "background-color: yellow; border: 1px solid gray; border-radius:2px; padding: 2 px;")
        # set temoerature control to default
        self.ui.doubleSpinBox_Temperature.setValue(Constants.Temperature_Set_Value)

        print("Temperature Control OFF ")
        # set the temperature control ONLY in NOT measuring mode
        self._my_serial.port = self.ui.cBox_Port.currentText()
        self._my_serial.baudrate = Constants.serial_default_speed  # 115200
        self._my_serial.stopbits = serial.STOPBITS_ONE
        self._my_serial.bytesize = serial.EIGHTBITS
        self._my_serial.timeout = Constants.serial_timeout_ms
        self._my_serial.writetimeout = Constants.serial_writetimeout_ms

        # Gets the state of the serial port
        # if not self._my_serial.isOpen():

        # DEV
        # check if a process is NOT running to verify the serial is available
        if (self.worker.is_running() == False):
            # open the serial port
            self._my_serial.open()
            var = 0
            cmd = 'X' + str(int(var)) + '\n'
            self._my_serial.write(cmd.encode())
            self._my_serial.close()

        # =============================================================================
        #         elif  ( self.worker.is_running() == True ):
        #             #DEV
        #             print ("waring: unable to send disable temperature control message")
        # =============================================================================

        # DEV version 0.1.1.d
        # save temperature control boolean to config file
        _path = Constants.manual_frequencies_path
        param = loadtxt(Constants.manual_frequencies_path)
        _ctrl_bool = 0
        np.savetxt(_path, np.row_stack([param[0], param[1], param[2], param[3], param[4], param[5], _ctrl_bool]),
                   fmt='%d')

        # enable - disable UI control
        self._Temperature_PID_Setting_isEnabled(False)

    def _Temperature_PID_Setting_isEnabled(self, my_bool):
        # PID set button
        self.ui.pButton_PID_Set.setEnabled(my_bool)
        # Tempeature set button
        self.ui.pButton_Temperature_Set.setEnabled(my_bool)
        # PID param control
        self.ui.spinBox_Cycling_Time.setEnabled(my_bool)
        self.ui.spinBox_P_Share.setEnabled(my_bool)
        self.ui.spinBox_I_Share.setEnabled(my_bool)
        self.ui.spinBox_D_Share.setEnabled(my_bool)
        # temperature param control
        self.ui.doubleSpinBox_Temperature.setEnabled(my_bool)
        # default parameter selection
        self.ui.cBox_PID.setEnabled(my_bool)

    def _Overtone_radioBtn_isEnabled(self, my_bool):
        self.ui.radioBtn_F0.setEnabled(my_bool)
        self.ui.radioBtn_F3.setEnabled(my_bool)
        self.ui.radioBtn_F5.setEnabled(my_bool)
        self.ui.radioBtn_F7.setEnabled(my_bool)
        self.ui.radioBtn_F9.setEnabled(my_bool)

    #########################################################################
    # INSA Functions
    #########################################################################
    # Injection
    # TODO : Emergency Stop Injection button
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

    def Stop_Injection(self):
        global Stop_Injection_Bool
        Stop_Injection_Bool = True
        fgt_set_pressure(0,0)

    def Pressure_Reset(self):
        fgt_set_pressure(0,0)
        fgt_close()

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
        #Injection()
        # Wait injection_time (in seconds) before setting pressure to 0
        #time.sleep(injection_time / 1000)

    def Injection(self):
        global injection_time
        y = True
        t_end = time.time() + injection_time/1000
        while time.time() < t_end:
            if Stop_Injection_Bool and y:
                fgt_set_pressure(0, 0)
                y = False
        fgt_set_pressure(0, 0)

    # Camera
    def start_camera(self):
        # TODO : Fix error when fist port chosen isn't available
        camera_port = int(self.ui.cBox_Camera_Port.currentText())
        if self.capture is None:
            self.capture = cv2.VideoCapture(camera_port) # argument may change with camera
            self.capture.set(cv2.CAP_PROP_FRAME_HEIGHT, 480)    # height of saved image   [480]
            self.capture.set(cv2.CAP_PROP_FRAME_WIDTH, 640)     # width of saved image    [640]
        self.timer.start()

    def screenshot(self):
        flag, frame = self.capture.read()
        path = 'screenshots'
        if flag:
            name = "Screenshot_{}.png".format(datetime.now().strftime("%Y-%m-%d_%H-%M-%S-%f"))
            cv2.imwrite(os.path.join(path, name), frame)

    def update_frame(self):
        ret, image = self.capture.read()
        self.displayImage(image, True)

    def capture_image(self):
        capture_rate = self.ui.spinBox_Capture_Rate.value()
        capture_time = self.ui.spinBox_Capture_Time.value()/1000  # capture time in ms on the ui but in sec in functions
        dt = 1/capture_rate

        dirname = "Recorded_images_{}".format(datetime.now().strftime("%Y-%m-%d_%H-%M-%S-%f"))   # creates new directory for every recording
        os.makedirs('recorded_images/{}'.format(dirname))

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

    def displayImage(self, img, window=True):

        """
        # resize image
        scale_percent = 80  # percent of original size
        width = int(img.shape[1] * scale_percent / 100)
        height = int(img.shape[0] * scale_percent / 100)
        """

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
            #self.ui.imgLabel.setPixmap(QtGui.QPixmap.scaled(width,height).fromImage(resized))

    # Simultaneous Injection + Image Capture
    def Injection_Capture(self):
        capture_rate = self.ui.spinBox_Capture_Rate.value()
        injection_time = self.ui.spinBox_Injection_Time.value() / 1000                  # injection time
        time_before_injection = self.ui.spinBox_Time_Before_Injection.value() / 1000    # image capture time before injection
        time_after_injection = self.ui.spinBox_Time_After_Injection.value() / 1000      # image capture time after injection

        dt = 1 / capture_rate

        dirname = "Recorded_images_{}".format(datetime.now().strftime("%Y-%m-%d_%H-%M-%S-%f"))  # creates new directory for every recording
        os.makedirs('recorded_images/{}'.format(dirname))

        def record_image():
            #print("recording")

            flag, frame = self.capture.read()

            path = 'recorded_images/{}'.format(
                dirname)
            if flag:
                name = "Photoluminescence_{}.png".format(datetime.now().strftime("%Y-%m-%d_%H-%M-%S-%f"))
                cv2.imwrite(os.path.join(path, name), frame)


        def injection():
            #print("Hello world!")


            pressure = self.ui.spinBox_Pressure.value()
            fgt_init()
            fgt_set_pressure(0, pressure)


        thread1 = perpetualTimer(dt, record_image)
        #thread2 = perpetualTimer(injection_time, injection)

        print("go capture")
        thread1.start()
        T.sleep(time_before_injection)
        print("go injection")
        #thread2.start()
        injection()

        T.sleep(injection_time)
        #thread2.cancel()
        fgt_set_pressure(0, 0)
        print("fini injection")

        T.sleep(time_after_injection)
        thread1.cancel()
        print("fini capture")

    #TODO : Enable injection / camera settings with buttons
    #def _Injection_isEnabled(self, my_bool):


    # PID CONTROL FUNCTION
    ##########################################################################

    def PID_Set(self):
        # TODO PID SET HERE
        print("Setting PID Parameter")
        self._get_PID()

        # serial port parameter
        self._my_serial.port = self.ui.cBox_Port.currentText()
        self._my_serial.baudrate = Constants.serial_default_speed  # 115200
        self._my_serial.stopbits = serial.STOPBITS_ONE
        self._my_serial.bytesize = serial.EIGHTBITS
        self._my_serial.timeout = Constants.serial_timeout_ms
        self._my_serial.writetimeout = Constants.serial_writetimeout_ms

        # check if a process is NOT running to verify the serial is available
        if (self.worker.is_running() == False):

            # open the serial port
            self._my_serial.open()

            # get PID paramter from UI
            _var_cycling_time = self.ui.spinBox_Cycling_Time.value()
            _var_P_share = self.ui.spinBox_P_Share.value()
            _var_I_Share = self.ui.spinBox_I_Share.value()
            _var_D_Share = self.ui.spinBox_D_Share.value()

            # set PID parameters
            cycling_time_msg = 'C' + str(int(_var_cycling_time)) + '\n'
            # VER 0.1.2 add a short sleep for communication
            sleep(0.1)
            self._my_serial.write(cycling_time_msg.encode())
            sleep(0.1)

            P_Share_msg = 'P' + str(int(_var_P_share)) + '\n'
            self._my_serial.write(P_Share_msg.encode())
            sleep(0.1)

            I_Share_msg = 'I' + str(int(_var_I_Share)) + '\n'
            self._my_serial.write(I_Share_msg.encode())
            sleep(0.1)

            D_Share_msg = 'D' + str(int(_var_D_Share)) + '\n'
            self._my_serial.write(D_Share_msg.encode())
            sleep(0.1)
            # close the serial port
            self._my_serial.close()

        # VER 0.1.2 TODO
        else:
            print("the worker is still running ")

    def _get_PID(self):
        # TODO get pid parameters from main gui
        _var_cycling_time = self.ui.spinBox_Cycling_Time.value()
        _var_P_share = self.ui.spinBox_P_Share.value()
        _var_I_Share = self.ui.spinBox_I_Share.value()
        _var_D_Share = self.ui.spinBox_D_Share.value()
        # get temperature param
        _var = self.ui.doubleSpinBox_Temperature.value() * 1000
        # change the setting boolean variable
        _var_bool = 1
        # get current value of control temperature boolean
        param = loadtxt(Constants.manual_frequencies_path)
        _ctrl_bool = param[6]

        _path = Constants.manual_frequencies_path
        np.savetxt(_path, np.row_stack(
            [_var, _var_cycling_time, _var_P_share, _var_I_Share, _var_D_Share, _var_bool, _ctrl_bool]), fmt='%d')

    def _set_PID_T_default(self):
        _path = Constants.manual_frequencies_path
        # save default file in config.ini file
        np.savetxt(_path, np.row_stack([Constants.Temperature_Set_Value * 1000,
                                        Constants.cycling_time_default,
                                        Constants.P_share_default,
                                        Constants.I_share_default,
                                        Constants.D_share_default,
                                        Constants.PID_boolean_default,
                                        Constants.CTRL_boolean_default]), fmt='%d')

        # update indicator to default value
        self.ui.doubleSpinBox_Temperature.setValue(Constants.Temperature_Set_Value)
        self.ui.spinBox_Cycling_Time.setValue(Constants.cycling_time_default)
        self.ui.spinBox_P_Share.setValue(Constants.P_share_default)
        self.ui.spinBox_I_Share.setValue(Constants.I_share_default)
        self.ui.spinBox_D_Share.setValue(Constants.D_share_default)

    def _PID_setting_changed(self):
        print("PID setting changed ")
        _my_index = self.ui.cBox_PID.currentIndex()

        self.ui.spinBox_Cycling_Time.setValue(Constants.cycling_time_setting[_my_index])
        self.ui.spinBox_P_Share.setValue(Constants.P_share_setting[_my_index])
        self.ui.spinBox_I_Share.setValue(Constants.I_share_setting[_my_index])
        self.ui.spinBox_D_Share.setValue(Constants.D_share_setting[_my_index])

        # TODO add automatic pid setting push button

    ###########################################################################
    # Overrides the QTCloseEvent,is connected to the close button of the window
    ###########################################################################
    def closeEvent(self, evnt):

        #:param evnt: QT evnt.
        if self.worker.is_running():
            print(TAG, 'Window closed without stopping the capture, application will stop...')
            Log.i(TAG, "Window closed without stopping the capture, application will stop...")
            self.stop()
            # self.ControlsWin.close()
            # self.PlotsWin.close()
            # self.InfoWin.close()
            # evnt.accept()

        # QtGui.QApplication.quit()

    # =============================================================================
    #         res = PopUp.question(self, Constants.app_title, "Are you sure you want to quit openQCM application now?")
    #         if res:
    #            # self.close()
    #            # evnt.accept()
    #            QtGui.QApplication.quit()
    #         else:
    #            evnt.ignore()
    # =============================================================================

    ###########################################################################
    # Enables or disables the UI elements of the window.
    ###########################################################################
    def _enable_ui(self, enabled):

        #:param enabled: The value to be set for the UI elements :type enabled: bool
        self.ui.cBox_Port.setEnabled(enabled)
        self.ui.cBox_Speed.setEnabled(enabled)
        self.ui.pButton_Start.setEnabled(enabled)

        # TODO delete or implement export txt file
        # self.ui.chBox_export.setEnabled(enabled)

        self.ui.cBox_Source.setEnabled(enabled)
        self.ui.pButton_Stop.setEnabled(not enabled)

        # TODO delete the sample
        # self.ui.sBox_Samples.setEnabled(not enabled) #insert

        self.ui.pButton_Clear.setEnabled(not enabled)
        self.ui.pButton_Reference.setEnabled(not enabled)
        self.ui.pButton_Reference_Not.setEnabled(not enabled)

        # DEV

    # =============================================================================
    #         self.ui.pButton_Tswitch_OFF.setEnabled(enabled)
    #         self.ui.pButton_Tswitch_ON.setEnabled(enabled)
    # =============================================================================

    ###########################################################################
    # Configures specific elements of the PyQtGraph plots.
    ###########################################################################
    def _configure_plot(self):

        # ----------------------------------------------------------------------
        # set background color background="#0c2c36"
        # TEMPERATURE and SWEEP PLOT #0c2c36
        self.ui.plt.setBackground(background=Constants.plot_background_color)
        # FREQUENY PLOT
        self.ui.pltB.setBackground(background=Constants.plot_background_color)
        # DISSIPATION PLOT
        self.ui.pltD.setBackground(background=Constants.plot_background_color)
        # ----------------------------------------------------------------------

        # defines the graph title
        title1 = "Real-Time Plot: Amplitude / Phase"
        title2 = "Real-Time Plot: Resonance Frequency "
        title3 = "Real-Time Plot: Temperature"
        # ----------------------------------------------------------------------

        # Configures elements of the PyQtGraph plots: amplitude
        self.ui.plt.setAntialiasing(True)
        self.ui.pltB.setAntialiasing(True)
        self.ui.pltD.setAntialiasing(True)

        self._xaxis_sweep = NonScientificAxis(orientation='bottom')
        self._xaxis_sweep.enableAutoSIPrefix(False)

        '''
        -----------------------------------------------------------------------
        SWEEP PLOT AMPLITUDE PHASE
        -----------------------------------------------------------------------
        '''
        self._plt0 = self.ui.plt.addPlot(row=0, col=0, title=title1, **{'font-size': '10pt'},
                                         axisItems={"bottom": self._xaxis_sweep})
        # self._plt0.showGrid(x=True, y=True)
        self._plt0.setLabel('bottom', 'Frequency', units='Hz')
        self._plt0.setLabel('left', 'Amplitude', units='dB', color=Constants.plot_title_color, **{'font-size': '10pt'})

        '''
        # Configures elements of the PyQtGraph plots: phase
        self._plt1 = self.u2.plt.addPlot(row=1, col=1, title= "Real-Time Plot: Phase", **{'font-size':'12pt'})
        self._plt1.showGrid(x=True, y=True)
        self._plt1.setLabel('bottom', 'Samples', units='n')
        self._plt1.setLabel('left', 'Phase', units='deg')
        '''
        # --------------------------------------------------------------------------------------------------------------
        # Configures elements of the PyQtGraph plots: Multiple Plot amplitude and phase
        self._plt1 = pg.ViewBox()
        self._plt0.showAxis('right')
        self._plt0.scene().addItem(self._plt1)
        self._plt0.getAxis('right').linkToView(self._plt1)
        self._plt1.setXLink(self._plt0)
        self._plt0.enableAutoRange(axis='y', enable=True)
        self._plt1.enableAutoRange(axis='y', enable=True)
        self._plt0.setLabel('right', 'Phase', units='deg', color=Constants.plot_title_color, **{'font-size': '10pt'})

        # VER 0.1.2
        # editing pyqtgraph context menu
        # https://groups.google.com/g/pyqtgraph/c/h-dyr0l6yZU/m/NpMQxh-jf5cJ
        # get rid of 'Plot Options'
        self._plt0.ctrlMenu = None
        self._plt1.ctrlMenu = None
        # get rid of 'Export'
        self._plt0.scene().contextMenu = None
        self._plt1.scene().contextMenu = None

        '''
        -----------------------------------------------------------------------
        FREQUENCY PLOT
        -----------------------------------------------------------------------
        '''
        # --------------------------------------------------------------------------------------------------------------
        # Configures elements of the PyQtGraph plots: resonance
        self._yaxis = NonScientificAxis(orientation='left')
        self._yaxis.enableAutoSIPrefix(False)
        # self._yaxis.setTickSpacing(levels=[(280, 0),(25, 0), (10, 0)]) #(20,1, None)
        self._xaxis = DateAxis(orientation='bottom')

        '''
        TODO 2m
        '''
        # Configures elements of the PyQtGraph plots: dissipatin
        self._yaxisD = NonScientificAxis(orientation='left')
        self._yaxisD.enableAutoSIPrefix(False)
        # self._yaxis.setTickSpacing(levels=[(280, 0),(25, 0), (10, 0)]) #(20,1, None)
        self._xaxisD = DateAxis(orientation='bottom')

        '''
        self._plt2 = self.PlotsWin.ui2.pltB.addPlot(row=0, col=2, title= title2, **{'font-size':'12pt'}, axisItems={"bottom":self._xaxis, 'left':self._yaxis})
        '''
        '''
        TODO 2m
        FREQUENCY PLOT
        '''
        self._plt2 = self.ui.pltB.addPlot(row=0, col=2, title=title2, **{'font-size': '12pt'},
                                          axisItems={"bottom": self._xaxis, 'left': self._yaxis})

        # self._plt2.showGrid(x=True, y=True)
        self._plt2.setLabel('bottom', 'Time', units='s')
        self._plt2.setLabel('left', 'Resonance Frequency', units='Hz', color=Constants.plot_title_color,
                            **{'font-size': '10pt'})

        # https://pyqtgraph.readthedocs.io/en/latest/graphicsItems/plotitem.html#pyqtgraph.PlotItem.addLegend
        self._legend_f = self._plt2.addLegend()

        # change the orientation of all items of the legend
        # self._legend_f.setColumnCount(5)
        # TODO LEGEND
        # self.setColumnCount(self._legend_f, 5)

        # VER 0.1.2
        # editing pyqtgraph context menu
        # https://groups.google.com/g/pyqtgraph/c/h-dyr0l6yZU/m/NpMQxh-jf5cJ
        # get rid of 'Plot Options'
        self._plt2.ctrlMenu = None
        # get rid of 'Export'
        self._plt2.scene().contextMenu = None

        '''
        TODO 2m
        DISSIPATION PLOT
        '''
        # self._pltD = self.ui.pltD.addPlot(row=0, col=1, title= "Real-Time Plot: Dissipation", **{'font-size':'12pt'}, axisItems={"bottom":self._xaxisD, 'left':self._yaxisD})
        self._pltD = self.ui.pltD.addPlot(row=0, col=1, title="Real-Time Plot: Dissipation", **{'font-size': '12pt'},
                                          axisItems={"bottom": self._xaxisD})
        self._pltD.setLabel('bottom', 'Time', units='s')
        self._pltD.setLabel('left', 'Dissipation', units='', color=Constants.plot_title_color, **{'font-size': '10pt'})

        # https://pyqtgraph.readthedocs.io/en/latest/graphicsItems/plotitem.html#pyqtgraph.PlotItem.addLegend
        self._legend_D = self._pltD.addLegend()
        # change the orientation of all items of the legend
        # self._legend_D.setColumnCount(5)
        # TODO LEGEND
        # self.setColumnCount(self._legend_D, 5)

        # VER 0.1.2
        # editing pyqtgraph context menu
        # https://groups.google.com/g/pyqtgraph/c/h-dyr0l6yZU/m/NpMQxh-jf5cJ
        # get rid of 'Plot Options'
        self._pltD.ctrlMenu = None
        # get rid of 'Export'
        self._pltD.scene().contextMenu = None

        # CONNECT PLOT
        # ---------------------------------------------------------------------
        self._pltD.setXLink(self._plt2)
        '''
        TODO 2m
        it is not TRIVIAL to connect y axis with different range
        # self._pltD.setYLink(self._plt2)
        '''

        # --------------------------------------------------------------------------------------------------------------
        # Configures elements of the PyQtGraph plots: Multiple Plot resonance frequency and dissipation

        # Configures elements of the PyQtGraph plots: temperature
        self._plt4 = self.ui.plt.addPlot(row=0, col=1, title=title3,
                                         axisItems={'bottom': DateAxis(orientation='bottom')})
        # self._plt4.showGrid(x=True, y=True)

        # do not autoscale y axis
        self._plt4.enableAutoRange(axis='y', enable=True)

        # VER 0.1.2
        # change the Temperature Y-range to 5 - 45 Â°C
        self._plt4.setYRange(5, 45, padding=0)

        self._plt4.setLabel('bottom', 'Time', units='s')
        self._plt4.setLabel('left', 'Temperature', units='°C', color=Constants.plot_title_color,
                            **{'font-size': '10pt'})

        # VER 0.1.2
        # editing pyqtgraph context menu
        # https://groups.google.com/g/pyqtgraph/c/h-dyr0l6yZU/m/NpMQxh-jf5cJ
        # get rid of 'Plot Options'
        self._plt4.ctrlMenu = None
        # get rid of 'Export'
        self._plt4.scene().contextMenu = None

    ###########################################################################
    # Configures specific elements of the QTimers
    ###########################################################################
    def _configure_timers(self):

        self._timer_plot = QtCore.QTimer(self)
        # self._timer_plot.timeout.connect(self._update_plot) #moved to start method

    ###########################################################################
    # Configures the connections between signals and UI elements
    ###########################################################################
    def _configure_signals(self):

        self.ui.pButton_Start.clicked.connect(self.start)
        self.ui.pButton_Stop.clicked.connect(self.stop)
        self.ui.pButton_Clear.clicked.connect(self.clear)
        self.ui.pButton_Reference.clicked.connect(self.reference)
        self.ui.pButton_Reference_Not.clicked.connect(self.reference_not)

        # TODO delete sample box
        # self.ui.sBox_Samples.valueChanged.connect(self._update_sample_size)

        self.ui.cBox_Source.currentIndexChanged.connect(self._source_changed)
        self.ui.pButton_Temperature_Set.clicked.connect(self.temperatureSet)
        self.ui.pButton_PID_Set.clicked.connect(self.PID_Set)

        self.ui.cBox_PID.currentIndexChanged.connect(self._PID_setting_changed)

        # Temperature control button
        self.ui.pButton_Tswitch_ON.clicked.connect(self.Temperature_Control_ON)
        self.ui.pButton_Tswitch_OFF.clicked.connect(self.Temperature_Control_OFF)

        ##############################################################################
        ########################## INSA Buttons ######################################
        ##############################################################################
        # Injection
        self.ui.pButton_InstrumentsInfo.clicked.connect(self.InstrumentsInfo)
        self.ui.pButton_Stop_Injection.clicked.connect(self.Stop_Injection)
        self.ui.pButton_Start_Injection.clicked.connect(self.Start_Injection)

        # Camera
        self.ui.pButton_Start_Camera.clicked.connect(self.start_camera)
        self.ui.pButton_Screenshot.clicked.connect(self.screenshot)
        self.ui.pButton_Capture_Image.clicked.connect(self.capture_image)
        self.ui.imgLabel.setScaledContents(True)
        self.capture = None
        self.timer = QtCore.QTimer(self, interval=5)    # Timer for camera refresh rate (in ms)
        self.timer.timeout.connect(self.update_frame)
        self._image_counter = 0

        # Camera + Injection
        self.ui.pButton_Injection_Capture.clicked.connect(self.Injection_Capture)

        # self.ui.pButton_Switch_ON.clicked.connect(self.Switch_Valve_ON)
        # self.ui.OFF.clicked.connect(self.Switch_Valve_OFF)




        ##########################################################################

        # Frequency scan selector
        self.ui.radioBtn_F0.clicked.connect(self._update_scan_selector)
        self.ui.radioBtn_F3.clicked.connect(self._update_scan_selector)
        self.ui.radioBtn_F5.clicked.connect(self._update_scan_selector)
        self.ui.radioBtn_F7.clicked.connect(self._update_scan_selector)
        self.ui.radioBtn_F9.clicked.connect(self._update_scan_selector)

        '''
        #--------
        self.InfoWin.ui3.pButton_Download.clicked.connect(self.start_download)
        '''
        '''
        TODO 2m
        connect start download
        '''

    ###########################################################################
    # Updates the sample size of the plot (now not used)
    ###########################################################################
    def _update_sample_size(self):

        # This function is connected to the valueChanged signal of the sample Spin Box.
        if self.worker is not None:
            self.worker.reset_buffers(Constants.argument_default_samples)

    ###########################################################################
    # Updates and redraws the graphics in the plot.
    ###########################################################################

    def _update_plot(self):
        global Injection_Bool

        # This function is connected to the timeout signal of a QTimer
        self.worker.consume_queue1()
        self.worker.consume_queue2()
        self.worker.consume_queue3()
        self.worker.consume_queue4()
        # TODO note that data is logged here, when self.worker.consume_queue5() is called
        self.worker.consume_queue5()
        # general error queue
        self.worker.consume_queue6()

        self.worker.consume_queue_F_multi()
        self.worker.consume_queue_D_multi()

        self.worker.consume_queue_A_multi()

        # SINGLE
        # --------------------------------------------------------------------
        if self._get_source() == SourceType.serial:
            vector1 = self.worker.get_d1_buffer()
            vector2 = self.worker.get_d2_buffer()
            vectortemp = self.worker.get_d3_buffer()

            # TODO changed the get error number of elemts
            # self._ser_error1,self._ser_error2, self._ser_control,self._ser_err_usb = self.worker.get_ser_error()
            self._ser_error1, self._ser_error2, self._ser_control, self._ser_err_usb, self._overtone_number = self.worker.get_ser_error()

            # print(self._ser_err_usb, end='\r')
            # if self._ser_err_usb <=1:
            if vector1.any:
                # progressbar
                if self._ser_control <= Constants.environment:
                    self._completed = self._ser_control * 2

                if str(vector1[0]) == 'nan' and not self._ser_error1 and not self._ser_error2:
                    label1 = 'processing...'
                    label2 = 'processing...'
                    label3 = 'processing...'
                    labelstatus = 'Processing'
                    '''
                  self.ControlsWin.ui1.infostatus.setStyleSheet('background: #ffff00; padding: 1px; border: 1px solid #cccccc') #ff8000
                  '''
                    '''
                  TODO 2m
                  '''
                    self.ui.infostatus.setStyleSheet(
                        'background: #ffff00; padding: 1px; border: 1px solid #cccccc')  # ff8000

                    color_err = '#000000'
                    labelbar = 'Please wait, processing early data...'

                elif (str(vector1[0]) == 'nan' and (self._ser_error1 or self._ser_error2)):
                    if self._ser_error1 and self._ser_error2:
                        label1 = ""
                        label2 = ""
                        label3 = ""
                        labelstatus = 'Warning'
                        color_err = '#ff0000'
                        labelbar = 'Warning: unable to apply half-power bandwidth method, lower and upper cut-off frequency not found'
                        '''
                        self.ControlsWin.ui1.infostatus.setStyleSheet('background: #ff0000; padding: 1px; border: 1px solid #cccccc')
                        '''
                        '''
                        TODO 2m
                        '''
                        self.ui.infostatus.setStyleSheet('background: #ff0000; padding: 1px; border: 1px solid #cccccc')

                    elif self._ser_error1:
                        label1 = ""
                        label2 = ""
                        label3 = ""
                        labelstatus = 'Warning'
                        color_err = '#ff0000'
                        labelbar = 'Warning: unable to apply half-power bandwidth method, lower cut-off frequency (left side) not found'
                        '''
                        self.ControlsWin.ui1.infostatus.setStyleSheet('background: #ff0000; padding: 1px; border: 1px solid #cccccc')
                        '''
                        '''
                        TODO 2m
                        '''
                        # self.ControlsWin.ui1.infostatus.setStyleSheet('background: #ff0000; padding: 1px; border: 1px solid #cccccc')
                    elif self._ser_error2:
                        label1 = ""
                        label2 = ""
                        label3 = ""
                        labelstatus = 'Warning'
                        color_err = '#ff0000'
                        labelbar = 'Warning: unable to apply half-power bandwidth method, upper cut-off frequency (right side) not found'
                        '''
                        self.ControlsWin.ui1.infostatus.setStyleSheet('background: #ff0000; padding: 1px; border: 1px solid #cccccc')
                        '''
                        '''
                        TODO 2m
                        '''
                        self.ui.infostatus.setStyleSheet('background: #ff0000; padding: 1px; border: 1px solid #cccccc')
                else:
                    if not self._ser_error1 and not self._ser_error2:
                        if not self._reference_flag:
                            d1 = float("{0:.2f}".format(vector1[0]))
                            d2 = float("{0:.4f}".format(vector2[0] * 1e6))
                            d3 = float("{0:.2f}".format(vectortemp[0]))
                        else:
                            a1 = vector1[0] - self._reference_value_frequency
                            a2 = vector2[0] - self._reference_value_dissipation
                            d1 = float("{0:.2f}".format(a1))
                            d2 = float("{0:.4f}".format(a2 * 1e6))
                            d3 = float("{0:.2f}".format(vectortemp[0]))
                        label1 = str(d1) + " Hz"
                        label2 = str(d2) + "e-06"
                        label3 = str(d3) + " °C"
                        labelstatus = 'Monitoring'
                        color_err = '#000000'
                        labelbar = 'Monitoring!'
                        '''
                      self.ControlsWin.ui1.infostatus.setStyleSheet('background: #00ff72; padding: 1px; border: 1px solid #cccccc')
                      '''
                        '''
                      TODO 2m
                      '''
                        self.ui.infostatus.setStyleSheet('background: #00ff72; padding: 1px; border: 1px solid #cccccc')

                    else:
                        if self._ser_error1 and self._ser_error2:
                            label1 = "-"
                            label2 = "-"
                            label3 = "-"
                            labelstatus = 'Warning'
                            color_err = '#ff0000'
                            labelbar = 'Warning: unable to apply half-power bandwidth method, lower and upper cut-off frequency not found'
                            '''
                        self.ControlsWin.ui1.infostatus.setStyleSheet('background: #ff0000; padding: 1px; border: 1px solid #cccccc')
                        '''
                            '''
                        TODO 2m
                        '''
                            self.ui.infostatus.setStyleSheet(
                                'background: #ff0000; padding: 1px; border: 1px solid #cccccc')

                        elif self._ser_error1:
                            label1 = "-"
                            label2 = "-"
                            label3 = "-"
                            labelstatus = 'Warning'
                            color_err = '#ff0000'
                            labelbar = 'Warning: unable to apply half-power bandwidth method, lower cut-off frequency (left side) not found'
                            '''
                        self.ControlsWin.ui1.infostatus.setStyleSheet('background: #ff0000; padding: 1px; border: 1px solid #cccccc')
                        '''
                            '''
                        TODO 2m
                        '''
                            self.ui.infostatus.setStyleSheet(
                                'background: #ff0000; padding: 1px; border: 1px solid #cccccc')

                        elif self._ser_error2:
                            label1 = "-"
                            label2 = "-"
                            label3 = "-"
                            labelstatus = 'Warning'
                            color_err = '#ff0000'
                            labelbar = 'Warning: unable to apply half-power bandwidth method, upper cut-off frequency (right side) not found'
                            '''
                        self.ControlsWin.ui1.infostatus.setStyleSheet('background: #ff0000; padding: 1px; border: 1px solid #cccccc')
                        '''
                            '''
                        TODO 2m
                        '''
                            self.ui.infostatus.setStyleSheet(
                                'background: #ff0000; padding: 1px; border: 1px solid #cccccc')
                '''
               self.InfoWin.ui3.l6a.setText("<font color=#0000ff > Temperature </font>" + label3)
               self.InfoWin.ui3.l6.setText("<font color=#0000ff > Dissipation </font>" + label2)
               self.InfoWin.ui3.l7.setText("<font color=#0000ff > Resonance Frequency </font>" + label1)
               '''
                '''
               TODO 2m set info
               '''

                '''
               self.ControlsWin.ui1.infostatus.setText("<font color=#000000 > Program Status </font>" + labelstatus)
               self.ControlsWin.ui1.infobar.setText("<font color=#0000ff> Infobar </font><font color={}>{}</font>".format(color_err,labelbar))
               # progressbar
               self.ControlsWin.ui1.progressBar.setValue(self._completed+2)
               '''
                '''
               TODO 2m
               '''
                self.ui.infostatus.setText("<font color=#000000 > Program Status </font>" + labelstatus)
                self.ui.infobar.setText(
                    "<font color=#0000ff> Infobar </font><font color={}>{}</font>".format(color_err, labelbar))
                # progressbar
                self.ui.progressBar.setValue(self._completed + 2)

            # elif self._ser_err_usb >1:
            # PopUp.warning(self, Constants.app_title, "Warning: USB cable device disconnected!")
            # self.stop()

        # CALIBRATION: dynamic info in infobar at run-time
        # ---------------------------------------------------------------------
        elif self._get_source() == SourceType.calibration:
            # flag for terminating calibration
            stop_flag = 0
            '''
            self.ControlsWin.ui1.pButton_Stop.setEnabled(False)
            '''
            '''
            TODO 2m
            '''
            self.ui.pButton_Stop.setEnabled(False)

            vector1 = self.worker.get_value1_buffer()
            # vector2[0] and vector3[0] flag error
            vector2 = self.worker.get_t3_buffer()
            vector3 = self.worker.get_d3_buffer()
            # print(vector1[0],vector2[0],vector3[0])
            label1 = 'not available'
            label2 = 'not available'
            label3 = 'not available'
            labelstatus = 'Calibration Processing'
            color_err = '#000000'
            labelbar = 'The operation might take just over a minute to complete... please wait...'

            '''
            self.ControlsWin.ui1.infostatus.setStyleSheet('background: #ffff00; padding: 1px; border: 1px solid #cccccc')
            '''
            '''
            TODO 2m
            '''
            self.ui.infostatus.setStyleSheet('background: #ffff00; padding: 1px; border: 1px solid #cccccc')

            # progressbar
            error1, error2, error3, self._ser_control, self._overtone_number = self.worker.get_ser_error()
            if self._ser_control < (Constants.calib_sections):
                self._completed = (self._ser_control / (Constants.calib_sections)) * 100
            # calibration buffer empty
            # if vector1[0]== 0 and vector3[0]==1:
            if error1 == 1 and vector3[0] == 1:
                label1 = 'not available'
                label2 = 'not available'
                label3 = 'not available'
                color_err = '#ff0000'
                labelstatus = 'Calibration Warning'
                '''
              self.ControlsWin.ui1.infostatus.setStyleSheet('background: #ff0000; padding: 1px; border: 1px solid #cccccc')
              '''
                '''
              TODO 2m
              '''
                self.ui.infostatus.setStyleSheet('background: #ff0000; padding: 1px; border: 1px solid #cccccc')

                labelbar = 'Calibration Warning: empty buffer! Please, repeat the Calibration after disconnecting/reconnecting Device!'
                stop_flag = 1
            # calibration buffer empty and ValueError from the serial port
            # elif vector1[0]== 0 and vector2[0]==1:
            elif error1 == 1 and vector2[0] == 1:
                label1 = 'not available'
                label2 = 'not available'
                label3 = 'not available'
                color_err = '#ff0000'
                labelstatus = 'Calibration Warning'
                '''
              self.ControlsWin.ui1.infostatus.setStyleSheet('background: #ff0000; padding: 1px; border: 1px solid #cccccc')
              '''
                '''
              TODO 2m
              '''
                self.ui.infostatus.setStyleSheet('background: #ff0000; padding: 1px; border: 1px solid #cccccc')

                labelbar = 'Calibration Warning: empty buffer/ValueError! Please, repeat the Calibration after disconnecting/reconnecting Device!'
                stop_flag = 1
            # calibration buffer not empty
            # elif vector1[0]!= 0:
            elif error1 == 0:
                label1 = 'not available'
                label2 = 'not available'
                label3 = 'not available'
                labelstatus = 'Calibration Processing'
                color_err = '#000000'
                labelbar = 'The operation might take just over a minute to complete... please wait...'
                if vector2[0] == 0 and vector3[0] == 0:
                    labelstatus = 'Calibration Success'
                    '''
                 self.ControlsWin.ui1.infostatus.setStyleSheet('background: #00ff72; padding: 1px; border: 1px solid #cccccc')
                 '''
                    '''
                 TODO 2m
                 '''
                    self.ui.infostatus.setStyleSheet('background: #00ff72; padding: 1px; border: 1px solid #cccccc')

                    color_err = '#000000'
                    labelbar = 'Calibration Success for baseline correction!'
                    stop_flag = 1
                    # print(self._k) #progressbar value 143
                elif vector2[0] == 1 or vector3[0] == 1:
                    color_err = '#ff0000'
                    labelstatus = 'Calibration Warning'
                    '''
                 self.ControlsWin.ui1.infostatus.setStyleSheet('background: #ff0000; padding: 1px; border: 1px solid #cccccc')
                 '''
                    '''
                 TODO 2m
                 '''
                    self.ui.infostatus.setStyleSheet('background: #ff0000; padding: 1px; border: 1px solid #cccccc')

                    if vector2[0] == 1:
                        labelbar = 'Calibration Warning: ValueError or generic error during signal acquisition. Please, repeat the Calibration'
                        stop_flag = 1  ##
                    elif vector3[0] == 1:
                        labelbar = 'Calibration Warning: unable to identify fundamental peak or apply peak detection algorithm. Please, repeat the Calibration!'
                        stop_flag = 1  ##
            '''
            self.InfoWin.ui3.l6a.setText("<font color=#0000ff>  Dissipation </font>" + label3)
            self.InfoWin.ui3.l6.setText("<font color=#0000ff>  Dissipation </font>" + label2)
            self.InfoWin.ui3.l7.setText("<font color=#0000ff>  Resonance Frequency </font>" + label1)
            '''
            '''
            TODO 2m set info label in main ui
            '''

            '''
            self.ControlsWin.ui1.infostatus.setText("<font color=#000000> Program Status </font>" + labelstatus)
            self.ControlsWin.ui1.infobar.setText("<font color=#0000ff> Infobar </font><font color={}>{}</font>".format(color_err,labelbar))
            # progressbar -------------
            self.ControlsWin.ui1.progressBar.setValue(self._completed+10)
            '''
            '''
            TODO 2m
            '''
            self.ui.infostatus.setText("<font color=#000000> Program Status </font>" + labelstatus)
            self.ui.infobar.setText(
                "<font color=#0000ff> Infobar </font><font color={}>{}</font>".format(color_err, labelbar))
            # progressbar -------------
            self.ui.progressBar.setValue(self._completed + 10)

            # terminate the  calibration (simulate clicked stop)
            if stop_flag == 1:
                self._timer_plot.stop()
                self._enable_ui(True)
                self.worker.stop()
            '''
            # Amplitude plot
            self._plt0.clear()
            #self._plt0.plot(list(self._xdict.keys()),self.worker.get_value1_buffer(),pen=Constants.plot_colors[0])
            self._plt0.plot(self.worker.get_value1_buffer(),pen=Constants.plot_colors[0])

            # Phase plot
            self._plt1.clear()
            self._plt1.plot(self.worker.get_value2_buffer(),pen=Constants.plot_colors[1])
            '''

        # MULTISCAN:
        # ---------------------------------------------------------------------
        elif self._get_source() == SourceType.multiscan:
            vector1 = self.worker.get_d1_buffer()

            # TODO update plot
            self._ser_error1, self._ser_error2, self._ser_control, self._ser_err_usb, self._overtone_number = self.worker.get_ser_error()

            # =============================================================================
            #             # TODO DEBUG STOP SOFTWARE
            #             if self._ser_error1:
            #                 print ("WARNING:: general error of kind #1 ")
            #             if self._ser_error2:
            #                 print ("WARNING:: general error of kind #2 ")
            #             if self._ser_err_usb:
            #                 print("WARNING: Main Window general error on serial ")
            # =============================================================================

            if vector1.any:
                # progressbar
                if self._ser_control <= Constants.environment:
                    self._completed = self._ser_control * 2

                    # VER 0.1.2
                    # Optimize and update infobar and infostatus in multiscan mode
                    labelstatus = 'Processing'
                    color_err = '#000000'
                    labelbar = 'Please wait, processing early data...'
                    self.ui.infostatus.setText("<font color=#000000> Program Status </font>" + labelstatus)
                    self.ui.infobar.setText(
                        "<font color=#0000ff> Infobar </font><font color={}>{}</font>".format(color_err, labelbar))

                # progressbar
                self.ui.progressBar.setValue(self._completed + 2)

            if self._ser_control == Constants.environment:
                # clear plt reset buffer at the end of processing early data
                self.clear()

                # VER 0.1.2
                # Optimize and update infobar and infostatus in multiscan mode
                labelstatus = 'Monitoring'
                color_err = '#000000'
                labelbar = 'Monitoring multiscan frequency and disspation '
                self.ui.infostatus.setText("<font color=#000000> Program Status </font>" + labelstatus)
                self.ui.infobar.setText(
                    "<font color=#0000ff> Infobar </font><font color={}>{}</font>".format(color_err, labelbar))

        # REFERENCE SET
        # ---------------------------------------------------------------------
        if self._reference_flag:

            # SINGLE
            # -----------------------------------------------------------------
            if self._get_source() == SourceType.serial:

                # =============================================================================
                #                 self._plt2.setLabel('left', 'Resonance Frequency', units='Hz', color=Constants.plot_colors[6], **{'font-size':'10pt'})
                #                 self._plt2.setLabel('right', 'Dissipation', units='', color=Constants.plot_colors[7], **{'font-size':'10pt'})
                # =============================================================================
                '''
                self.InfoWin.ui3.inforef1.setText("<font color=#0000ff > Ref. Frequency </font>" + self._labelref1)
                self.InfoWin.ui3.inforef2.setText("<font color=#0000ff > Ref. Dissipation </font>" + self._labelref2)
                '''
                '''
                TODO 2m set reset frequency abnd dissipation
                '''

                # AMPLITUDE and PHASE
                # -------------------------------------------------------------
                def updateViews1():
                    self._plt0.clear()
                    # VER 0.1.2
                    # software freeze when interacting with software when resizing GUI window
                    if self._get_source() != SourceType.multiscan:
                        self._plt1.clear()
                    self._plt1.setGeometry(self._plt0.vb.sceneBoundingRect())
                    self._plt1.linkedViewChanged(self._plt0.vb, self._plt1.XAxis)

                # updates for multiple plot y-axes
                updateViews1()
                self._plt0.vb.sigResized.connect(updateViews1)
                self._plt0.plot(x=self._readFREQ, y=self.worker.get_value1_buffer(), pen=Constants.plot_colors[0])
                self._plt1.addItem(
                    pg.PlotCurveItem(x=self._readFREQ, y=self.worker.get_value2_buffer(), pen=Constants.plot_colors[1]))

                # frequency and dissipation update view
                def updateViews2():
                    self._plt2.clear()
                    self._pltD.clear()

                updateViews2()
                self._plt2.vb.sigResized.connect(updateViews2)

                #  TODO set the legend in single mode
                overtone_selected = self._overtones_number_all - self.ui.cBox_Speed.currentIndex() - 1

                # FREQUENCY and DISSIPATION
                # -------------------------------------------------------------
                self._vector_1 = np.array(self.worker.get_d1_buffer()) - self._reference_value_frequency
                self._vector_2 = np.array(self.worker.get_d2_buffer()) - self._reference_value_dissipation

                # VER 0.1.2
                # get y_freq and y_dissipation max value
                y_freq_single_max = np.nanmax(self._vector_1)
                y_freq_single_min = np.nanmin(self._vector_1)
                y_diss_single_max = np.nanmax(self._vector_2)
                y_diss_single_min = np.nanmin(self._vector_2)

                # VER 0.1.2 
                # set the y-range of dissipation and frequency axis 
                try:
                    self._plt2.setYRange(y_freq_single_min - 100, y_freq_single_max + 100, padding=0)
                    self._pltD.setYRange(y_diss_single_min - 0.000001, y_diss_single_max + 0.000001, padding=0)
                except:
                    pass

                self._plt2.plot(x=self.worker.get_t1_buffer(), y=self._vector_1,
                                pen=pg.mkPen(color=Constants.plot_color_multi[overtone_selected],
                                             width=Constants.plot_line_width))
                self._pltD.plot(x=self.worker.get_t1_buffer(), y=self._vector_2,
                                pen=pg.mkPen(color=Constants.plot_color_multi[overtone_selected],
                                             width=Constants.plot_line_width))

                # update frequency and dissipation indicator
                self._update_indicator_F_single(overtone_selected, self._vector_1)
                self._update_indicator_D_single(overtone_selected, self._vector_2)

                # Prevent the user from zooming/panning out of this specified region
                # TODO set limit of y range axis
                # =============================================================================
                #
                #                 if self._get_source() == SourceType.serial:
                #                    #dy = [(value, str(value)) for value in (range(int(min(self._readFREQ)), int(max(self._readFREQ)+1)))]
                #                    #self._yaxis.setTicks([dy])
                #                    #tickBottom = {self._readFREQ[250]:self._readFREQ[0],self._readFREQ[-1]:self._readFREQ[250]}
                #                    #self._yaxis.setTicks([tickBottom.items()])
                #                    self._plt2.setLimits(yMax=self._vector_reference_frequency[-1],yMin=self._vector_reference_frequency[0], minYRange=5)
                #     #               self._plt3.setLimits(yMax=self._vector_reference_dissipation[-1],yMin=self._vector_reference_dissipation[0], minYRange=1e-7)
                #                    self._plt4.setLimits(yMax=50,yMin=-10)
                #     #            self._plt3.addItem(pg.PlotCurveItem(self.worker.get_t2_buffer(),self._vector_2, pen=Constants.plot_colors[7]))
                # =============================================================================

                # TEMPERATURE
                # -------------------------------------------------------------
                self._plt4.clear()
                # do not autoscale y
                self._plt4.enableAutoRange(axis='y', enable=True)

                # set temperature y range
                # VER 0.1.2
                self._plt4.setYRange(5, 45, padding=0)

                # get temperature buffer
                y_temperature = self.worker.get_d3_buffer()
                self._plt4.plot(x=self.worker.get_t3_buffer(), y=y_temperature, pen=Constants.plot_colors[4])

                # set temperature current value
                label_indicator_temperature = float("{0:.2f}".format(y_temperature[0]))
                self.ui.indicator_temperature.setText(str(label_indicator_temperature))

            # MULTISCAN
            # -----------------------------------------------------------------
            elif self._get_source() == SourceType.multiscan:

                def updateViews_multi():
                    self._plt1.setGeometry(self._plt0.vb.sceneBoundingRect())
                    self._plt1.linkedViewChanged(self._plt0.vb, self._plt1.XAxis)

                ''' -----------------------------------------------------------
                # AMPLITUDE
                # ------------------------------------------------------------
                # Loads frequencies from calibration file
                peaks_mag = self.load_frequencies_file()
                # get sweep frequency range for each overtone
                x_sweep_axis = self.worker.get_frequency_range_multi(Constants.argument_default_samples, self._overtone_number)
                # center the sweep frequency axis around each overtone frequency peak
                x_sweep_axis = x_sweep_axis - peaks_mag[self._overtone_number]

                # 10M
                if (peaks_mag[0] >9e+06 and peaks_mag[0]<11e+06):
                   # set XY constant in multiscan mode
                   self._plt0.setXRange(-(Constants.L10_5th_overtone + 1000), Constants.R10_5th_overtone + 1000, padding = 0)
                   self._plt0.setYRange(-5, 30, padding = 0)

                # TODO 5M
                if (peaks_mag[0] >4e+06 and peaks_mag[0]<6e+06):
                   self._plt0.setXRange( -(Constants.L5_7th_overtone + 1000), Constants.R5_7th_overtone + 1000, padding = 0 )


                # updates for multiple plot y-axes
                updateViews_multi()
                self._plt0.vb.sigResized.connect(updateViews_multi)

                # TODO AMPLITUDE PLOT
                # clear plot at fundamental sweep
                if self._overtone_number == 0:
                   self._plt0.clear()

                # TODO AMPLITUDE PLOT
                if self.scan_selector[self._overtone_number] == True:
                    self._plt0.plot ( x = x_sweep_axis, y = self.worker.get_value1_buffer(), pen = Constants.plot_color_multi[self._overtone_number] )
                 -----------------------------------------------------------'''

                # AMPLITUDE
                # -------------------------------------------------------------
                # Loads frequencies from calibration file
                peaks_mag = self.load_frequencies_file()

                # VER 0.1.2 do not set XY limit
                # =============================================================================
                #                 # 10M
                #                 if (peaks_mag[0] >9e+06 and peaks_mag[0]<11e+06):
                #                     # set XY constant in multiscan mode
                #                     self._plt0.setXRange(-(Constants.L10_5th_overtone + 1000), Constants.R10_5th_overtone + 1000, padding = 0)
                #                     # self._plt0.setYRange(-10, 40, padding = 0)
                #                 # TODO 5M
                #                 if (peaks_mag[0] >4e+06 and peaks_mag[0]<6e+06):
                #                     # set XY constant in multiscan mode
                #                     self._plt0.setXRange( -(Constants.L5_7th_overtone + 1000), Constants.R5_7th_overtone + 1000, padding = 0 )
                #                     # self._plt0.setYRange(-5, 30, padding = 0)
                # =============================================================================

                # updates for multiple plot y-axes
                updateViews_multi()
                self._plt0.vb.sigResized.connect(updateViews_multi)

                self._plt0.clear()

                # loop on the lines
                # TODO introduce this control because I've some trouble in reset buffer ???
                if self._ser_control > Constants.environment:
                    for idx in range(self._overtones_number_all):
                        # get and scale frequency axis
                        x_sweep_axis = self.worker.get_F_Sweep_values_buffer(idx) - peaks_mag[idx]
                        # get amplitude axis
                        y_sweep_axis = self.worker.get_A_values_buffer(idx)

                        # plot sweep
                        # TODO there is something stange in worker.reset_buffers()
                        if self.scan_selector[idx] == True and isinstance(x_sweep_axis, np.ndarray):
                            self._plt0.plot(x=x_sweep_axis, y=y_sweep_axis, pen=Constants.plot_color_multi[idx])

                # FREQUENCY and DISSIPATION
                # ------------------------------------------------------------

                self._plt2.clear()

                self._pltD.clear()

                for idx in range(self._overtones_number_all):
                    if self.scan_selector[idx] == True:
                        # get time axis
                        time_axis_new = self.worker.get_time_values_buffer(idx)

                        # get y frequency and dissipation axis
                        y_freq = np.array(self.worker.get_F_values_buffer(idx)) - self._reference_value_frequency_array[
                            idx]
                        y_diss = np.array(self.worker.get_D_values_buffer(idx)) - \
                                 self._reference_value_dissipation_array[idx]

                        # VER 0.1.2
                        # get y_freq and y_dissipation max and min values array
                        self._y_freq_max[idx] = np.nanmax(y_freq)
                        self._y_freq_min[idx] = np.nanmin(y_freq)
                        self._y_diss_max[idx] = np.nanmax(y_diss)
                        self._y_diss_min[idx] = np.nanmin(y_diss)

                        # get the frequency and dissipation max and min on overtones
                        y_freq_max = max(self._y_freq_max)
                        y_freq_min = min(self._y_freq_min)
                        y_diss_max = max(self._y_diss_max)
                        y_diss_min = min(self._y_diss_min)

                        # print ("GET MAXIMUM VALUES FREQ AND DISSIP = ")
                        # print (y_freq_max), print (y_diss_max)

                        try:
                            # VER 0.1.2
                            # set the y-range of dissipation and frequency axis
                            self._plt2.setYRange(y_freq_min - 100, y_freq_max + 100, padding=0)
                            self._pltD.setYRange(y_diss_min - 0.000005, y_diss_max + 0.000005, padding=0)
                        except:
                            pass

                        # plot frequency and dissipation data
                        self._plt2.plot(x=time_axis_new, y=y_freq, pen=pg.mkPen(color=Constants.plot_color_multi[idx],
                                                                                width=Constants.plot_line_width))
                        self._pltD.plot(x=time_axis_new, y=y_diss, pen=pg.mkPen(color=Constants.plot_color_multi[idx],
                                                                                width=Constants.plot_line_width))

                        self._update_indicator_F(idx, y_freq)
                        self._update_indicator_D(idx, y_diss)

                    else:
                        dummy = [0]
                        self._update_indicator_F(idx, dummy)
                        self._update_indicator_D(idx, dummy)

                        # VER 0.1.2
                        # reset limits to get the correct y-ranges
                        self._y_freq_max[idx] = 0
                        self._y_freq_min[idx] = 0
                        self._y_diss_max[idx] = 0
                        self._y_diss_min[idx] = 0

                # TEMPERATURE
                # ------------------------------------------------------------
                self._plt4.clear()
                # do not autoscale y
                self._plt4.enableAutoRange(axis='y', enable=True)

                # set temperature y range
                # VER 0.1.2
                self._plt4.setYRange(5, 45, padding=0)

                # get temperature buffer
                y_temperature = self.worker.get_d3_buffer()

                self._plt4.plot(x=self.worker.get_t3_buffer(), y=y_temperature, pen=Constants.plot_colors[4])

                # set temperature current value
                label_indicator_temperature = float("{0:.2f}".format(y_temperature[0]))
                self.ui.indicator_temperature.setText(str(label_indicator_temperature))

        # REFERENCE NOT SET
        # ---------------------------------------------------------------------
        else:

            self._plt2.setLabel('left', 'Resonance Frequency', units='Hz', color=Constants.plot_title_color,
                                **{'font-size': '10pt'})

            # define update views for amplitude plt
            def updateViews1():
                self._plt0.clear()
                # VER 0.1.2
                # software freeze when interacting with software when resizing GUI window
                if self._get_source() != SourceType.multiscan:
                    self._plt1.clear()
                self._plt1.setGeometry(self._plt0.vb.sceneBoundingRect())
                self._plt1.linkedViewChanged(self._plt0.vb, self._plt1.XAxis)

            def updateViews_multi():
                self._plt1.setGeometry(self._plt0.vb.sceneBoundingRect())
                self._plt1.linkedViewChanged(self._plt0.vb, self._plt1.XAxis)

            # CALIBRATION
            # -----------------------------------------------------------------
            if self._get_source() == SourceType.calibration:

                # updates for multiple plot y-axes
                updateViews1()
                self._plt0.vb.sigResized.connect(updateViews1)

                calibration_readFREQ = np.arange(len(self.worker.get_value1_buffer())) * (
                    Constants.calib_fStep) + Constants.calibration_frequency_start
                self._plt0.plot(x=calibration_readFREQ, y=self.worker.get_value1_buffer(), pen=Constants.plot_colors[0])
                self._plt1.addItem(pg.PlotCurveItem(x=calibration_readFREQ, y=self.worker.get_value2_buffer(),
                                                    pen=Constants.plot_colors[1]))

            # SINGLE
            # -----------------------------------------------------------------
            elif self._get_source() == SourceType.serial:

                # AMPLITUDE and PHASE
                # -------------------------------------------------------------

                # updates for multiple plot y-axes
                updateViews1()
                self._plt0.vb.sigResized.connect(updateViews1)

                self._plt0.plot(x=self._readFREQ, y=self.worker.get_value1_buffer(), pen=Constants.plot_colors[0])
                self._plt1.addItem(
                    pg.PlotCurveItem(x=self._readFREQ, y=self.worker.get_value2_buffer(), pen=Constants.plot_colors[1]))

                #  TODO set the legend in single mode
                overtone_selected = self._overtones_number_all - self.ui.cBox_Speed.currentIndex() - 1

                y_freq = self.worker.get_d1_buffer()
                y_diss = self.worker.get_d2_buffer()

                # VER 0.1.2
                # get y_freq and y_dissipation max value
                y_freq_single_max = np.nanmax(y_freq)
                y_freq_single_min = np.nanmin(y_freq)
                y_diss_single_max = np.nanmax(y_diss)
                y_diss_single_min = np.nanmin(y_diss)

                # FREQUENCY and DISSIPATION
                # --------------------------------------------------------------
                self._plt2.clear()
                time_x = self.worker.get_t1_buffer()
                self._pltD.clear()
                time_x_D = self.worker.get_t3_buffer()

                # VER 0.1.2
                # check if t y-axis limit is not a nan
                if not (np.isnan(y_freq_single_max)):
                    try:
                        self._plt2.setYRange(y_freq_min - 100, y_freq_max + 100, padding=0)
                        self._pltD.setYRange(y_diss_min - 0.000001, y_diss_max + 0.000001, padding=0)
                    except:
                        pass

                self._plt2.plot(x=time_x, y=y_freq, pen=pg.mkPen(color=Constants.plot_color_multi[overtone_selected],
                                                                 width=Constants.plot_line_width))
                self._pltD.plot(x=time_x_D, y=y_diss, pen=pg.mkPen(color=Constants.plot_color_multi[overtone_selected],
                                                                   width=Constants.plot_line_width))

                # update frequency and dissipatuon indicator
                self._update_indicator_F_single(overtone_selected, y_freq)
                self._update_indicator_D_single(overtone_selected, y_diss)

                # TEMPERATURE
                # -----------------------------------------------------------------

                self._plt4.clear()
                # do not autoscale y
                self._plt4.enableAutoRange(axis='y', enable=True)
                # set temperature y range
                # VER 0.1.2
                # change the Temperature Y-range to 5 - 45 °C
                self._plt4.setYRange(5, 45, padding=0)

                # get temperature buffer
                y_temperature = self.worker.get_d3_buffer()
                self._plt4.plot(x=time_x_D, y=y_temperature, pen=Constants.plot_colors[4])

                # set temperature current value
                label_indicator_temperature = float("{0:.2f}".format(y_temperature[0]))
                self.ui.indicator_temperature.setText(str(label_indicator_temperature))

            # MULTI
            # -----------------------------------------------------------------
            elif self._get_source() == SourceType.multiscan:
                '''
               # AMPLITUDE
               # -------------------------------------------------------------
               # Loads frequencies from calibration file
               peaks_mag = self.load_frequencies_file()
               # get sweep frequency range for each overtone
               x_sweep_axis = self.worker.get_frequency_range_multi(Constants.argument_default_samples, self._overtone_number)
               # center the sweep frequency axis around each overtone frequency peak
               x_sweep_axis = x_sweep_axis - peaks_mag[self._overtone_number]

               # 10M
               if (peaks_mag[0] >9e+06 and peaks_mag[0]<11e+06):
                   # set XY constant in multiscan mode
                   self._plt0.setXRange(-(Constants.L10_5th_overtone + 1000), Constants.R10_5th_overtone + 1000, padding = 0)
                   self._plt0.setYRange(-5, 30, padding = 0)

               # TODO 5M
               if (peaks_mag[0] >4e+06 and peaks_mag[0]<6e+06):
                   self._plt0.setXRange( -(Constants.L5_7th_overtone + 1000), Constants.R5_7th_overtone + 1000, padding = 0 )

               # updates for multiple plot y-axes
               updateViews_multi()
               self._plt0.vb.sigResized.connect(updateViews_multi)

               # TODO clear plot at fundamental sweep
               if self._overtone_number == 0:
                   self._plt0.clear()

               # plot
               if self.scan_selector[self._overtone_number] == True:
                   self._plt0.plot ( x = x_sweep_axis, y = self.worker.get_value1_buffer(), pen = Constants.plot_color_multi[self._overtone_number] )
               '''

                # AMPLITUDE
                # -------------------------------------------------------------
                # Loads frequencies from calibration file
                peaks_mag = self.load_frequencies_file()

                # VER 0.1.2
                # =============================================================================
                #                # 10M
                #                if (peaks_mag[0] >9e+06 and peaks_mag[0]<11e+06):
                #                    # set XY constant in multiscan mode
                #                    self._plt0.setXRange(-(Constants.L10_5th_overtone + 1000), Constants.R10_5th_overtone + 1000, padding = 0)
                #                    # self._plt0.setYRange(-5, 35, padding = 0)
                #                # TODO 5M
                #                if (peaks_mag[0] >4e+06 and peaks_mag[0]<6e+06):
                #                    # set XY constant in multiscan mode
                #                    self._plt0.setXRange( -(Constants.L5_7th_overtone + 1000), Constants.R5_7th_overtone + 1000, padding = 0 )
                #                    # self._plt0.setYRange(-5, 30, padding = 0)
                # =============================================================================

                # updates for multiple plot y-axes
                updateViews_multi()
                self._plt0.vb.sigResized.connect(updateViews_multi)

                self._plt0.clear()

                # loop on the lines
                # TODO introduce this control because I've some trouble in reset buffer ???
                if self._ser_control > Constants.environment:
                    for idx in range(self._overtones_number_all):
                        # get and scale frequency axis
                        x_sweep_axis = self.worker.get_F_Sweep_values_buffer(idx) - peaks_mag[idx]
                        # get amplitude axis
                        y_sweep_axis = self.worker.get_A_values_buffer(idx)
                        # plot sweep
                        if self.scan_selector[idx] == True and isinstance(x_sweep_axis, np.ndarray):
                            self._plt0.plot(x=x_sweep_axis, y=y_sweep_axis, pen=Constants.plot_color_multi[idx])

                # FREQUENCY and DISSIPATION
                # --------------------------------------------------------------
                # clear plot
                self._plt2.clear()
                self._pltD.clear()

                # loop on the lines
                for idx in range(self._overtones_number_all):
                    if self.scan_selector[idx] == True:
                        # get time axis
                        time_axis_new = self.worker.get_time_values_buffer(idx)

                        # get y frequency and dissipation axis
                        y_freq = self.worker.get_F_values_buffer(idx)
                        y_diss = self.worker.get_D_values_buffer(idx)
                        # frequency
                        self._plt2.plot(x=time_axis_new, y=y_freq, pen=pg.mkPen(color=Constants.plot_color_multi[idx],
                                                                                width=Constants.plot_line_width))
                        # dissipation
                        self._pltD.plot(x=time_axis_new, y=y_diss, pen=Constants.plot_color_multi[idx],
                                        width=Constants.plot_line_width)

                        self._update_indicator_F(idx, y_freq)
                        self._update_indicator_D(idx, y_diss)

                    else:
                        dummy = [0]
                        self._update_indicator_F(idx, dummy)
                        self._update_indicator_D(idx, dummy)

                # TEMPERATURE
                # --------------------------------------------------------------
                self._plt4.clear()
                # do not autoscale y
                self._plt4.enableAutoRange(axis='y', enable=True)
                # VER 0.1.2
                # change the Temperature Y-range to 5 - 45 °C
                self._plt4.setYRange(5, 45, padding=0)

                # get temperarre buffer
                y_temperature = self.worker.get_d3_buffer()

                self._plt4.plot(x=self.worker.get_t3_buffer(), y=y_temperature, pen=Constants.plot_colors[4])

                # set temperature current value
                label_indicator_temperature = float("{0:.2f}".format(y_temperature[0]))
                self.ui.indicator_temperature.setText(str(label_indicator_temperature))


        """
        # INSA : INJECTION: When the balance begins monitoring, injection is started
        # --------------------------------------------------------------------
        # TODO : Emergency Stop Injection button
        pressure = self.ui.spinBox_Pressure.value()
        injection_time = self.ui.spinBox_Injection_Time.value()
        if labelbar == 'Monitoring!' and Injection_Bool:
            # pressure = self.ui.spinBox_Pressure.value()
            # injection_time = self.ui.spinBox_Injection_Time.value()

            ## Initialize the session
            # This step is optional, if not called session will be automatically created
            fgt_init()

            # Set pressure to the given pressure on first pressure channel of the list
            # mbar is the default unit at initialization
            fgt_set_pressure(0, pressure)

            # Wait injection_time (in seconds) before setting pressure to 0
            # injection time is in ms and time.sleep is in sec
            #time.sleep(injection_time / 1000)

            ## Close the session
            # Set pressure to 0 before closing
            #fgt_set_pressure(0, 0)

            Injection_Bool = False
            i = 0
            i += 1
            print(i)
            y = True
            t_end = time.time() + injection_time / 1000
            while time.time() < t_end:
                if Stop_Injection_Bool and y:
                    fgt_set_pressure(0, 0)
                    y = False
            fgt_set_pressure(0, 0)

            fgt_close()


#            ## Initialize the session
#            # This step is optional, if not called session will be automatically created
#            fgt_init()
#
#            # Set pressure to the given pressure on first pressure channel of the list
#            # mbar is the default unit at initialization
#            fgt_set_pressure(0, pressure)
#
#            while not exit.is_set():
#                # Wait injection_time (in seconds) before setting pressure to 0
#                # injection time is in ms and time.sleep is in sec
#                exit.wait(injection_time / 1000)
#
#            ## Close the session
#            # Set pressure to 0 before closing
#            fgt_set_pressure(0, 0)
#            fgt_close()
#
#            Injection_Bool = False
        """

    def _update_indicator_F(self, index, value):

        label = float("{0:.1f}".format(value[0]))

        if (index == 0):
            if (self.scan_selector[index] == True):
                self.ui.F0.setText(str(label))
            else:
                self.ui.F0.setText("nan")
        elif (index == 1):
            if (self.scan_selector[index] == True):
                self.ui.F3.setText(str(label))
            else:
                self.ui.F3.setText("nan")
        elif (index == 2):
            if (self.scan_selector[index] == True):
                self.ui.F5.setText(str(label))
            else:
                self.ui.F5.setText("nan")
        elif (index == 3):
            if (self.scan_selector[index] == True):
                self.ui.F7.setText(str(label))
            else:
                self.ui.F7.setText("nan")
        elif (index == 4):
            if (self.scan_selector[index] == True):
                self.ui.F9.setText(str(label))
            else:
                self.ui.F9.setText("nan")

    def _update_indicator_F_single(self, index, value):

        label = float("{0:.1f}".format(value[0]))

        if (index == 0):
            self.ui.F0.setText(str(label))
        else:
            self.ui.F0.setText("nan")

        if (index == 1):
            self.ui.F3.setText(str(label))
        else:
            self.ui.F3.setText("nan")

        if (index == 2):
            self.ui.F5.setText(str(label))
        else:
            self.ui.F5.setText("nan")

        if (index == 3):
            self.ui.F7.setText(str(label))
        else:
            self.ui.F7.setText("nan")

        if (index == 4):
            self.ui.F9.setText(str(label))
        else:
            self.ui.F9.setText("nan")

    def _update_indicator_D(self, index, value):
        value_multiplied = value[0] * 1e6
        label = float("{0:.3f}".format(value_multiplied))

        if (index == 0):
            if (self.scan_selector[index] == True):
                self.ui.D0.setText(str(label))
            else:
                self.ui.D0.setText("nan")
        elif (index == 1):
            if (self.scan_selector[index] == True):
                self.ui.D3.setText(str(label))
            else:
                self.ui.D3.setText("nan")
        elif (index == 2):
            if (self.scan_selector[index] == True):
                self.ui.D5.setText(str(label))
            else:
                self.ui.D5.setText("nan")
        elif (index == 3):
            if (self.scan_selector[index] == True):
                self.ui.D7.setText(str(label))
            else:
                self.ui.D7.setText("nan")
        elif (index == 4):
            if (self.scan_selector[index] == True):
                self.ui.D9.setText(str(label))
            else:
                self.ui.D9.setText("nan")

    def _update_indicator_D_single(self, index, value):
        value_multiplied = value[0] * 1e6
        label = float("{0:.3f}".format(value_multiplied))

        if (index == 0):
            self.ui.D0.setText(str(label))
        else:
            self.ui.D0.setText("nan")

        if (index == 1):
            self.ui.D3.setText(str(label))
        else:
            self.ui.D3.setText("nan")
        if (index == 2):
            self.ui.D5.setText(str(label))
        else:
            self.ui.D5.setText("nan")

        if (index == 3):
            self.ui.D7.setText(str(label))
        else:
            self.ui.D7.setText("nan")

        if (index == 4):
            self.ui.D9.setText(str(label))
        else:
            self.ui.D9.setText("nan")

    def _update_scan_selector(self):
        self.scan_selector[0] = self.ui.radioBtn_F0.isChecked()
        self.scan_selector[1] = self.ui.radioBtn_F3.isChecked()
        self.scan_selector[2] = self.ui.radioBtn_F5.isChecked()
        self.scan_selector[3] = self.ui.radioBtn_F7.isChecked()
        self.scan_selector[4] = self.ui.radioBtn_F9.isChecked()

    ###########################################################################
    # Updates the source and depending boxes on change
    ###########################################################################
    def _source_changed(self):

        # It is connected to the indexValueChanged signal of the Source ComboBox.

        # single frequency measurement
        if self._get_source() == SourceType.serial:
            print(TAG, "Scanning the source: {}".format(Constants.app_sources[1]))  # self._get_source().name
            Log.i(TAG, "Scanning the source: {}".format(Constants.app_sources[1]))

            # show - hide overtone radio button selector
            self._Overtone_radioBtn_isEnabled(False)

        # calibration
        elif self._get_source() == SourceType.calibration:
            print(TAG, "Scanning the source: {}".format(Constants.app_sources[0]))  # self._get_source().name
            Log.i(TAG, "Scanning the source: {}".format(Constants.app_sources[0]))

            # show - hide overtone radio button selector
            self._Overtone_radioBtn_isEnabled(False)


        # multi frequency measurement
        elif self._get_source() == SourceType.multiscan:
            print(TAG, "Scanning the source: {}".format(Constants.app_sources[2]))  # self._get_source().name
            Log.i(TAG, "Scanning the source: {}".format(Constants.app_sources[2]))

            # show - hide overtone radio button selector
            self._Overtone_radioBtn_isEnabled(True)

            # init radio button
            self.ui.radioBtn_F0.setChecked(True)
            self.ui.radioBtn_F3.setChecked(True)
            self.ui.radioBtn_F5.setChecked(True)
            self.ui.radioBtn_F7.setChecked(True)
            self.ui.radioBtn_F9.setChecked(True)

            self._update_scan_selector()

        '''
        # Clears boxes before adding new
        self.ControlsWin.ui1.cBox_Port.clear()
        self.ControlsWin.ui1.cBox_Speed.clear()
        '''
        '''
        TODO 2m
        '''
        # Clears boxes before adding new
        self.ui.cBox_Port.clear()
        self.ui.cBox_Speed.clear()

        # Gets the current source type
        source = self._get_source()
        ports = self.worker.get_source_ports(source)
        speeds = self.worker.get_source_speeds(source)

        '''
        if ports is not None:
            self.ControlsWin.ui1.cBox_Port.addItems(ports)
        if speeds is not None:
            self.ControlsWin.ui1.cBox_Speed.addItems(speeds)
        if self._get_source() == SourceType.serial:
            self.ControlsWin.ui1.cBox_Speed.setCurrentIndex(len(speeds) - 1)
        '''
        '''
        TODO 2m
        '''
        # set COM port
        if ports is not None:
            self.ui.cBox_Port.addItems(ports)
        if speeds is not None:
            self.ui.cBox_Speed.addItems(speeds)

        # populates the drop-down menu with detected freqeuencies
        if self._get_source() == SourceType.serial:
            self.ui.cBox_Speed.setCurrentIndex(len(speeds) - 1)

    ###########################################################################
    # Gets the current source type
    ###########################################################################
    def _get_source(self):

        '''
        #:rtype: SourceType.
        return SourceType(self.ControlsWin.ui1.cBox_Source.currentIndex())
        '''
        '''
        TODO 2m
        '''
        return SourceType(self.ui.cBox_Source.currentIndex())

    # LOAD FREQUENCENCIES FILE
    @staticmethod
    def load_frequencies_file():
        data = loadtxt(Constants.cvs_peakfrequencies_path)
        peaks_mag = data[:, 0]
        # peaks_phase = data[:,1] #unused at the moment
        return peaks_mag

    ###########################################################################
    # Cleans history plot
    ###########################################################################
    def clear(self):
        support = self.worker.get_d1_buffer()
        if support.any:
            if str(support[0]) != 'nan':
                print(TAG, "All Plots Cleared!", end='\r')
                self._update_sample_size()
                self._plt2.clear()
                self._pltD.clear()
                self._plt4.clear()

                # VER 0.1.2
                # clear amplitude and phase sweep plot
                self._plt0.clear()

    ###########################################################################
    # Reference set/reset
    ###########################################################################
    # =============================================================================
    #     def reference(self):
    #         import numpy as np
    #         #import sys
    #         support=self.worker.get_d1_buffer()
    #         if support.any:
    #             if str(support[0])!='nan':
    #                 ref_vector1 = [c for c in self.worker.get_d1_buffer() if ~np.isnan(c)]
    #                 ref_vector2 = [c for c in self.worker.get_d2_buffer() if ~np.isnan(c)]
    #                 self._reference_value_frequency = ref_vector1[0]
    #                 self._reference_value_dissipation = ref_vector2[0]
    #                 #sys.stdout.write("\033[K") #clear line
    #                 if self._reference_flag:
    #                     self._reference_flag = False
    #                     print(TAG, "Reference reset!   ", end='\r')
    #                     self._labelref1 = "not set"
    #                     self._labelref2 = "not set"
    #                 else:
    #                     self._reference_flag = True
    #                     d1=float("{0:.2f}".format(self._reference_value_frequency))
    #                     d2=float("{0:.4f}".format(self._reference_value_dissipation*1e6))
    #                     self._labelref1 = str(d1)+ "Hz"
    #                     self._labelref2 = str(d2)+ "e-06"
    #                     print(TAG, "Reference set!     ", end='\r')
    #                     # TODO minor changes: it calculates (in a unelegant way) the frequency y - range axis
    #                     self._vector_reference_frequency[:] = [s - self._reference_value_frequency for s in self._readFREQ]
    #                     # TODO minor changes: it calculates (in a unelegant way) the dissipation y - range axis
    #                     xs = np.array(np.linspace(0, ((self._readFREQ[-1]-self._readFREQ[0])/self._readFREQ[0]), len(self._readFREQ)))
    #                     self._vector_reference_dissipation = xs-self._reference_value_dissipation
    # =============================================================================

    ###########################################################################
    # Set reference
    ###########################################################################
    def reference(self):
        import numpy as np

        # SINGLE
        # ---------------------------------------------------------------------
        if self._get_source() == SourceType.serial:
            support = self.worker.get_d1_buffer()
            if support.any:
                if str(support[0]) != 'nan':
                    ref_vector1 = [c for c in self.worker.get_d1_buffer() if ~np.isnan(c)]
                    ref_vector2 = [c for c in self.worker.get_d2_buffer() if ~np.isnan(c)]

                    # get frequency reference value
                    self._reference_value_frequency = ref_vector1[0]
                    # get dissipation reference value
                    self._reference_value_dissipation = ref_vector2[0]

                    # =============================================================================
                    #                     if self._reference_flag:
                    #                         self._reference_flag = False
                    #                         print(TAG, "Reference reset!   ", end='\r')
                    #                         self._labelref1 = "not set"
                    #                         self._labelref2 = "not set"
                    #
                    #                         # clear all
                    #                         self.clear()
                    # =============================================================================

                    # VER 0.1.2
                    # changed the function set reference in single mode

                    # set flag reference true
                    self._reference_flag = True

                    # get frequency and dissipation reference value
                    d1 = float("{0:.2f}".format(self._reference_value_frequency))
                    d2 = float("{0:.4f}".format(self._reference_value_dissipation * 1e6))

                    self._labelref1 = str(d1) + "Hz"
                    self._labelref2 = str(d2) + "e-06"

                    print(TAG, "Reference set!     ", end='\r')
                    # TODO minor changes: it calculates (in a unelegant way) the frequency y - range axis
                    self._vector_reference_frequency[:] = [s - self._reference_value_frequency for s in self._readFREQ]
                    # TODO minor changes: it calculates (in a unelegant way) the dissipation y - range axis
                    xs = np.array(np.linspace(0, ((self._readFREQ[-1] - self._readFREQ[0]) / self._readFREQ[0]),
                                              len(self._readFREQ)))
                    self._vector_reference_dissipation = xs - self._reference_value_dissipation

                    # clear all
                    # self.clear()
        # MULTI
        # ---------------------------------------------------------------------
        elif self._get_source() == SourceType.multiscan:

            # =============================================================================
            #             if self._reference_flag:
            #                 print(TAG, "Reference reset!   ", end='\r')
            #                 self._reference_flag = False
            #                 # clear all
            #                 self.clear()
            # =============================================================================

            # VER 0.1.2
            # changed the function set reference in multiscan mode

            print(TAG, "Reference set!     ", end='\r')

            # set flag reference true
            self._reference_flag = True

            # get current frequencies overtone array
            for idx in range(self._overtones_number_all):
                frequency_reference = self.worker.get_F_values_buffer(idx)
                self._reference_value_frequency_array[idx] = frequency_reference[0]

            # get current dissipation overtone array
            for idx in range(self._overtones_number_all):
                dissipation_reference = self.worker.get_D_values_buffer(idx)
                self._reference_value_dissipation_array[idx] = dissipation_reference[0]

    ###########################################################################
    # Reset reference
    ###########################################################################
    # VER 0.1.2
    # add reset reference button / function
    def reference_not(self):
        print(TAG, "Reference reset!   ", end='\r')
        self._reference_flag = False
        # clear all
        self.clear()

        self._plt2.enableAutoRange(axis='y', enable=True)
        self._pltD.enableAutoRange(axis='y', enable=True)

    ###########################################################################
    # Checking internet connection
    ###########################################################################
    def internet_on(self):
        from urllib.request import urlopen
        try:
            url = "https://openqcm.com/shared/news.html"
            urlopen(url, timeout=10)
            return True
        except:
            return False

    ########################################################################################################
    # Gets information from openQCM webpage and enables download button if new version software is available
    ########################################################################################################
    def get_web_info(self):
        import pandas as pd
        # check if an Internet connection is active
        self._internet_connected = self.internet_on()
        # Get latest info from openQCM webpage
        c_types = {
            '1': '1',
            '2': '2',
            '3': '3', }
        r_types = {
            '1': 'A',
            '2': 'B',
            '3': 'C', }
        if self._internet_connected:
            color = '#00c600'
            labelweb2 = 'ONLINE'
            print(TAG, 'Checking your internet connection {} '.format(labelweb2))
            tables = pd.read_html('https://openqcm.com/shared/news.html', index_col=0, header=0, match='1')
            df = tables[0]
            # create empty list of string
            self._webinfo = ["" for x in range(len(df.columns) * len(df.index))]  # len(df.columns)*len(df.index)=9
            # row acess mode to Pandas dataframe
            k = 0
            for j in [1, 2, 3]:
                for i in [1, 2, 3]:
                    self._webinfo[k] = str(df.loc[r_types[str(j)], c_types[str(i)]])
                    k += 1
            # check for update
            if self._webinfo[0] == Constants.app_version:
                labelweb3 = 'last version installed!'
            else:
                labelweb3 = 'version {} available!'.format(self._webinfo[0])
                '''
              self.InfoWin.ui3.pButton_Download.setEnabled(True)
              '''
                '''
              TODO 2m set enabled infowin buton download
              '''
        else:
            color = '#ff0000'
            labelweb2 = 'OFFLINE'
            labelweb3 = 'Offline, unable to check'
            print(TAG, 'Checking your internet connection {} '.format(labelweb2))

        '''
        self.InfoWin.ui3.lweb2.setText("<font color=#0000ff > Checking your internet connection &nbsp;&nbsp;&nbsp;&nbsp;</font><font size=4 color={}>{}</font>".format(color,labelweb2))
        self.InfoWin.ui3.lweb3.setText("<font color=#0000ff > Software update status </font>" + labelweb3)
        '''
        '''
        TODO 2m set infowin internet connection adn software update
        '''

    ###########################################################################
    # Opens webpage for download
    ###########################################################################
    def start_download(self):
        import webbrowser
        url_download = 'https://openqcm.com/shared/q-1/openQCM_Q-1_py_v{}.zip '.format(self._webinfo[0])
        webbrowser.open(url_download)
