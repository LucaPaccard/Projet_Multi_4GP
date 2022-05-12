
from openQCM.ui.mainWindow_ui import Ui_Controls, Ui_Info, Ui_Plots
#from openQCM.ui.ui_controls import Ui_Controls
#from openQCM.ui.ui_info import Ui_Info
#from openQCM.ui.ui_plots import Ui_Plots

from pyqtgraph import AxisItem
import pyqtgraph as pg
from PyQt5 import QtCore, QtGui
from openQCM.core.worker import Worker
from openQCM.core.constants import Constants, SourceType, DateAxis, NonScientificAxis
from openQCM.ui.popUp import PopUp
from openQCM.common.logger import Logger as Log
import numpy as np
import sys

TAG = ""#"[MainWindow]"

##########################################################################################
# Package that handles the UIs elements and connects to worker service to execute processes
##########################################################################################

class PlotsWindow(QtGui.QMainWindow):
     
    def __init__(self, samples=Constants.argument_default_samples):
        super().__init__()  
        self.ui2 = Ui_Plots()
        self.ui2.setupUi(self)
    '''
    def closeEvent(self, event):
        #print(" Exit Real-Time Plot GUI")
        res =PopUp.question(self, Constants.app_title, "Are you sure you want to quit openQCM application now?")
        if res: 
           #self.close()
           QtGui.QApplication.quit()
        else:
           event.ignore()#pass
        #sys.exit(0)
    '''
    def closeEvent(self, event):
        #print(" Exit Real-Time Plot GUI")
           event.ignore()

    
#------------------------------------------------------------------------------
class InfoWindow(QtGui.QMainWindow):    
    
    def __init__(self, samples=Constants.argument_default_samples):
        super().__init__()
        self.ui3 = Ui_Info()
        self.ui3.setupUi(self)
    '''    
    def closeEvent(self, event):
        #print(" Exit Information GUI") 
        res =PopUp.question(self, Constants.app_title, "Are you sure you want to quit openQCM application now?")
        if res: 
           #self.close()
           QtGui.QApplication.quit()
        else:
           event.ignore()
    '''
    def closeEvent(self, event):
        #print(" Exit Real-Time Plot GUI")
           event.ignore()    

#------------------------------------------------------------------------------
class ControlsWindow(QtGui.QMainWindow):    
    
    def __init__(self, samples=Constants.argument_default_samples):
        super().__init__()
        self.ui1 = Ui_Controls()
        self.ui1.setupUi(self)
        
    def closeEvent(self, event):
        #print(" Exit Setup/Control GUI")
        res =PopUp.question(self, Constants.app_title, "Are you sure you want to quit openQCM application now?")
        if res: 
           #self.close()
           QtGui.QApplication.quit()
        else:
           event.ignore()

#------------------------------------------------------------------------------

class MainWindow(QtGui.QMainWindow):
    
    ###########################################################################
    # Initializes methods, values and sets the UI
    ###########################################################################   
    def __init__(self, samples=Constants.argument_default_samples):
        
        #:param samples: Default samples shown in the plot :type samples: int.
        # to be always placed at the beginning, initializes some important methods
        QtGui.QMainWindow.__init__(self)
        
        # Sets up the user interface from the generated Python script using Qt Designer
        # Instantiates Ui classes

        # Calls setupUi method of created instance
        self.ControlsWin = ControlsWindow()
        self.ControlsWin.move(530, 568)
        self.PlotsWin = PlotsWindow()
        self.PlotsWin.move(530,10) #GUI position (x,y) on the screen 
        self.InfoWin = InfoWindow()
        self.InfoWin.move(1242, 10)
        self.ControlsWin.show()
        self.PlotsWin.show()
        self.InfoWin.show()

        # Shared variables, initial values
        self._plt0 = None
        self._plt1 = None
        self._plt2 = None
        self._plt3 = None
        self._plt4 = None
        self._timer_plot = None
        self._readFREQ = None
        self._QCS_installed = None
        self._ser_control = 0
        self._ser_error1 = 0
        self._ser_error2 = 0
        self._ser_err_usb= 0
        
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
        self.ControlsWin.ui1.cBox_Source.addItems(Constants.app_sources)

        # Configures specific elements of the PyQtGraph plots
        self._configure_plot()
       
        # Configures specific elements of the QTimers
        self._configure_timers()
        
        # Configures the connections between signals and UI elements
        self._configure_signals()

        # Populates combo box for serial ports
        self._source_changed()
        self.ControlsWin.ui1.cBox_Source.setCurrentIndex(SourceType.serial.value)
        self.ControlsWin.ui1.sBox_Samples.setValue(samples)  #samples
        
        # enable ui
        self._enable_ui(True)
        ###################################################################################################################################
        self.get_web_info()
        # Gets the QCS installed on the device (not used now)
        # self._QCS_installed = PopUp.question_QCM(self, Constants.app_title, "Please choose the Quartz Crystal Resonator installed on the openQCM-1 Device (default 5MHz if exit)")

    ###########################################################################
    # Starts the acquisition of the selected serial port
    ###########################################################################      
    def start(self):
        
        import os
        os.system('cls' if os.name == 'nt' else 'clear')
        
        # This function is connected to the clicked signal of the Start button.
        #print("")
        print(TAG, 'Clicked START')
        Log.i(TAG, "Clicked START")
        # Instantiates process
        self.worker = Worker(QCS_on = self._QCS_installed,
                             port = self.ControlsWin.ui1.cBox_Port.currentText(),
                             speed = self.ControlsWin.ui1.cBox_Speed.currentText(),
                             samples = self.ControlsWin.ui1.sBox_Samples.value(),
                             source = self._get_source(),
                             export_enabled = self.ControlsWin.ui1.chBox_export.isChecked())     
        
        if self.worker.start():
            # Gets frequency range 
            self._readFREQ = self.worker.get_frequency_range()
            # Duplicate frequencies
            self._reference_flag = False
            self._vector_reference_frequency = list(self._readFREQ)
            self._reference_value_frequency = 0
            self._reference_value_dissipation = 0
            self._labelref1 = "not set"
            self._labelref2 = "not set"
            # progressbar variables
            self._completed=0
            self._ser_control = 0
            # error variables
            self._ser_error1 = 0
            self._ser_error2 = 0
            self._ser_err_usb= 0
            ##### other useful location #########
            #self.get_web_info()
            #####

            if self._get_source() == SourceType.serial:
                overtones_number = len(self.worker.get_source_speeds(SourceType.serial))
                if overtones_number==5:
                   label_quartz = "@5MHz_QCM"
                elif overtones_number==3:
                   label_quartz = "@10MHz_QCM"
                self.InfoWin.ui3.info1a.setText("<font color=#0000ff > Device Setup </font>" + label_quartz)
                label11= "Measurement openQCM Q-1"
                self.InfoWin.ui3.info11.setText("<font color=#0000ff > Operation Mode </font>" + label11)
                self._overtone_name,self._overtone_value, self._fStep = self.worker.get_overtone()
                label6= str(int(self._overtone_value))+" Hz"
                self.InfoWin.ui3.info6.setText("<font color=#0000ff > Frequency Value </font>" + label6)
                label2= str(self._overtone_name)# +" "+ str(self._overtone_value)+"Hz"
                self.InfoWin.ui3.info2.setText("<font color=#0000ff > Selected Frequency </font>" + label2)
                label3= str(int(self._readFREQ[0]))+" Hz"
                self.InfoWin.ui3.info3.setText("<font color=#0000ff > Start Frequency </font>" + label3)
                label4= str(int(self._readFREQ[-1]))+" Hz"
                self.InfoWin.ui3.info4.setText("<font color=#0000ff > Stop Frequency </font>" + label4)
                label4a= str(int(self._readFREQ[-1]-self._readFREQ[0]))+" Hz"
                self.InfoWin.ui3.info4a.setText("<font color=#0000ff > Frequency Range </font>" + label4a)                      
                label5= str(int(self._fStep))+" Hz"
                self.InfoWin.ui3.info5.setText("<font color=#0000ff > Sample Rate </font>" + label5) 
                label7= str(Constants.argument_default_samples-1)
                self.InfoWin.ui3.info7.setText("<font color=#0000ff > Sample Number </font>" + label7)
                                     
            elif self._get_source() == SourceType.calibration:
                label_quartz = self.ControlsWin.ui1.cBox_Speed.currentText()
                self.InfoWin.ui3.info1a.setText("<font color=#0000ff > Device Setup </font>" + label_quartz)
                label11= "Calibration openQCM Q-1"
                self.InfoWin.ui3.info11.setText("<font color=#0000ff > Operation Mode </font>" + label11)
                label6= "Overall Frequency Range"
                self.InfoWin.ui3.info6.setText("<font color=#0000ff > Frequency Value </font>" + label6)
                label2= "Overall Frequency Range"
                self.InfoWin.ui3.info2.setText("<font color=#0000ff > Selected Frequency </font>" + label2)
                label3= str(Constants.calibration_frequency_start)+" Hz"
                self.InfoWin.ui3.info3.setText("<font color=#0000ff > Start Frequency </font>" + label3)
                label4= str(Constants.calibration_frequency_stop)+" Hz"
                self.InfoWin.ui3.info4.setText("<font color=#0000ff > Stop Frequency </font>" + label4) 
                label4a= str(int(Constants.calibration_frequency_stop - Constants.calibration_frequency_start))+" Hz"
                self.InfoWin.ui3.info4a.setText("<font color=#0000ff > Frequency Range </font>" + label4a)                        
                label5= str(int(Constants.calibration_fStep))+" Hz"
                self.InfoWin.ui3.info5.setText("<font color=#0000ff > Sample Rate </font>" + label5)  
                label7= str(Constants.calibration_default_samples-1)
                self.InfoWin.ui3.info7.setText("<font color=#0000ff > Sample Number </font>" + label7)  
            #
            self._timer_plot.start(Constants.plot_update_ms)
            self._enable_ui(False)
            self.ControlsWin.ui1.sBox_Samples.setEnabled(False) #insert
            
            if self._get_source() == SourceType.calibration:
               self.ControlsWin.ui1.pButton_Clear.setEnabled(False) #insert
               self.ControlsWin.ui1.pButton_Reference.setEnabled(False) #insert
        else:
            print(TAG, "Warning: port is not available!")
            Log.i(TAG, "Warning: port is not available")
            PopUp.warning(self, Constants.app_title, "Warning: Selected Port [{}] is not available!".format(self.ControlsWin.ui1.cBox_Port.currentText()))

        
    ###########################################################################
    # Stops the acquisition of the selected serial port
    ###########################################################################    
    def stop(self):
        
        # This function is connected to the clicked signal of the Stop button.
        self.ControlsWin.ui1.infostatus.setStyleSheet('background: white; padding: 1px; border: 1px solid #cccccc')
        self.ControlsWin.ui1.infostatus.setText("<font color=#000000 > Program Status Stanby</font>")
        self.ControlsWin.ui1.infobar.setText("<font color=#0000ff > Infobar </font>")
        self.InfoWin.ui3.inforef1.setText("<font color=#0000ff > Ref. Frequency </font>")               
        self.InfoWin.ui3.inforef2.setText("<font color=#0000ff > Ref. Dissipation </font>")
        print("")
        print(TAG, "Clicked STOP")
        Log.i(TAG, "Clicked STOP")
        self._timer_plot.stop()
        self._enable_ui(True)
        self.worker.stop()
    
    ###########################################################################
    # Overrides the QTCloseEvent,is connected to the close button of the window
    ###########################################################################      
    def closeEvent(self, evnt):
        
        #:param evnt: QT evnt.
        if self.worker.is_running():
          print(TAG, 'Window closed without stopping the capture, application will stop...')
          Log.i(TAG, "Window closed without stopping the capture, application will stop...")
          self.stop()
          #self.ControlsWin.close()
          #self.PlotsWin.close()
          #self.InfoWin.close()
          #evnt.accept()
    
          
    ###########################################################################
    # Enables or disables the UI elements of the window.
    ########################################################################### 
    def _enable_ui(self, enabled):
        
        #:param enabled: The value to be set for the UI elements :type enabled: bool
        self.ControlsWin.ui1.cBox_Port.setEnabled(enabled)
        self.ControlsWin.ui1.cBox_Speed.setEnabled(enabled)
        self.ControlsWin.ui1.pButton_Start.setEnabled(enabled)
        self.ControlsWin.ui1.chBox_export.setEnabled(enabled)
        self.ControlsWin.ui1.cBox_Source.setEnabled(enabled)
        self.ControlsWin.ui1.pButton_Stop.setEnabled(not enabled)
        self.ControlsWin.ui1.sBox_Samples.setEnabled(not enabled) #insert
        self.ControlsWin.ui1.pButton_Clear.setEnabled(not enabled)
        self.ControlsWin.ui1.pButton_Reference.setEnabled(not enabled)


    ###########################################################################
    # Configures specific elements of the PyQtGraph plots.
    ########################################################################### 
    def _configure_plot(self):
        
        #----------------------------------------------------------------------
        # set background color 
        self.PlotsWin.ui2.plt.setBackground(background='#0e1c20')
        self.PlotsWin.ui2.pltB.setBackground(background='#0e1c20')                          
        
        #----------------------------------------------------------------------
        # defines the graph title
        title1 = "Real-Time Plot: Amplitude / Phase"
        title2 = "Real-Time Plot: Resonance Frequency / Dissipation"
        title3 = "Real-Time Plot: Temperature"
        #--------------------------------------------------------------------------------------------------------------
        # Configures elements of the PyQtGraph plots: amplitude
        self.PlotsWin.ui2.plt.setAntialiasing(True)
        self.PlotsWin.ui2.pltB.setAntialiasing(True)
        self._plt0 = self.PlotsWin.ui2.plt.addPlot(row=0, col=1, title= title1, **{'font-size':'12pt'})
        self._plt0.showGrid(x=True, y=True)
        self._plt0.setLabel('bottom', 'Frequency', units='Hz')
        self._plt0.setLabel('left', 'Amplitude', units='dB', color=Constants.plot_colors[0], **{'font-size':'10pt'})
        
        '''
        # Configures elements of the PyQtGraph plots: phase
        self._plt1 = self.u2.plt.addPlot(row=1, col=1, title= "Real-Time Plot: Phase", **{'font-size':'12pt'})
        self._plt1.showGrid(x=True, y=True) 
        self._plt1.setLabel('bottom', 'Samples', units='n')
        self._plt1.setLabel('left', 'Phase', units='deg')
        '''
        #--------------------------------------------------------------------------------------------------------------
        # Configures elements of the PyQtGraph plots: Multiple Plot amplitude and phase
        self._plt1 = pg.ViewBox()  
        self._plt0.showAxis('right')
        self._plt0.scene().addItem(self._plt1)
        self._plt0.getAxis('right').linkToView(self._plt1)
        self._plt1.setXLink(self._plt0)
        self._plt0.enableAutoRange(axis= 'y', enable = True)
        self._plt1.enableAutoRange(axis= 'y', enable = True)
        self._plt0.setLabel('right', 'Phase', units='deg', color=Constants.plot_colors[1], **{'font-size':'10pt'})
        
        #--------------------------------------------------------------------------------------------------------------
        # Configures elements of the PyQtGraph plots: resonance
        self._yaxis = AxisItem(orientation='left')  #NonScientificAxis(orientation='left')#
        #self._yaxis.setTickSpacing(levels=[(280, 0),(25, 0), (10, 0)]) #(20,1, None)
        self._xaxis = DateAxis(orientation='bottom')
        self._plt2 = self.PlotsWin.ui2.pltB.addPlot(row=0, col=2, title= title2, **{'font-size':'12pt'}, axisItems={"bottom":self._xaxis, "left":self._yaxis})
        self._plt2.showGrid(x=True, y=True) 
        self._plt2.setLabel('bottom', 'Time',units='s')
        self._plt2.setLabel('left', 'Resonance Frequency', units='Hz', color=Constants.plot_colors[2], **{'font-size':'10pt'})
        
        #--------------------------------------------------------------------------------------------------------------
        # Configures elements of the PyQtGraph plots: Multiple Plot resonance frequency and dissipation
        self._plt3 = pg.ViewBox()  
        self._plt2.showAxis('right')
        self._plt2.scene().addItem(self._plt3)
        self._plt2.getAxis('right').linkToView(self._plt3)
        self._plt3.setXLink(self._plt2)
        self._plt2.enableAutoRange(axis= 'y', enable = True) 
        self._plt3.enableAutoRange(axis= 'y', enable = True)
        self._plt2.setLabel('bottom', 'Time',units='s')
        self._plt2.setLabel('right', 'Dissipation', units='', color=Constants.plot_colors[3], **{'font-size':'10pt'})
        '''
        # indipentent plots
        # Configures elements of the PyQtGraph plots: Q-factor/dissipation
        self._plt3 = self.PlotsWin.ui2.plt.addPlot(row=1, col=2, title="Real-Time Plot: Dissipation", axisItems={'bottom': DateAxis(orientation='bottom')})
        self._plt3.showGrid(x=True, y=True) 
        self._plt3.setLabel('bottom', 'Time (hh:mm:ss)',units='')
        self._plt3.setLabel('left', 'Dissipation', units='')
        '''
        #---------------------------------------------------------------------------------------------------------------
        # Configures text comment and info
        self._plt2.clear()
        self._text1 = pg.TextItem('',(255, 255, 255), anchor=(0.5, 0.5))
        self._text1.setHtml("<span style='font-size: 14pt'>Welcome to the real-time openQCM Q-1 device GUI </span>")
        self._text2 = pg.TextItem('',(250, 250, 250), anchor=(0.5, 0.5))
        self._text2.setHtml("<span style='font-size: 10pt'>Please do not forget, if needed, to calibrate (in air) your device with the properly selected quartz before starting. </span>")
        self._text3 = pg.TextItem('',(250, 250, 250), anchor=(0.5, 0.5))
        self._text3.setHtml("<span style='font-size: 10pt'>Please check the command line on the side for useful information. </span>")
        self._text4 = pg.TextItem('',(250, 250, 250), anchor=(0.5, 0.5))
        self._text4.setHtml("<span style='font-size: 10pt'>(https://openqcm.com - info@openqcm.com) </span>")
        self._text1.setPos(0.5,0.6)
        self._text2.setPos(0.5,0.3)
        self._text3.setPos(0.5,0.2)
        self._text4.setPos(0.5,0.1)
        self._plt2.addItem(self._text1, ignoreBounds=True)
        self._plt2.addItem(self._text2, ignoreBounds=True)
        self._plt2.addItem(self._text3, ignoreBounds=True)
        self._plt2.addItem(self._text4, ignoreBounds=True)
        #######################
        #self._lr = pg.LinearRegionItem(orientation=pg.LinearRegionItem.Horizontal)
        #self._inf1 = pg.InfiniteLine(movable=True, angle=0, label='x={value:0.2f}', 
        #             labelOpts={'color': (100,200,100), 'fill': (200,200,200,50), 'movable': True})
        #######################
        #-----------------------------------------------------------------------------------------------------------------
        # Configures elements of the PyQtGraph plots: temperature
        self._plt4 = self.PlotsWin.ui2.plt.addPlot(row=0, col=3, title= title3, axisItems={'bottom': DateAxis(orientation='bottom')})
        self._plt4.showGrid(x=True, y=True) 
        self._plt4.setLabel('bottom', 'Time',units='s')
        self._plt4.setLabel('left', 'Temperature', units='°C', color=Constants.plot_colors[4], **{'font-size':'10pt'})
      
        
    ###########################################################################
    # Configures specific elements of the QTimers
    ########################################################################### 
    def _configure_timers(self):
        
        self._timer_plot = QtCore.QTimer(self)
        self._timer_plot.timeout.connect(self._update_plot)

    
    ###########################################################################
    # Configures the connections between signals and UI elements
    ########################################################################### 
    def _configure_signals(self):
        
        self.ControlsWin.ui1.pButton_Start.clicked.connect(self.start)
        self.ControlsWin.ui1.pButton_Stop.clicked.connect(self.stop)
        self.ControlsWin.ui1.pButton_Clear.clicked.connect(self.clear)
        self.ControlsWin.ui1.pButton_Reference.clicked.connect(self.reference)
        self.ControlsWin.ui1.sBox_Samples.valueChanged.connect(self._update_sample_size)
        self.ControlsWin.ui1.cBox_Source.currentIndexChanged.connect(self._source_changed)
        #--------
        self.InfoWin.ui3.pButton_Download.clicked.connect(self.start_download)
    
    ###########################################################################
    # Updates the sample size of the plot (now not used)
    ########################################################################### 
    def _update_sample_size(self):
        
        # This function is connected to the valueChanged signal of the sample Spin Box.
        if self.worker is not None:
            #Log.i(TAG, "Changing sample size")
            self.worker.reset_buffers(self.ControlsWin.ui1.sBox_Samples.value())

    
    ###########################################################################
    # Updates and redraws the graphics in the plot.
    ########################################################################### 
    def _update_plot(self):
        
        # This function is connected to the timeout signal of a QTimer
        self.worker.consume_queue1()
        self.worker.consume_queue2()
        self.worker.consume_queue3()
        self.worker.consume_queue4()
        self.worker.consume_queue5() 
        self.worker.consume_queue6() 
        
        # MEASUREMENT: dynamic frequency and dissipation labels at run-time
        ###################################################################
        if  self._get_source() == SourceType.serial:
            vector1 = self.worker.get_d1_buffer()
            vector2 = self.worker.get_d2_buffer()
            vectortemp = self.worker.get_d3_buffer()
            self._ser_error1,self._ser_error2, self._ser_control,self._ser_err_usb = self.worker.get_ser_error()
            #print(self._ser_err_usb, end='\r')
            if self._ser_err_usb <=1:
                if vector1.any:
                   # progressbar
                   if self._ser_control<=Constants.environment:
                      self._completed = self._ser_control*2

                   if str(vector1[0])=='nan' and not self._ser_error1 and not self._ser_error2:
                      label1 = 'processing...'
                      label2 = 'processing...'
                      label3 = 'processing...' 
                      labelstatus = 'Processing'
                      self.ControlsWin.ui1.infostatus.setStyleSheet('background: #ffff00; padding: 1px; border: 1px solid #cccccc') #ff8000
                      color_err = '#000000'   
                      labelbar = 'Please wait, processing early data...'

                   elif (str(vector1[0])=='nan' and (self._ser_error1 or self._ser_error2)):
                          if self._ser_error1 and self._ser_error2:
                            label1= ""
                            label2= ""
                            label3= ""
                            labelstatus = 'Warning'
                            color_err = '#ff0000'
                            labelbar = 'Warning: unable to apply half-power bandwidth method, lower and upper cut-off frequency not found'
                            self.ControlsWin.ui1.infostatus.setStyleSheet('background: #ff0000; padding: 1px; border: 1px solid #cccccc')
                          elif self._ser_error1:
                            label1= ""
                            label2= ""
                            label3= ""
                            labelstatus = 'Warning'
                            color_err = '#ff0000'
                            labelbar = 'Warning: unable to apply half-power bandwidth method, lower cut-off frequency (left side) not found'
                            self.ControlsWin.ui1.infostatus.setStyleSheet('background: #ff0000; padding: 1px; border: 1px solid #cccccc')
                          elif self._ser_error2:
                            label1= ""
                            label2= ""
                            label3= ""
                            labelstatus = 'Warning'
                            color_err = '#ff0000'
                            labelbar = 'Warning: unable to apply half-power bandwidth method, upper cut-off frequency (right side) not found'
                            self.ControlsWin.ui1.infostatus.setStyleSheet('background: #ff0000; padding: 1px; border: 1px solid #cccccc')
                   else:
                      if not self._ser_error1 and not self._ser_error2:
                          if not self._reference_flag:
                              d1=float("{0:.2f}".format(vector1[0]))
                              d2=float("{0:.4f}".format(vector2[0]*1e6))
                              d3=float("{0:.2f}".format(vectortemp[0]))
                          else:
                              a1= vector1[0]-self._reference_value_frequency
                              a2= vector2[0]-self._reference_value_dissipation
                              d1=float("{0:.2f}".format(a1))
                              d2=float("{0:.4f}".format(a2*1e6))
                              d3=float("{0:.2f}".format(vectortemp[0])) 
                          label1= str(d1)+ " Hz"
                          label2= str(d2)+ "e-06"
                          label3= str(d3)+ " °C" 
                          labelstatus = 'Monitoring'
                          color_err = '#000000'
                          labelbar = 'Monitoring!'
                          self.ControlsWin.ui1.infostatus.setStyleSheet('background: #00ff72; padding: 1px; border: 1px solid #cccccc')
                      else:
                          if self._ser_error1 and self._ser_error2:
                            label1= "-"
                            label2= "-"
                            label3= "-"
                            labelstatus = 'Warning'
                            color_err = '#ff0000'
                            labelbar = 'Warning: unable to apply half-power bandwidth method, lower and upper cut-off frequency not found'
                            self.ControlsWin.ui1.infostatus.setStyleSheet('background: #ff0000; padding: 1px; border: 1px solid #cccccc')
                          elif self._ser_error1:
                            label1= "-"
                            label2= "-"
                            label3= "-"
                            labelstatus = 'Warning'
                            color_err = '#ff0000'
                            labelbar = 'Warning: unable to apply half-power bandwidth method, lower cut-off frequency (left side) not found'
                            self.ControlsWin.ui1.infostatus.setStyleSheet('background: #ff0000; padding: 1px; border: 1px solid #cccccc')
                          elif self._ser_error2:
                            label1= "-"
                            label2= "-"
                            label3= "-"
                            labelstatus = 'Warning'
                            color_err = '#ff0000'
                            labelbar = 'Warning: unable to apply half-power bandwidth method, upper cut-off frequency (right side) not found'
                            self.ControlsWin.ui1.infostatus.setStyleSheet('background: #ff0000; padding: 1px; border: 1px solid #cccccc')
                             
                   self.InfoWin.ui3.l6a.setText("<font color=#0000ff > Temperature </font>" + label3) 
                   self.InfoWin.ui3.l6.setText("<font color=#0000ff > Dissipation </font>" + label2)
                   self.InfoWin.ui3.l7.setText("<font color=#0000ff > Resonance Frequency </font>" + label1)
                   self.ControlsWin.ui1.infostatus.setText("<font color=#000000 > Program Status </font>" + labelstatus)
                   self.ControlsWin.ui1.infobar.setText("<font color=#0000ff> Infobar </font><font color={}>{}</font>".format(color_err,labelbar))
                   # progressbar 
                   self.ControlsWin.ui1.progressBar.setValue(self._completed+2)
            
            elif self._ser_err_usb >1:
                 PopUp.warning(self, Constants.app_title, "Warning: USB cable device disconnected!")  
                 self.stop() 
        
        # CALIBRATION: dynamic info in infobar at run-time
        ##################################################
        elif self._get_source() == SourceType.calibration:
            # flag for terminating calibration
            stop_flag=0
            self.ControlsWin.ui1.pButton_Stop.setEnabled(False)
            vector1 = self.worker.get_value1_buffer()
            # vector2[0] and vector3[0] flag error
            vector2 = self.worker.get_t3_buffer()
            vector3 = self.worker.get_d3_buffer()
            #print(vector1[0],vector2[0],vector3[0])
            label1 = 'not available'
            label2 = 'not available'
            label3 = 'not available' 
            labelstatus = 'Calibration Processing'
            color_err = '#000000'
            labelbar = 'The operation might take just over a minute to complete... please wait...'
            self.ControlsWin.ui1.infostatus.setStyleSheet('background: #ffff00; padding: 1px; border: 1px solid #cccccc')
            
            # progressbar
            error1,error2,error3,len_buffer = self.worker.get_ser_error()
            #print(len_buffer)
            if len_buffer <830000:
               self._completed = (len_buffer/830000)*100
               #print(self._completed)

            # calibration buffer empty
            if vector1[0]== 0 and vector3[0]==1:
              label1 = 'not available'
              label2 = 'not available'
              label3 = 'not available' 
              color_err = '#ff0000'
              labelstatus = 'Calibration Warning'
              self.ControlsWin.ui1.infostatus.setStyleSheet('background: #ff0000; padding: 1px; border: 1px solid #cccccc')
              labelbar = 'Calibration Warning: empty buffer! Please, repeat the Calibration after disconnecting/reconnecting Device!'
              stop_flag=1
            # calibration buffer empty and ValueError from the serial port 
            elif vector1[0]== 0 and vector2[0]==1:
              label1 = 'not available'
              label2 = 'not available'
              label3 = 'not available' 
              color_err = '#ff0000'
              labelstatus = 'Calibration Warning'
              self.ControlsWin.ui1.infostatus.setStyleSheet('background: #ff0000; padding: 1px; border: 1px solid #cccccc')
              labelbar = 'Calibration Warning: empty buffer/ValueError! Please, repeat the Calibration after disconnecting/reconnecting Device!'
              stop_flag=1 
            # calibration buffer not empty
            elif vector1[0]!= 0:
              label1 = 'not available'
              label2 = 'not available'
              label3 = 'not available' 
              labelstatus = 'Calibration Processing'
              color_err = '#000000'
              labelbar = 'The operation might take just over a minute to complete... please wait...'
              if vector2[0]== 0 and vector3[0]== 0:
                 labelstatus = 'Calibration Success'
                 self.ControlsWin.ui1.infostatus.setStyleSheet('background: #00ff72; padding: 1px; border: 1px solid #cccccc')
                 color_err = '#000000' 
                 labelbar = 'Calibration Success for baseline correction!' 
                 stop_flag=1
                 #print(self._k) #progressbar value 143
              elif vector2[0]== 1 or vector3[0]== 1:
                 color_err = '#ff0000'
                 labelstatus = 'Calibration Warning'
                 self.ControlsWin.ui1.infostatus.setStyleSheet('background: #ff0000; padding: 1px; border: 1px solid #cccccc')
                 if vector2[0]== 1:
                   labelbar = 'Calibration Warning: ValueError! Please, repeat the Calibration after disconnecting/reconnecting Device!'
                   stop_flag=1 ##
                 elif vector3[0]== 1: 
                   labelbar = 'Calibration Warning: Error during Peak Detection, incompatible peaks number or frequencies! Please, repeat the Calibration!'
                   stop_flag=1 ##

            self.InfoWin.ui3.l6a.setText("<font color=#0000ff>  Dissipation </font>" + label3) 
            self.InfoWin.ui3.l6.setText("<font color=#0000ff>  Dissipation </font>" + label2)
            self.InfoWin.ui3.l7.setText("<font color=#0000ff>  Resonance Frequency </font>" + label1)
            self.ControlsWin.ui1.infostatus.setText("<font color=#000000> Program Status </font>" + labelstatus)  
            self.ControlsWin.ui1.infobar.setText("<font color=#0000ff> Infobar </font><font color={}>{}</font>".format(color_err,labelbar)) 
            # progressbar -------------
            self.ControlsWin.ui1.progressBar.setValue(self._completed+1)
            
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
        ############################################################################################################################
        # REFERENCE SET 
        ############################################################################################################################
        if self._reference_flag:
            self._plt2.setLabel('left', 'Resonance Frequency', units='Hz', color=Constants.plot_colors[6], **{'font-size':'10pt'})
            self._plt2.setLabel('right', 'Dissipation', units='', color=Constants.plot_colors[7], **{'font-size':'10pt'})
            self.InfoWin.ui3.inforef1.setText("<font color=#0000ff > Ref. Frequency </font>" + self._labelref1)               
            self.InfoWin.ui3.inforef2.setText("<font color=#0000ff > Ref. Dissipation </font>" + self._labelref2)
            
            ###################################################################
            # Amplitude and phase multiple Plot
            def updateViews1():
                self._plt0.clear()
                self._plt1.clear()
                self._plt1.setGeometry(self._plt0.vb.sceneBoundingRect())
                self._plt1.linkedViewChanged(self._plt0.vb, self._plt1.XAxis)
            
            # updates for multiple plot y-axes
            updateViews1()
            self._plt0.vb.sigResized.connect(updateViews1) 
            self._plt0.plot(x=self._readFREQ,y=self.worker.get_value1_buffer(),pen=Constants.plot_colors[0])
            self._plt1.addItem(pg.PlotCurveItem(x=self._readFREQ,y=self.worker.get_value2_buffer(),pen=Constants.plot_colors[1]))
            
            ###################################################################
            # Resonance frequency and dissipation multiple Plot
            def updateViews2():
                self._plt2.clear()
                self._plt3.clear()
                self._plt3.setGeometry(self._plt2.vb.sceneBoundingRect())
                self._plt3.linkedViewChanged(self._plt2.vb, self._plt3.XAxis)
            
            # updates for multiple plot y-axes
            updateViews2()
            self._plt2.vb.sigResized.connect(updateViews2) 
            self._vector_1 = np.array(self.worker.get_d1_buffer())-self._reference_value_frequency          
            self._plt2.plot(x=self.worker.get_t1_buffer(),y=self._vector_1,pen=Constants.plot_colors[6]) #2
            self._vector_2 = np.array(self.worker.get_d2_buffer())-self._reference_value_dissipation 
            #Prevent the user from zooming/panning out of this specified region
            if self._get_source() == SourceType.serial:
               #dy = [(value, str(value)) for value in (range(int(min(self._readFREQ)), int(max(self._readFREQ)+1)))]
               #self._yaxis.setTicks([dy])
               #tickBottom = {self._readFREQ[250]:self._readFREQ[0],self._readFREQ[-1]:self._readFREQ[250]}
               #self._yaxis.setTicks([tickBottom.items()])
               self._plt2.setLimits(yMax=self._vector_reference_frequency[-1],yMin=self._vector_reference_frequency[0], minYRange=5)
               self._plt3.setLimits(yMax=self._vector_reference_dissipation[-1],yMin=self._vector_reference_dissipation[0], minYRange=1e-7) 
               self._plt4.setLimits(yMax=50,yMin=-10)
            self._plt3.addItem(pg.PlotCurveItem(self.worker.get_t2_buffer(),self._vector_2,pen=Constants.plot_colors[7]))
            
            ###################################################################
            # Temperature plot
            self._plt4.clear() 
            self._plt4.plot(x= self.worker.get_t3_buffer(),y=self.worker.get_d3_buffer(),pen=Constants.plot_colors[4])
            
        ###########################################################################################################################
        # REFERENCE NOT SET 
        ###########################################################################################################################  
        else:   
            self._plt2.setLabel('left', 'Resonance Frequency', units='Hz', color=Constants.plot_colors[2], **{'font-size':'10pt'})
            self._plt2.setLabel('right', 'Dissipation', units='', color=Constants.plot_colors[3], **{'font-size':'10pt'})
            self.InfoWin.ui3.inforef1.setText("<font color=#0000ff > Ref. Frequency </font>" + self._labelref1)               
            self.InfoWin.ui3.inforef2.setText("<font color=#0000ff > Ref. Dissipation </font>" + self._labelref2)
            
            ###################################################################
            # Amplitude and phase multiple Plot
            def updateViews1():
                self._plt0.clear()
                self._plt1.clear()
                self._plt1.setGeometry(self._plt0.vb.sceneBoundingRect())
                self._plt1.linkedViewChanged(self._plt0.vb, self._plt1.XAxis)
            # updates for multiple plot y-axes
            updateViews1()
            self._plt0.vb.sigResized.connect(updateViews1) 
            self._plt0.plot(x=self._readFREQ,y=self.worker.get_value1_buffer(),pen=Constants.plot_colors[0])
            self._plt1.addItem(pg.PlotCurveItem(x=self._readFREQ,y=self.worker.get_value2_buffer(),pen=Constants.plot_colors[1]))
            ###
            #img = pg.QtGui.QGraphicsPixmapItem(pg.QtGui.QPixmap('favicon.png'))
            #img.scale(1, -1)
            #self._plt1.addItem(img)
            
            ###################################################################
            # Resonance frequency and dissipation multiple Plot
            def updateViews2():
                self._plt2.clear()
                self._plt3.clear()
                self._plt3.setGeometry(self._plt2.vb.sceneBoundingRect())
                self._plt3.linkedViewChanged(self._plt2.vb, self._plt3.XAxis)
            # updates for multiple plot y-axes
            updateViews2()
            #text = pg.TextItem(html='<div style="text-align: center"><span style="color: #ff0000;"> %s</span></div>' % "prova",anchor=(0.0, 0.0)) 
            #text.setPos(0.0, 0.0) #text########
            #self._plt3.addItem(text)
            self._plt2.vb.sigResized.connect(updateViews2) 
            self._plt2.plot(x= self.worker.get_t1_buffer(),y=self.worker.get_d1_buffer(),pen=Constants.plot_colors[2])
            #Prevent the user from zooming/panning out of this specified region
            if self._get_source() == SourceType.serial:
               #dy = [(value, str(value)) for value in (range(int(min(self._readFREQ)), int(max(self._readFREQ)+1)))]
               #self._yaxis.setTicks([dy])
               #tickBottom = {self._readFREQ[250]:self._readFREQ[0],self._readFREQ[-1]:self._readFREQ[250]}
               #self._yaxis.setTicks([tickBottom.items()])
               self._plt2.setLimits(yMax=self._readFREQ[-1],yMin=self._readFREQ[0], minYRange=10, maxYRange=10000) 
               self._plt3.setLimits(yMax=(self._readFREQ[-1]-self._readFREQ[0])/self._readFREQ[0],yMin=0, minYRange=1e-6)
               self._plt4.setLimits(yMax=50,yMin=-10, minYRange=0.5)
            self._plt3.addItem(pg.PlotCurveItem(self.worker.get_t2_buffer(),self.worker.get_d2_buffer(),pen=Constants.plot_colors[3]))
            
            ##############################
            # Add  lines with labels
            #inf1 = pg.InfiniteLine(movable=True, angle=90, label='x={value:0.2f}', 
            #self._plt2.addItem(self._inf1)
            #self._plt2.addItem(self._lr)
            ##############################
            '''
            # Resonance frequency plot
            self._plt2.clear()
            self._plt2.plot(x= self.worker.get_t1_buffer(),y=self.worker.get_d1_buffer(),pen=Constants.plot_colors[2])
            # dissipation plot
            self._plt3.clear()
            self._plt3.plot(x= self.worker.get_t2_buffer(),y=self.worker.get_d2_buffer(),pen=Constants.plot_colors[3])
            '''
            ###################################################################
            # Temperature plot
            self._plt4.clear() 
            self._plt4.plot(x= self.worker.get_t3_buffer(),y=self.worker.get_d3_buffer(),pen=Constants.plot_colors[4])     
          
    ###########################################################################################################################################    
        
        
    ###########################################################################
    # Updates the source and depending boxes on change
    ###########################################################################  
    def _source_changed(self):
        
        # It is connected to the indexValueChanged signal of the Source ComboBox.
        if self._get_source() == SourceType.serial:
           print(TAG, "Scanning the source: {}".format(Constants.app_sources[0]))  # self._get_source().name
           Log.i(TAG, "Scanning the source: {}".format(Constants.app_sources[0]))
        elif self._get_source() == SourceType.calibration:
           print(TAG, "Scanning the source: {}".format(Constants.app_sources[1]))  # self._get_source().name
           Log.i(TAG, "Scanning the source: {}".format(Constants.app_sources[1]))
        
        # Clears boxes before adding new
        self.ControlsWin.ui1.cBox_Port.clear()
        self.ControlsWin.ui1.cBox_Speed.clear()
        
        # Gets the current source type
        source = self._get_source()
        ports = self.worker.get_source_ports(source)
        speeds = self.worker.get_source_speeds(source)

        if ports is not None:
            self.ControlsWin.ui1.cBox_Port.addItems(ports)
        if speeds is not None:
            self.ControlsWin.ui1.cBox_Speed.addItems(speeds)
        if self._get_source() == SourceType.serial:
            self.ControlsWin.ui1.cBox_Speed.setCurrentIndex(len(speeds) - 1)

    
    ###########################################################################
    # Gets the current source type
    ########################################################################### 
    def _get_source(self):
        
        #:rtype: SourceType.
        return SourceType(self.ControlsWin.ui1.cBox_Source.currentIndex())
    
    
    ###########################################################################
    # Cleans history plot 
    ###########################################################################    
    def clear(self):
        support=self.worker.get_d1_buffer()
        if support.any:
            if str(support[0])!='nan':
                print(TAG, "All Plots Cleared!", end='\r')
                self._update_sample_size()
                self._plt2.clear()
                self._plt3.clear()
                self._plt4.clear()
        
        
    ###########################################################################
    # Reference set/reset
    ###########################################################################     
    def reference(self):
        import numpy as np
        #import sys
        support=self.worker.get_d1_buffer()
        if support.any:
            if str(support[0])!='nan':
                ref_vector1 = [c for c in self.worker.get_d1_buffer() if ~np.isnan(c)]
                ref_vector2 = [c for c in self.worker.get_d2_buffer() if ~np.isnan(c)]
                self._reference_value_frequency = ref_vector1[0]
                self._reference_value_dissipation = ref_vector2[0]
                #sys.stdout.write("\033[K") #clear line
                if self._reference_flag:
                    self._reference_flag = False
                    print(TAG, "Reference reset!   ", end='\r')
                    self._labelref1 = "not set"
                    self._labelref2 = "not set"
                else:
                    self._reference_flag = True
                    d1=float("{0:.2f}".format(self._reference_value_frequency))
                    d2=float("{0:.4f}".format(self._reference_value_dissipation*1e6))
                    self._labelref1 = str(d1)+ "Hz"
                    self._labelref2 = str(d2)+ "e-06"
                    print(TAG, "Reference set!     ", end='\r')
                    self._vector_reference_frequency[:] = [s - self._reference_value_frequency for s in self._readFREQ]
                    xs = np.array(np.linspace(0, ((self._readFREQ[-1]-self._readFREQ[0])/self._readFREQ[0]), len(self._readFREQ)))
                    self._vector_reference_dissipation = xs-self._reference_value_dissipation
                    
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
                   '3': '3',}
        r_types = {
                   '1': 'A',
                   '2': 'B',
                   '3': 'C',}
        if self._internet_connected:
           color = '#00c600'
           labelweb2 = 'ONLINE'
           print (TAG,'Checking your internet connection {} '.format(labelweb2))
           tables = pd.read_html('https://openqcm.com/shared/news.html', index_col=0, header=0, match='1')
           df = tables[0]
           # create empty list of string  
           self._webinfo = ["" for x in range(len(df.columns)*len(df.index))] #len(df.columns)*len(df.index)=9
           # row acess mode to Pandas dataframe
           k=0
           for j in [1,2,3]:
              for i in [1,2,3]:
                  self._webinfo[k]= str(df.loc[r_types[str(j)], c_types[str(i)]])
                  k+=1
            # check for update
           if self._webinfo[0] == Constants.app_version:
              labelweb3 = 'last version installed!' 
           else:
              labelweb3 = 'version {} available!'.format(self._webinfo[0]) 
              self.InfoWin.ui3.pButton_Download.setEnabled(True)
        else:
           color = '#ff0000'
           labelweb2 = 'OFFLINE'
           labelweb3 = 'Offline, unable to check'
           print (TAG,'Checking your internet connection {} '.format(labelweb2)) 
           
        self.InfoWin.ui3.lweb2.setText("<font color=#0000ff > Checking your internet connection &nbsp;&nbsp;&nbsp;&nbsp;</font><font size=4 color={}>{}</font>".format(color,labelweb2)) 
        self.InfoWin.ui3.lweb3.setText("<font color=#0000ff > Software update status </font>" + labelweb3)  

    ########################################################################### 
    # Opens webpage for download
    ###########################################################################
    def start_download(self):
        import webbrowser
        webbrowser.open('https://openqcm.com/shared/q-1/openQCM_Q-1_py_v2.0.zip')
        
       