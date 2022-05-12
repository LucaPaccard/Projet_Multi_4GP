
from PyQt5 import QtCore, QtGui, QtWidgets

class Ui_Controls(object):
    
    def setupUi(self, MainWindow):
        MainWindow.setObjectName("MainWindow")
        MainWindow.resize(1090, 70)
        MainWindow.setMinimumSize(QtCore.QSize(0, 0))
        MainWindow.setStyleSheet("")
        MainWindow.setTabShape(QtWidgets.QTabWidget.Rounded)
        self.centralwidget = QtWidgets.QWidget(MainWindow)
        self.centralwidget.setObjectName("centralwidget")
        self.gridLayout = QtWidgets.QGridLayout(self.centralwidget)
        self.gridLayout.setObjectName("gridLayout")
        self.Layout_controls = QtWidgets.QGridLayout()
        self.Layout_controls.setObjectName("Layout_controls")
               
        # frequency/quartz combobox -------------------------------------------
        self.cBox_Speed = QtWidgets.QComboBox(self.centralwidget)
        self.cBox_Speed.setEditable(False)
        self.cBox_Speed.setObjectName("cBox_Speed")
        self.Layout_controls.addWidget(self.cBox_Speed, 4, 1, 1, 1)
        
        # stop button ---------------------------------------------------------
        self.pButton_Stop = QtWidgets.QPushButton(self.centralwidget)
        self.pButton_Stop.setObjectName("pButton_Stop")
        self.Layout_controls.addWidget(self.pButton_Stop, 3, 4, 1, 1)
        
        # COM port combobox ---------------------------------------------------
        self.cBox_Port = QtWidgets.QComboBox(self.centralwidget)
        self.cBox_Port.setEditable(False)
        self.cBox_Port.setObjectName("cBox_Port")
        self.Layout_controls.addWidget(self.cBox_Port, 2, 1, 1, 1)
        
        # Operation mode - source ---------------------------------------------
        self.cBox_Source = QtWidgets.QComboBox(self.centralwidget)
        self.cBox_Source.setObjectName("cBox_Source")
        self.Layout_controls.addWidget(self.cBox_Source, 2, 0, 1, 1)
        
        # start button --------------------------------------------------------
        self.pButton_Start = QtWidgets.QPushButton(self.centralwidget)
        self.pButton_Start.setIcon(QtGui.QIcon(QtGui.QPixmap("favicon.ico"))) #.png
        self.pButton_Start.setMinimumSize(QtCore.QSize(0, 0))
        self.pButton_Start.setObjectName("pButton_Start")
        self.Layout_controls.addWidget(self.pButton_Start, 2, 4, 1, 1)
        
        # clear plots button --------------------------------------------------
        self.pButton_Clear = QtWidgets.QPushButton(self.centralwidget)
        self.pButton_Clear.setMinimumSize(QtCore.QSize(0, 0))
        self.pButton_Clear.setObjectName("pButton_Clear")
        self.Layout_controls.addWidget(self.pButton_Clear, 2, 3, 1, 1)
        
        # reference button ----------------------------------------------------
        self.pButton_Reference = QtWidgets.QPushButton(self.centralwidget)
        self.pButton_Reference.setMinimumSize(QtCore.QSize(0, 0))
        self.pButton_Reference.setObjectName("pButton_Reference")
        self.Layout_controls.addWidget(self.pButton_Reference, 3, 3, 1, 1)
        
        # samples SpinBox -----------------------------------------------------
        self.sBox_Samples = QtWidgets.QSpinBox(self.centralwidget)
        self.sBox_Samples.setMinimum(1)
        self.sBox_Samples.setMaximum(100000)
        self.sBox_Samples.setProperty("value", 500)
        self.sBox_Samples.setObjectName("sBox_Samples")
        #self.sBox_Samples.setEnabled(False)
        self.Layout_controls.addWidget(self.sBox_Samples, 2, 2, 1, 1)
        
        # export file CheckBox ------------------------------------------------
        self.chBox_export = QtWidgets.QCheckBox(self.centralwidget)
        self.chBox_export.setEnabled(True)
        self.chBox_export.setObjectName("chBox_export")
        self.Layout_controls.addWidget(self.chBox_export, 4, 2, 1, 1)
        
        # Samples Number / History Buffer Size---------------------------------
        self.l5 = QtWidgets.QLabel()
        self.l5.setStyleSheet('background: #008EC0; padding: 1px;')
        self.l5.setText("<font color=#ffffff > Samples Number / History Buffer Size </font>")
        self.Layout_controls.addWidget(self.l5, 1, 2, 1, 1) 
        
        # Control Buttons------------------------------------------------------
        self.l = QtWidgets.QLabel()
        self.l.setStyleSheet('background: #008EC0; padding: 1px;')
        self.l.setText("<font color=#ffffff > Control Buttons </font>")
        self.Layout_controls.addWidget(self.l, 1, 3, 1, 2)     
        
        # Operation Mode ------------------------------------------------------
        self.l0 = QtWidgets.QLabel()
        self.l0.setStyleSheet('background: #008EC0; padding: 1px;')
        self.l0.setText("<font color=#ffffff >Operation Mode</font> </a>")
        self.Layout_controls.addWidget(self.l0, 1, 0, 1, 1)
        
        # Resonance Frequency / Quartz Sensor ---------------------------------
        self.l2 = QtWidgets.QLabel()
        self.l2.setStyleSheet('background: #008EC0; padding: 1px;')
        self.l2.setText("<font color=#ffffff > Resonance Frequency / Quartz Sensor </font>")
        self.l2.setFixedHeight(15)
        self.Layout_controls.addWidget(self.l2, 3, 1, 1, 1)
        
        # Serial COM Port -----------------------------------------------------
        self.l1 = QtWidgets.QLabel()
        self.l1.setStyleSheet('background: #008EC0; padding: 1px;')
        self.l1.setText("<font color=#ffffff > Serial COM Port </font>")
        self.Layout_controls.addWidget(self.l1, 1, 1, 1, 1)
        
        #######################################################################
        # openQCM logo---------------------------------------------------------
        self.l3 = QtWidgets.QLabel()
        self.l3.setAlignment(QtCore.Qt.AlignRight)
        self.Layout_controls.addWidget(self.l3, 4, 5, 1, 1)
        self.l3.setPixmap(QtGui.QPixmap("openqcm-logo.png"))
        
        # openQCM link --------------------------------------------------------
        self.l4 = QtWidgets.QLabel()
        self.Layout_controls.addWidget(self.l4, 3, 5, 1, 1)
        def link(linkStr):
            QtGui.QDesktopServices.openUrl(QtCore.QUrl(linkStr))
        self.l4.linkActivated.connect(link)
        self.l4.setAlignment(QtCore.Qt.AlignRight)
        self.l4.setText('<a href="https://openqcm.com/openqcm-q-1-software"> <font size=4 color=#008EC0 >openqcm.com</font>') #&nbsp;
                        
        # info@openqcm.com Mail -----------------------------------------------
        self.lmail = QtWidgets.QLabel()
        self.Layout_controls.addWidget(self.lmail, 2, 5, 1, 1)
        def linkmail(linkStr):
            QtGui.QDesktopServices.openUrl(QtCore.QUrl(linkStr))
        self.lmail.linkActivated.connect(linkmail)
        self.lmail.setAlignment(QtCore.Qt.AlignRight)
        self.lmail.setText('<a href="mailto:info@openqcm.com"> <font color=#008EC0 >info@openqcm.com</font>')

        # Save file -----------------------------------------------------------
        self.infosave = QtWidgets.QLabel()
        self.infosave.setStyleSheet('background: #008EC0; padding: 1px;')
        #self.infosave.setAlignment(QtCore.Qt.AlignCenter) 
        self.infosave.setFixedHeight(15)
        self.infosave.setText("<font color=#ffffff > Save file </font>")
        self.Layout_controls.addWidget(self.infosave, 3, 2, 1, 1)
        
        # Program Status standby ----------------------------------------------  
        self.infostatus = QtWidgets.QLabel()
        self.infostatus.setStyleSheet('background: white; padding: 1px; border: 1px solid #cccccc')
        self.infostatus.setAlignment(QtCore.Qt.AlignCenter) 
        self.infostatus.setText("<font color=#000000 > Program Status standby </font>")
        self.Layout_controls.addWidget(self.infostatus, 4, 3, 1, 2)
        
        # Infobar -------------------------------------------------------------
        self.infobar = QtWidgets.QLabel()
        self.infobar.setStyleSheet('background: white; padding: 1px; border: 1px solid #cccccc')
        #self.infobar.setAlignment(QtCore.Qt.AlignCenter) 
        self.infobar.setText("<font color=#0000ff > Infobar </font>")
        self.Layout_controls.addWidget(self.infobar, 0, 0, 1, 5)
        ##
        #self.progressBar = QtWidgets.QProgressBar()
        #self.progressBar.setProperty("value", 0)
        #self.progressBar.setObjectName("progressBar")
        #self.Layout_controls.addWidget(self.progressBar, 5, 2, 1, 2)
        #self.l1a = QtWidgets.QLabel()
        ##
        # ---------------------------------------------------------------------
        
        self.gridLayout.addLayout(self.Layout_controls, 7, 1, 1, 1)
        MainWindow.setCentralWidget(self.centralwidget)

        self.retranslateUi(MainWindow)
        QtCore.QMetaObject.connectSlotsByName(MainWindow)

    def retranslateUi(self, MainWindow):
        _translate = QtCore.QCoreApplication.translate
        MainWindow.setWindowTitle(_translate("MainWindow", "OPENQCM Q-1 - Setup/Control GUI"))
        MainWindow.setWindowIcon(QtGui.QIcon('favicon.ico')) #.png
        self.pButton_Stop.setText(_translate("MainWindow", "STOP"))
        self.pButton_Start.setText(_translate("MainWindow", "START"))
        self.pButton_Clear.setText(_translate("MainWindow", "Clear Plots"))
        self.pButton_Reference.setText(_translate("MainWindow", "Set/Reset Reference"))
        self.sBox_Samples.setSuffix(_translate("MainWindow", " / 180 min"))
        self.sBox_Samples.setPrefix(_translate("MainWindow", ""))
        self.chBox_export.setText(_translate("MainWindow", "Txt Export Sweep File"))

'''
if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    MainWindow = QtWidgets.QMainWindow()
    ui = Ui_Controls()
    ui.setupUi(MainWindow)
    MainWindow.show()
    sys.exit(app.exec_())
'''
