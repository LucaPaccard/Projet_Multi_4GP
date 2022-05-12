
from PyQt5 import QtCore, QtGui, QtWidgets
from pyqtgraph import GraphicsLayoutWidget

class Ui_Controls(object): #QtWidgets.QMainWindow
    
    def setupUi(self, MainWindow1):
        MainWindow1.setObjectName("MainWindow1")
        #MainWindow1.setGeometry(50, 50, 975, 70)
        # MainWindow1.setFixedSize(980, 150)
        #MainWindow1.resize(550, 50)
        MainWindow1.setMinimumSize(QtCore.QSize(980, 150))
        MainWindow1.resize(980, 150)
        MainWindow1.setStyleSheet("")
        MainWindow1.setTabShape(QtWidgets.QTabWidget.Rounded) 
        self.centralwidget = QtWidgets.QWidget(MainWindow1)
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
        self.pButton_Stop.setIcon(QtGui.QIcon(QtGui.QPixmap("stop_icon.ico"))) #.png
        self.pButton_Stop.setMinimumSize(QtCore.QSize(0, 0))
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
        self.pButton_Start.setIcon(QtGui.QIcon(QtGui.QPixmap("start_icon.ico")))
        self.pButton_Start.setMinimumSize(QtCore.QSize(0, 0))
        self.pButton_Start.setObjectName("pButton_Start")
        self.Layout_controls.addWidget(self.pButton_Start, 2, 4, 1, 1)
        
        # clear plots button --------------------------------------------------
        self.pButton_Clear = QtWidgets.QPushButton(self.centralwidget)
        self.pButton_Clear.setIcon(QtGui.QIcon(QtGui.QPixmap("clear_icon.ico")))
        self.pButton_Clear.setMinimumSize(QtCore.QSize(0, 0))
        self.pButton_Clear.setObjectName("pButton_Clear")
        self.Layout_controls.addWidget(self.pButton_Clear, 2, 3, 1, 1)
        
        # reference button ----------------------------------------------------
        self.pButton_Reference = QtWidgets.QPushButton(self.centralwidget)
        #self.pButton_Reference.setIcon(QtGui.QIcon(QtGui.QPixmap("ref_icon.ico")))
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
        self.l5.setFixedHeight(15)
        self.Layout_controls.addWidget(self.l5, 1, 2, 1, 1) 
        
        # Control Buttons------------------------------------------------------
        self.l = QtWidgets.QLabel()
        self.l.setStyleSheet('background: #008EC0; padding: 1px;')
        self.l.setText("<font color=#ffffff > Control Buttons </font>")
        self.l.setFixedHeight(15)
        self.Layout_controls.addWidget(self.l, 1, 3, 1, 2)     
        
        # Operation Mode ------------------------------------------------------
        self.l0 = QtWidgets.QLabel()
        self.l0.setStyleSheet('background: #008EC0; padding: 1px;')
        self.l0.setText("<font color=#ffffff >Operation Mode</font> </a>")
        self.l0.setFixedHeight(15)
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
        self.l1.setFixedHeight(15)
        self.Layout_controls.addWidget(self.l1, 1, 1, 1, 1)
        
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
        self.Layout_controls.addWidget(self.lmail, 2, 5, 1, 1) #25 40
        def linkmail(linkStr):
            QtGui.QDesktopServices.openUrl(QtCore.QUrl(linkStr))
        self.lmail.linkActivated.connect(linkmail)
        self.lmail.setAlignment(QtCore.Qt.AlignRight)
        #self.lmail.setAlignment(QtCore.Qt.AlignLeft)
        self.lmail.setText('<a href="mailto:info@openqcm.com"> <font color=#008EC0 >info@openqcm.com</font>')
        
        # software user guide --------------------------------------------------------
        
        self.lg = QtWidgets.QLabel()
        self.Layout_controls.addWidget(self.lg, 1, 5, 1, 1) #30
        def link(linkStr):
            QtGui.QDesktopServices.openUrl(QtCore.QUrl(linkStr))
        self.lg.linkActivated.connect(link)
        self.lg.setAlignment(QtCore.Qt.AlignRight)
        self.lg.setText('<a href="https://openqcm.com/shared/q-1/openQCM_Q-1-userguide-v2.0.pdf"> <font color=#008EC0 >User Guide</font>') #&nbsp;
        #####################################
        '''
        self.ico = QtWidgets.QLabel()
        self.Layout_controls.addWidget(self.ico, 2, 5, 1, 1)
        self.title = QtWidgets.QLabel()
        def link(linkStr):
            QtGui.QDesktopServices.openUrl(QtCore.QUrl(linkStr))
        self.title.linkActivated.connect(link)
        self.title.setText('<a href="https://openqcm.com/openqcm-q-1-software"> <font color=#008EC0 >user guide</font>')
        self.Layout_controls.addWidget(self.title, 2, 5, 1, 1)
        self.pixmap = QtGui.QPixmap("guide.ico")
        self.ico.setPixmap(self.pixmap)
        self.ico.setAlignment(QtCore.Qt.AlignRight)
        self.title.setMinimumHeight(self.pixmap.height())
        self.title.setAlignment(QtCore.Qt.AlignRight)
        '''
        #####################################
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
        
        # Progressbar -------------------------------------------------------------
        styleBar="""
                    QProgressBar
                    {
                     border: 0.5px solid grey;
                     border-radius: 1px;
                     text-align: center;
                    }
                     QProgressBar::chunk
                    {
                     background-color: qlineargradient(spread:pad, x1:0, y1:0, x2:1, y2:1, stop:0 rgba(182, 182, 182, 200), stop:1 rgba(209, 209, 209, 200));
                    }
                 """#background:url("openqcm-logo.png")
        self.progressBar = QtWidgets.QProgressBar()
        self.progressBar.setStyleSheet(styleBar) 
        self.progressBar.setProperty("value", 0)
        self.progressBar.setGeometry(QtCore.QRect(0, 0, 50, 10))
        self.progressBar.setObjectName("progressBar")
        self.Layout_controls.setColumnStretch(0,0)
        self.Layout_controls.setColumnStretch(1,2)
        self.Layout_controls.setColumnStretch(2,2)
        self.Layout_controls.setColumnStretch(3,2)
        self.Layout_controls.setColumnStretch(4,2)
        self.Layout_controls.addWidget(self.progressBar, 0, 5, 1, 1)
        # ---------------------------------------------------------------------
        
        self.gridLayout.addLayout(self.Layout_controls, 7, 1, 1, 1)
        MainWindow1.setCentralWidget(self.centralwidget) 

        self.retranslateUi(MainWindow1)
        QtCore.QMetaObject.connectSlotsByName(MainWindow1)

    def retranslateUi(self, MainWindow):
        _translate = QtCore.QCoreApplication.translate
        MainWindow.setWindowTitle(_translate("MainWindow", " Real-Time openQCM GUI - 2.1.1 BETA TEST - Setup/Control GUI"))
        MainWindow.setWindowIcon(QtGui.QIcon('favicon.ico')) #.png
        self.pButton_Stop.setText(_translate("MainWindow", " STOP"))
        self.pButton_Start.setText(_translate("MainWindow", "START"))
        self.pButton_Clear.setText(_translate("MainWindow", "Clear Plots"))
        self.pButton_Reference.setText(_translate("MainWindow", "Set/Reset Reference"))
        self.sBox_Samples.setSuffix(_translate("MainWindow", " / 180 min"))
        self.sBox_Samples.setPrefix(_translate("MainWindow", ""))
        self.chBox_export.setText(_translate("MainWindow", "Txt Export Sweep File"))


#######################################################################################################################

class Ui_Plots(object): 
    def setupUi(self, MainWindow2):
        MainWindow2.setObjectName("MainWindow2")
        #MainWindow2.setGeometry(100, 100, 890, 750)
        #MainWindow2.setFixedSize(1091, 770)
        #MainWindow2.resize(1091, 770)
        MainWindow2.setMinimumSize(QtCore.QSize(705, 518))
        MainWindow2.setStyleSheet("")
        MainWindow2.setTabShape(QtWidgets.QTabWidget.Rounded) 
        self.centralwidget = QtWidgets.QWidget(MainWindow2)
        self.centralwidget.setObjectName("centralwidget")
        self.gridLayout = QtWidgets.QGridLayout(self.centralwidget)
        self.gridLayout.setObjectName("gridLayout")
        self.Layout_graphs = QtWidgets.QGridLayout()
        self.Layout_graphs.setObjectName("Layout_graphs")
        
        self.plt = GraphicsLayoutWidget(self.centralwidget)
        self.pltB = GraphicsLayoutWidget(self.centralwidget)
        
        self.plt.setAutoFillBackground(False)
        self.plt.setStyleSheet("border: 0px;")
        self.plt.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.plt.setFrameShadow(QtWidgets.QFrame.Plain)
        self.plt.setLineWidth(0)
        self.plt.setObjectName("plt")
        
        self.pltB.setAutoFillBackground(False)
        self.pltB.setStyleSheet("border: 0px;")
        self.pltB.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.pltB.setFrameShadow(QtWidgets.QFrame.Plain)
        self.pltB.setLineWidth(0)
        self.pltB.setObjectName("pltB")
        
        self.label = QtWidgets.QLabel()
        self.Layout_graphs.addWidget(self.label, 0, 0, 1, 1)
        def link1(linkStr):
            QtGui.QDesktopServices.openUrl(QtCore.QUrl(linkStr))
        self.label.linkActivated.connect(link1)
        self.label.setText('<a href="https://openqcm.com/"> <font color=#000000 >Open-source Python application for displaying, processing and storing real-time data from openQCM Q-1 Device</font> </a>')
        
        self.Layout_graphs.addWidget(self.plt, 1, 0, 1, 1)
        self.Layout_graphs.addWidget(self.pltB, 2, 0, 1, 1)
        
        self.gridLayout.addLayout(self.Layout_graphs, 3, 1, 1, 1)
        MainWindow2.setCentralWidget(self.centralwidget)

        self.retranslateUi(MainWindow2)
        QtCore.QMetaObject.connectSlotsByName(MainWindow2)

    def retranslateUi(self, MainWindow2):
        _translate = QtCore.QCoreApplication.translate
        MainWindow2.setWindowIcon(QtGui.QIcon('favicon.ico')) #.png
        MainWindow2.setWindowTitle(_translate("MainWindow2", "OPENQCM Q-1 - Real-Time Plot GUI 2.1"))


############################################################################################################

class Ui_Info(object):
    def setupUi(self, MainWindow3):
        #MainWindow3.setObjectName("MainWindow3")
        #MainWindow3.setGeometry(500, 50, 100, 500)
        #MainWindow3.setFixedSize(100, 500)
        #MainWindow3.resize(100, 500)
        #MainWindow3.setMinimumSize(QtCore.QSize(0, 0))
        MainWindow3.setStyleSheet("")
        MainWindow3.setTabShape(QtWidgets.QTabWidget.Rounded) 
        MainWindow3.setMinimumSize(QtCore.QSize(268, 518))
        self.centralwidget = QtWidgets.QWidget(MainWindow3)

        self.centralwidget.setObjectName("centralwidget")
        self.gridLayout = QtWidgets.QGridLayout(self.centralwidget)
        self.gridLayout.setObjectName("gridLayout")
        self.gridLayout_2 = QtWidgets.QGridLayout()
        self.gridLayout_2.setObjectName("gridLayout_2")
        
        # Setup Information -------------------------------------------------------------------
        self.info1 = QtWidgets.QLabel()
        self.info1.setStyleSheet('background: #008EC0; padding: 1px;')
        self.info1.setText("<font color=#ffffff > Setup Information&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</font>")
        #self.info1.setFixedWidth(250)
        self.info1.setFixedHeight(15)
        self.gridLayout_2.addWidget(self.info1, 0, 0, 1, 1)
        
        # Device Setup -------------------------------------------------------------------------
        self.info1a = QtWidgets.QLabel()
        self.info1a.setStyleSheet('background: white; padding: 1px; border: 1px solid #cccccc')
        self.info1a.setText("<font color=#0000ff > Device Setup</font>")
        #self.info1a.setFixedWidth(250)
        #self.info1a.setFixedHeight(22)
        self.gridLayout_2.addWidget(self.info1a, 1, 0, 1, 1)
        
        # Operation Mode -----------------------------------------------------------------------
        self.info11 = QtWidgets.QLabel()
        self.info11.setStyleSheet('background: white; padding: 1px; border: 1px solid #cccccc')
        self.info11.setText("<font color=#0000ff > Operation Mode </font>")
        #self.info11.setFixedWidth(250)
        self.gridLayout_2.addWidget(self.info11, 2, 0, 1, 1)
        
        # Data Information ---------------------------------------------------------------------
        self.info = QtWidgets.QLabel()
        self.info.setStyleSheet('background: #008EC0; padding: 1px;')
        self.info.setText("<font color=#ffffff > Data Information&nbsp;</font>")
        #self.info.setFixedWidth(250)
        self.info.setFixedHeight(15)                  
        self.gridLayout_2.addWidget(self.info, 3, 0, 1, 1)
        
        # Selected Frequency -------------------------------------------------------------------
        self.info2 = QtWidgets.QLabel()
        self.info2.setStyleSheet('background: white; padding: 1px; border: 1px solid #cccccc')
        self.info2.setText("<font color=#0000ff > Selected Frequency </font>")
        #self.info2.setFixedWidth(250)                  
        self.gridLayout_2.addWidget(self.info2, 4, 0, 1, 1)
        
        # Frequency Value ----------------------------------------------------------------------
        self.info6 = QtWidgets.QLabel()
        self.info6.setStyleSheet('background: white; padding: 1px; border: 1px solid #cccccc')
        self.info6.setText("<font color=#0000ff > Frequency Value </font>")
        #self.info6.setFixedWidth(250)                  
        self.gridLayout_2.addWidget(self.info6, 5, 0, 1, 1)
        
        # Start Frequency ----------------------------------------------------------------------
        self.info3 = QtWidgets.QLabel()
        self.info3.setStyleSheet('background: white; padding: 1px; border: 1px solid #cccccc')
        self.info3.setText("<font color=#0000ff > Start Frequency </font>")
        #self.info3.setFixedWidth(250)                   
        self.gridLayout_2.addWidget(self.info3, 6, 0, 1, 1)
        
        # Stop Frequency -----------------------------------------------------------------------
        self.info4 = QtWidgets.QLabel()
        self.info4.setStyleSheet('background: white; padding: 1px; border: 1px solid #cccccc')
        self.info4.setText("<font color=#0000ff > Stop Frequency </font>")
        #self.info4.setFixedWidth(250)                   
        self.gridLayout_2.addWidget(self.info4, 7, 0, 1, 1)
        
        # Frequency Range----------------------------------------------------------------------
        self.info4a = QtWidgets.QLabel()
        self.info4a.setStyleSheet('background: white; padding: 1px; border: 1px solid #cccccc')
        self.info4a.setText("<font color=#0000ff > Frequency Range </font>")
        #self.info4a.setFixedWidth(250)                   
        self.gridLayout_2.addWidget(self.info4a, 8, 0, 1, 1)
        
        # Sample Rate----------------------------------------------------------------------
        self.info5 = QtWidgets.QLabel()
        self.info5.setStyleSheet('background: white; padding: 1px; border: 1px solid #cccccc')
        self.info5.setText("<font color=#0000ff > Sample Rate </font>")
        #self.info5.setFixedWidth(250)                   
        self.gridLayout_2.addWidget(self.info5, 9, 0, 1, 1)
        
        # Sample Number----------------------------------------------------------------------
        self.info7 = QtWidgets.QLabel()
        self.info7.setStyleSheet('background: white; padding: 1px; border: 1px solid #cccccc')
        self.info7.setText("<font color=#0000ff > Sample Number </font>")
        #self.info7.setFixedWidth(250)                 
        self.gridLayout_2.addWidget(self.info7, 10, 0, 1, 1)
        
        # Reference Settings -------------------------------------------------------------------
        self.inforef = QtWidgets.QLabel()
        self.inforef.setStyleSheet('background: #008EC0; padding: 1px;')
        #self.inforef1.setAlignment(QtCore.Qt.AlignCenter) 
        self.inforef.setText("<font color=#ffffff > Reference Settings </font>")
        self.inforef.setFixedHeight(15)
        #self.inforef.setFixedWidth(250)                      
        self.gridLayout_2.addWidget(self.inforef, 11, 0, 1, 1)
        
        # Ref. Frequency -----------------------------------------------------------------------
        self.inforef1 = QtWidgets.QLabel()
        self.inforef1.setStyleSheet('background: white; padding: 1px; border: 1px solid #cccccc')
        #self.inforef1.setAlignment(QtCore.Qt.AlignCenter) 
        self.inforef1.setText("<font color=#0000ff > Ref. Frequency </font>")
        #self.inforef1.setFixedWidth(250)
        self.gridLayout_2.addWidget(self.inforef1, 12, 0, 1, 1)
        # Ref. Dissipation -----------------------------------------------------------------------
        
        self.inforef2 = QtWidgets.QLabel()
        self.inforef2.setStyleSheet('background: white; padding: 1px; border: 1px solid #cccccc')
        #self.inforef2.setAlignment(QtCore.Qt.AlignCenter) 
        self.inforef2.setText("<font color=#0000ff > Ref. Dissipation </font>")
        #self.inforef2.setFixedWidth(250)                     
        self.gridLayout_2.addWidget(self.inforef2, 13, 0, 1, 1)
        
        # Current Data ---------------------------------------------------------------------------
        self.l8 = QtWidgets.QLabel()
        self.l8.setStyleSheet('background: #008EC0; padding: 1px;')
        self.l8.setText("<font color=#ffffff > Current Data </font>")
        self.l8.setFixedHeight(15)                
        #self.l8.setFixedWidth(250)
        self.gridLayout_2.addWidget(self.l8, 14, 0, 1, 1) 
        
        # Resonance Frequency -------------------------------------------------------------------
        self.l7 = QtWidgets.QLabel()
        self.l7.setStyleSheet('background: white; padding: 1px; border: 1px solid #cccccc')
        self.l7.setText("<font color=#0000ff >  Resonance Frequency </font>")
        #self.l7.setAlignment(QtCore.Qt.AlignCenter) 
        #self.l7.setFixedWidth(250)
        self.gridLayout_2.addWidget(self.l7, 15, 0, 1, 1) 
        
        # Dissipation ---------------------------------------------------------------------------
        self.l6 = QtWidgets.QLabel()
        self.l6.setStyleSheet('background: white; padding: 1px; border: 1px solid #cccccc')
        self.l6.setText("<font color=#0000ff > Dissipation  </font>")
        #self.l6.setFixedWidth(250)                
        self.gridLayout_2.addWidget(self.l6, 16, 0, 1, 1) 
        
        # Temperature ---------------------------------------------------------------------------
        self.l6a = QtWidgets.QLabel()
        self.l6a.setStyleSheet('background: white; padding: 1px; border: 1px solid #cccccc')
        self.l6a.setText("<font color=#0000ff >  Temperature </font>")
        #self.l6a.setFixedWidth(250)                
        self.gridLayout_2.addWidget(self.l6a, 17, 0, 1, 1)
     
        # Info from openQCM website -------------------------------------------------------------
        self.lweb = QtWidgets.QLabel()
        self.lweb.setStyleSheet('background: #008EC0; padding: 1px;')
        self.lweb.setText("<font color=#ffffff > Info from openQCM Website </font>")
        self.lweb.setFixedHeight(15)                
        #self.lweb.setFixedWidth(250)
        self.gridLayout_2.addWidget(self.lweb, 18, 0, 1, 1) 
        
        # Check internet connection -------------------------------------------------------------
        self.lweb2 = QtWidgets.QLabel()
        self.lweb2.setStyleSheet('background: white; padding: 1px; border: 1px solid #cccccc')
        self.lweb2.setText("<font color=#0000ff > Checking your internet connection </font>")
        # self.lweb2.setFixedHeight(20)                
        #self.lweb2.setFixedWidth(250)
        self.gridLayout_2.addWidget(self.lweb2, 19, 0, 1, 1)   
        
        # Software update status ----------------------------------------------------------------
        self.lweb3 = QtWidgets.QLabel()
        self.lweb3.setStyleSheet('background: white; padding: 1px; border: 1px solid #cccccc')
        self.lweb3.setText("<font color=#0000ff > Software update status </font>")
        # self.lweb3.setFixedHeight(16)                
        #self.lweb3.setFixedWidth(300)
        self.gridLayout_2.addWidget(self.lweb3, 20, 0, 1, 1) 
        
        # dovnload button -----------------------------------------------------------------------
        self.pButton_Download = QtWidgets.QPushButton(self.centralwidget)
        self.pButton_Download.setIcon(QtGui.QIcon(QtGui.QPixmap("download_icon.ico"))) #.png
        self.pButton_Download.setMinimumSize(QtCore.QSize(0, 0))
        self.pButton_Download.setObjectName("pButton_download")
        self.pButton_Download.setFixedWidth(145)
        self.pButton_Download.setEnabled(False)
        self.gridLayout_2.addWidget(self.pButton_Download, 21, 0, 1, 1, QtCore.Qt.AlignRight)
        ##########################################################################################
        
        self.gridLayout.addLayout(self.gridLayout_2, 3, 1, 1, 1)
        MainWindow3.setCentralWidget(self.centralwidget) 

        self.retranslateUi(MainWindow3)
        QtCore.QMetaObject.connectSlotsByName(MainWindow3)

    def retranslateUi(self, MainWindow3):
        _translate = QtCore.QCoreApplication.translate
        self.pButton_Download.setText(_translate("MainWindow3", " Download zip file"))
        MainWindow3.setWindowIcon(QtGui.QIcon('favicon.ico')) #.png 
        MainWindow3.setWindowTitle(_translate("MainWindow3", "Information GUI"))


