
from PyQt5 import QtCore, QtGui, QtWidgets

class Ui_Info(object):
    def setupUi(self, MainWindow):
        MainWindow.setObjectName("MainWindow")
        MainWindow.resize(100, 500)
        MainWindow.setMinimumSize(QtCore.QSize(0, 0))
        MainWindow.setStyleSheet("")
        MainWindow.setTabShape(QtWidgets.QTabWidget.Rounded)
        self.centralwidget = QtWidgets.QWidget(MainWindow)
        self.centralwidget.setObjectName("centralwidget")
        self.gridLayout = QtWidgets.QGridLayout(self.centralwidget)
        self.gridLayout.setObjectName("gridLayout")
        self.gridLayout_2 = QtWidgets.QGridLayout()
        self.gridLayout_2.setObjectName("gridLayout_2")
        
        # Setup Information -------------------------------------------------------------------
        self.info1 = QtWidgets.QLabel()
        self.info1.setStyleSheet('background: #008EC0; padding: 1px;')
        self.info1.setText("<font color=#ffffff > Setup Information&nbsp;</font>")
        self.info1.setFixedWidth(250)
        self.info1.setFixedHeight(15)
        self.gridLayout_2.addWidget(self.info1, 0, 0, 1, 1)
        
        # Device Setup -------------------------------------------------------------------------
        self.info1a = QtWidgets.QLabel()
        self.info1a.setStyleSheet('background: white; padding: 1px; border: 1px solid #cccccc')
        self.info1a.setText("<font color=#0000ff > Device Setup</font>")
        self.info1a.setFixedWidth(250)
        #self.info1a.setFixedHeight(22)
        self.gridLayout_2.addWidget(self.info1a, 1, 0, 1, 1)
        
        # Operation Mode -----------------------------------------------------------------------
        self.info11 = QtWidgets.QLabel()
        self.info11.setStyleSheet('background: white; padding: 1px; border: 1px solid #cccccc')
        self.info11.setText("<font color=#0000ff > Operation Mode </font>")
        self.info11.setFixedWidth(250)
        self.gridLayout_2.addWidget(self.info11, 2, 0, 1, 1)
        
        # Data Information ---------------------------------------------------------------------
        self.info = QtWidgets.QLabel()
        self.info.setStyleSheet('background: #008EC0; padding: 1px;')
        self.info.setText("<font color=#ffffff > Data Information&nbsp;</font>")
        self.info.setFixedWidth(250)
        self.info.setFixedHeight(15)                  
        self.gridLayout_2.addWidget(self.info, 3, 0, 1, 1)
        
        # Selected Frequency -------------------------------------------------------------------
        self.info2 = QtWidgets.QLabel()
        self.info2.setStyleSheet('background: white; padding: 1px; border: 1px solid #cccccc')
        self.info2.setText("<font color=#0000ff > Selected Frequency </font>")
        self.info2.setFixedWidth(250)                  
        self.gridLayout_2.addWidget(self.info2, 4, 0, 1, 1)
        
        # Frequency Value ----------------------------------------------------------------------
        self.info6 = QtWidgets.QLabel()
        self.info6.setStyleSheet('background: white; padding: 1px; border: 1px solid #cccccc')
        self.info6.setText("<font color=#0000ff > Frequency Value </font>")
        self.info6.setFixedWidth(250)                  
        self.gridLayout_2.addWidget(self.info6, 5, 0, 1, 1)
        
        # Start Frequency ----------------------------------------------------------------------
        self.info3 = QtWidgets.QLabel()
        self.info3.setStyleSheet('background: white; padding: 1px; border: 1px solid #cccccc')
        self.info3.setText("<font color=#0000ff > Start Frequency </font>")
        self.info3.setFixedWidth(250)                   
        self.gridLayout_2.addWidget(self.info3, 6, 0, 1, 1)
        
        # Stop Frequency -----------------------------------------------------------------------
        self.info4 = QtWidgets.QLabel()
        self.info4.setStyleSheet('background: white; padding: 1px; border: 1px solid #cccccc')
        self.info4.setText("<font color=#0000ff > Stop Frequency </font>")
        self.info4.setFixedWidth(250)                   
        self.gridLayout_2.addWidget(self.info4, 7, 0, 1, 1)
        
        # Frequency Range----------------------------------------------------------------------
        self.info4a = QtWidgets.QLabel()
        self.info4a.setStyleSheet('background: white; padding: 1px; border: 1px solid #cccccc')
        self.info4a.setText("<font color=#0000ff > Frequency Range </font>")
        self.info4a.setFixedWidth(250)                   
        self.gridLayout_2.addWidget(self.info4a, 8, 0, 1, 1)
        
        # Sample Rate----------------------------------------------------------------------
        self.info5 = QtWidgets.QLabel()
        self.info5.setStyleSheet('background: white; padding: 1px; border: 1px solid #cccccc')
        self.info5.setText("<font color=#0000ff > Sample Rate </font>")
        self.info5.setFixedWidth(250)                   
        self.gridLayout_2.addWidget(self.info5, 9, 0, 1, 1)
        
        # Sample Number----------------------------------------------------------------------
        self.info7 = QtWidgets.QLabel()
        self.info7.setStyleSheet('background: white; padding: 1px; border: 1px solid #cccccc')
        self.info7.setText("<font color=#0000ff > Sample Number </font>")
        self.info7.setFixedWidth(250)                 
        self.gridLayout_2.addWidget(self.info7, 10, 0, 1, 1)
        
        # Reference Settings -------------------------------------------------------------------
        self.inforef = QtWidgets.QLabel()
        self.inforef.setStyleSheet('background: #008EC0; padding: 1px;')
        #self.inforef1.setAlignment(QtCore.Qt.AlignCenter) 
        self.inforef.setText("<font color=#ffffff > Reference Settings </font>")
        self.inforef.setFixedHeight(15)
        self.inforef.setFixedWidth(250)                      
        self.gridLayout_2.addWidget(self.inforef, 11, 0, 1, 1)
        
        # Ref. Frequency -----------------------------------------------------------------------
        self.inforef1 = QtWidgets.QLabel()
        self.inforef1.setStyleSheet('background: white; padding: 1px; border: 1px solid #cccccc')
        #self.inforef1.setAlignment(QtCore.Qt.AlignCenter) 
        self.inforef1.setText("<font color=#0000ff > Ref. Frequency </font>")
        self.inforef1.setFixedWidth(250)
        self.gridLayout_2.addWidget(self.inforef1, 12, 0, 1, 1)
        # Ref. Dissipation -----------------------------------------------------------------------
        
        self.inforef2 = QtWidgets.QLabel()
        self.inforef2.setStyleSheet('background: white; padding: 1px; border: 1px solid #cccccc')
        #self.inforef2.setAlignment(QtCore.Qt.AlignCenter) 
        self.inforef2.setText("<font color=#0000ff > Ref. Dissipation </font>")
        self.inforef2.setFixedWidth(250)                     
        self.gridLayout_2.addWidget(self.inforef2, 13, 0, 1, 1)
        
        # Current Data ---------------------------------------------------------------------------
        self.l8 = QtWidgets.QLabel()
        self.l8.setStyleSheet('background: #008EC0; padding: 1px;')
        self.l8.setText("<font color=#ffffff > Current Data </font>")
        self.l8.setFixedHeight(15)                
        self.l8.setFixedWidth(250)
        self.gridLayout_2.addWidget(self.l8, 14, 0, 1, 1) 
        
        # Resonance Frequency -------------------------------------------------------------------
        self.l7 = QtWidgets.QLabel()
        self.l7.setStyleSheet('background: white; padding: 1px; border: 1px solid #cccccc')
        self.l7.setText("<font color=#0000ff >  Resonance Frequency </font>")
        #self.l7.setAlignment(QtCore.Qt.AlignCenter) 
        self.l7.setFixedWidth(250)
        self.gridLayout_2.addWidget(self.l7, 15, 0, 1, 1) 
        
        # Dissipation ---------------------------------------------------------------------------
        self.l6 = QtWidgets.QLabel()
        self.l6.setStyleSheet('background: white; padding: 1px; border: 1px solid #cccccc')
        self.l6.setText("<font color=#0000ff > Dissipation  </font>")
        self.l6.setFixedWidth(250)                
        self.gridLayout_2.addWidget(self.l6, 16, 0, 1, 1) 
        
        # Temperature ---------------------------------------------------------------------------
        self.l6a = QtWidgets.QLabel()
        self.l6a.setStyleSheet('background: white; padding: 1px; border: 1px solid #cccccc')
        self.l6a.setText("<font color=#0000ff >  Temperature </font>")
        self.l6a.setFixedWidth(250)                
        self.gridLayout_2.addWidget(self.l6a, 17, 0, 1, 1)
        #----------------------------------------------------------------------
        
        self.gridLayout.addLayout(self.gridLayout_2, 3, 1, 1, 1)
        MainWindow.setCentralWidget(self.centralwidget)

        self.retranslateUi(MainWindow)
        QtCore.QMetaObject.connectSlotsByName(MainWindow)

    def retranslateUi(self, MainWindow):
        _translate = QtCore.QCoreApplication.translate
        MainWindow.setWindowIcon(QtGui.QIcon('favicon.ico')) #.png
        MainWindow.setWindowTitle(_translate("MainWindow", "Information GUI"))

'''
if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    MainWindow = QtWidgets.QMainWindow()
    ui = Ui_Info()
    ui.setupUi(MainWindow)
    MainWindow.show()
    sys.exit(app.exec_())
'''
