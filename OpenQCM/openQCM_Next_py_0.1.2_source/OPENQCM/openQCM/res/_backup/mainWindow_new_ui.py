# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'mainWindow_new.ui'
#
# Created by: PyQt5 UI code generator 5.9.2
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets

class Ui_MainWindow(object):
    def setupUi(self, MainWindow):
        MainWindow.setObjectName("MainWindow")
        MainWindow.resize(1297, 720)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(MainWindow.sizePolicy().hasHeightForWidth())
        MainWindow.setSizePolicy(sizePolicy)
        MainWindow.setMinimumSize(QtCore.QSize(1297, 720))
        MainWindow.setMaximumSize(QtCore.QSize(16777215, 16777215))
        icon = QtGui.QIcon()
        icon.addPixmap(QtGui.QPixmap("favicon.png"), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        MainWindow.setWindowIcon(icon)
        MainWindow.setStyleSheet("")
        MainWindow.setTabShape(QtWidgets.QTabWidget.Rounded)
        self.centralwidget = QtWidgets.QWidget(MainWindow)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.centralwidget.sizePolicy().hasHeightForWidth())
        self.centralwidget.setSizePolicy(sizePolicy)
        self.centralwidget.setMinimumSize(QtCore.QSize(1056, 0))
        self.centralwidget.setObjectName("centralwidget")
        self.layoutWidget = QtWidgets.QWidget(self.centralwidget)
        self.layoutWidget.setGeometry(QtCore.QRect(380, 10, 731, 631))
        self.layoutWidget.setObjectName("layoutWidget")
        self.Layout_graphs = QtWidgets.QGridLayout(self.layoutWidget)
        self.Layout_graphs.setContentsMargins(0, 0, 0, 0)
        self.Layout_graphs.setObjectName("Layout_graphs")
        self.pltB = GraphicsLayoutWidget(self.layoutWidget)
        self.pltB.setObjectName("pltB")
        self.Layout_graphs.addWidget(self.pltB, 0, 0, 1, 1)
        self.pltD = GraphicsLayoutWidget(self.layoutWidget)
        self.pltD.setObjectName("pltD")
        self.Layout_graphs.addWidget(self.pltD, 1, 0, 1, 1)
        self.verticalLayoutWidget = QtWidgets.QWidget(self.centralwidget)
        self.verticalLayoutWidget.setGeometry(QtCore.QRect(10, 330, 361, 311))
        self.verticalLayoutWidget.setObjectName("verticalLayoutWidget")
        self.horizontalLayout = QtWidgets.QHBoxLayout(self.verticalLayoutWidget)
        self.horizontalLayout.setContentsMargins(0, 0, 0, 0)
        self.horizontalLayout.setObjectName("horizontalLayout")
        self.plt = GraphicsLayoutWidget(self.verticalLayoutWidget)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.plt.sizePolicy().hasHeightForWidth())
        self.plt.setSizePolicy(sizePolicy)
        self.plt.setAutoFillBackground(False)
        self.plt.setStyleSheet("border: 0px;")
        self.plt.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.plt.setFrameShadow(QtWidgets.QFrame.Plain)
        self.plt.setLineWidth(0)
        self.plt.setObjectName("plt")
        self.horizontalLayout.addWidget(self.plt)
        self.gridLayoutWidget = QtWidgets.QWidget(self.centralwidget)
        self.gridLayoutWidget.setGeometry(QtCore.QRect(10, 10, 361, 312))
        self.gridLayoutWidget.setObjectName("gridLayoutWidget")
        self.gridLayout = QtWidgets.QGridLayout(self.gridLayoutWidget)
        self.gridLayout.setContentsMargins(0, 0, 0, 0)
        self.gridLayout.setObjectName("gridLayout")
        self.pButton_Tswitch_ON = QtWidgets.QPushButton(self.gridLayoutWidget)
        self.pButton_Tswitch_ON.setObjectName("pButton_Tswitch_ON")
        self.gridLayout.addWidget(self.pButton_Tswitch_ON, 5, 0, 1, 1)
        self.label_Temperature = QtWidgets.QLabel(self.gridLayoutWidget)
        self.label_Temperature.setObjectName("label_Temperature")
        self.gridLayout.addWidget(self.label_Temperature, 12, 0, 1, 1)
        self.doubleSpinBox_Temperature = QtWidgets.QDoubleSpinBox(self.gridLayoutWidget)
        self.doubleSpinBox_Temperature.setDecimals(0)
        self.doubleSpinBox_Temperature.setMinimum(5.0)
        self.doubleSpinBox_Temperature.setMaximum(45.0)
        self.doubleSpinBox_Temperature.setProperty("value", 25.0)
        self.doubleSpinBox_Temperature.setObjectName("doubleSpinBox_Temperature")
        self.gridLayout.addWidget(self.doubleSpinBox_Temperature, 11, 1, 1, 1)
        self.cBox_Port = QtWidgets.QComboBox(self.gridLayoutWidget)
        self.cBox_Port.setEditable(True)
        self.cBox_Port.setObjectName("cBox_Port")
        self.gridLayout.addWidget(self.cBox_Port, 0, 1, 1, 1)
        self.cBox_Speed = QtWidgets.QComboBox(self.gridLayoutWidget)
        self.cBox_Speed.setEditable(True)
        self.cBox_Speed.setObjectName("cBox_Speed")
        self.gridLayout.addWidget(self.cBox_Speed, 2, 1, 1, 1)
        self.l1 = QtWidgets.QLabel(self.gridLayoutWidget)
        self.l1.setObjectName("l1")
        self.gridLayout.addWidget(self.l1, 0, 0, 1, 1)
        self.pButton_Temperature_Set = QtWidgets.QPushButton(self.gridLayoutWidget)
        self.pButton_Temperature_Set.setObjectName("pButton_Temperature_Set")
        self.gridLayout.addWidget(self.pButton_Temperature_Set, 11, 0, 1, 1)
        self.cBox_Source = QtWidgets.QComboBox(self.gridLayoutWidget)
        self.cBox_Source.setObjectName("cBox_Source")
        self.gridLayout.addWidget(self.cBox_Source, 1, 1, 1, 1)
        self.info11 = QtWidgets.QLabel(self.gridLayoutWidget)
        self.info11.setObjectName("info11")
        self.gridLayout.addWidget(self.info11, 1, 0, 1, 1)
        self.l2 = QtWidgets.QLabel(self.gridLayoutWidget)
        self.l2.setObjectName("l2")
        self.gridLayout.addWidget(self.l2, 2, 0, 1, 1)
        self.label_D_Share = QtWidgets.QLabel(self.gridLayoutWidget)
        self.label_D_Share.setObjectName("label_D_Share")
        self.gridLayout.addWidget(self.label_D_Share, 10, 0, 1, 1)
        self.spinBox_Cycling_Time = QtWidgets.QSpinBox(self.gridLayoutWidget)
        self.spinBox_Cycling_Time.setMinimum(1)
        self.spinBox_Cycling_Time.setMaximum(1000)
        self.spinBox_Cycling_Time.setProperty("value", 50)
        self.spinBox_Cycling_Time.setObjectName("spinBox_Cycling_Time")
        self.gridLayout.addWidget(self.spinBox_Cycling_Time, 7, 1, 1, 1)
        self.pButton_Tswitch_OFF = QtWidgets.QPushButton(self.gridLayoutWidget)
        self.pButton_Tswitch_OFF.setObjectName("pButton_Tswitch_OFF")
        self.gridLayout.addWidget(self.pButton_Tswitch_OFF, 5, 1, 1, 1)
        self.label_P_Share = QtWidgets.QLabel(self.gridLayoutWidget)
        self.label_P_Share.setObjectName("label_P_Share")
        self.gridLayout.addWidget(self.label_P_Share, 8, 0, 1, 1)
        self.label_Cycling_Time = QtWidgets.QLabel(self.gridLayoutWidget)
        self.label_Cycling_Time.setObjectName("label_Cycling_Time")
        self.gridLayout.addWidget(self.label_Cycling_Time, 7, 0, 1, 1)
        self.pButton_PID_Set = QtWidgets.QPushButton(self.gridLayoutWidget)
        self.pButton_PID_Set.setObjectName("pButton_PID_Set")
        self.gridLayout.addWidget(self.pButton_PID_Set, 6, 0, 1, 1)
        self.spinBox_P_Share = QtWidgets.QSpinBox(self.gridLayoutWidget)
        self.spinBox_P_Share.setMaximum(100000)
        self.spinBox_P_Share.setProperty("value", 1000)
        self.spinBox_P_Share.setObjectName("spinBox_P_Share")
        self.gridLayout.addWidget(self.spinBox_P_Share, 8, 1, 1, 1)
        self.label_I_Share = QtWidgets.QLabel(self.gridLayoutWidget)
        self.label_I_Share.setObjectName("label_I_Share")
        self.gridLayout.addWidget(self.label_I_Share, 9, 0, 1, 1)
        self.label_Temperature_state = QtWidgets.QLabel(self.gridLayoutWidget)
        self.label_Temperature_state.setObjectName("label_Temperature_state")
        self.gridLayout.addWidget(self.label_Temperature_state, 4, 0, 1, 2)
        self.spinBox_I_Share = QtWidgets.QSpinBox(self.gridLayoutWidget)
        self.spinBox_I_Share.setMaximum(100000)
        self.spinBox_I_Share.setSingleStep(0)
        self.spinBox_I_Share.setProperty("value", 200)
        self.spinBox_I_Share.setObjectName("spinBox_I_Share")
        self.gridLayout.addWidget(self.spinBox_I_Share, 9, 1, 1, 1)
        self.line = QtWidgets.QFrame(self.gridLayoutWidget)
        self.line.setFrameShape(QtWidgets.QFrame.HLine)
        self.line.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.line.setObjectName("line")
        self.gridLayout.addWidget(self.line, 3, 0, 1, 2)
        self.spinBox_D_Share = QtWidgets.QSpinBox(self.gridLayoutWidget)
        self.spinBox_D_Share.setMaximum(100000)
        self.spinBox_D_Share.setProperty("value", 100)
        self.spinBox_D_Share.setObjectName("spinBox_D_Share")
        self.gridLayout.addWidget(self.spinBox_D_Share, 10, 1, 1, 1)
        self.indicator_temperature = QtWidgets.QLabel(self.gridLayoutWidget)
        self.indicator_temperature.setObjectName("indicator_temperature")
        self.gridLayout.addWidget(self.indicator_temperature, 12, 1, 1, 1)
        self.cBox_PID = QtWidgets.QComboBox(self.gridLayoutWidget)
        self.cBox_PID.setObjectName("cBox_PID")
        self.gridLayout.addWidget(self.cBox_PID, 6, 1, 1, 1)
        self.line_2 = QtWidgets.QFrame(self.centralwidget)
        self.line_2.setGeometry(QtCore.QRect(10, 640, 1101, 16))
        self.line_2.setFrameShape(QtWidgets.QFrame.HLine)
        self.line_2.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.line_2.setObjectName("line_2")
        self.layoutWidget1 = QtWidgets.QWidget(self.centralwidget)
        self.layoutWidget1.setGeometry(QtCore.QRect(10, 650, 361, 43))
        self.layoutWidget1.setObjectName("layoutWidget1")
        self.verticalLayout = QtWidgets.QVBoxLayout(self.layoutWidget1)
        self.verticalLayout.setContentsMargins(0, 0, 0, 0)
        self.verticalLayout.setObjectName("verticalLayout")
        self.infobar = QtWidgets.QLabel(self.layoutWidget1)
        self.infobar.setObjectName("infobar")
        self.verticalLayout.addWidget(self.infobar)
        self.infostatus = QtWidgets.QLabel(self.layoutWidget1)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.infostatus.sizePolicy().hasHeightForWidth())
        self.infostatus.setSizePolicy(sizePolicy)
        self.infostatus.setObjectName("infostatus")
        self.verticalLayout.addWidget(self.infostatus)
        self.splitter = QtWidgets.QSplitter(self.centralwidget)
        self.splitter.setGeometry(QtCore.QRect(380, 660, 731, 31))
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.splitter.sizePolicy().hasHeightForWidth())
        self.splitter.setSizePolicy(sizePolicy)
        self.splitter.setOrientation(QtCore.Qt.Horizontal)
        self.splitter.setObjectName("splitter")
        self.pButton_Start = QtWidgets.QPushButton(self.splitter)
        self.pButton_Start.setMinimumSize(QtCore.QSize(0, 0))
        self.pButton_Start.setObjectName("pButton_Start")
        self.pButton_Stop = QtWidgets.QPushButton(self.splitter)
        self.pButton_Stop.setObjectName("pButton_Stop")
        self.pButton_Reference = QtWidgets.QPushButton(self.splitter)
        self.pButton_Reference.setObjectName("pButton_Reference")
        self.pButton_Clear = QtWidgets.QPushButton(self.splitter)
        self.pButton_Clear.setObjectName("pButton_Clear")
        self.progressBar = QtWidgets.QProgressBar(self.splitter)
        self.progressBar.setProperty("value", 0)
        self.progressBar.setAlignment(QtCore.Qt.AlignCenter)
        self.progressBar.setTextVisible(True)
        self.progressBar.setObjectName("progressBar")
        self.gridLayoutWidget_2 = QtWidgets.QWidget(self.centralwidget)
        self.gridLayoutWidget_2.setGeometry(QtCore.QRect(1120, 10, 121, 131))
        self.gridLayoutWidget_2.setObjectName("gridLayoutWidget_2")
        self.gridLayout_F = QtWidgets.QGridLayout(self.gridLayoutWidget_2)
        self.gridLayout_F.setContentsMargins(0, 0, 0, 0)
        self.gridLayout_F.setObjectName("gridLayout_F")
        self.label_F3 = QtWidgets.QLabel(self.gridLayoutWidget_2)
        self.label_F3.setObjectName("label_F3")
        self.gridLayout_F.addWidget(self.label_F3, 2, 0, 1, 1)
        self.label_F_Title = QtWidgets.QLabel(self.gridLayoutWidget_2)
        self.label_F_Title.setObjectName("label_F_Title")
        self.gridLayout_F.addWidget(self.label_F_Title, 0, 0, 1, 2)
        self.F9 = QtWidgets.QLabel(self.gridLayoutWidget_2)
        self.F9.setObjectName("F9")
        self.gridLayout_F.addWidget(self.F9, 5, 1, 1, 1)
        self.F3 = QtWidgets.QLabel(self.gridLayoutWidget_2)
        self.F3.setObjectName("F3")
        self.gridLayout_F.addWidget(self.F3, 2, 1, 1, 1)
        self.label_F5 = QtWidgets.QLabel(self.gridLayoutWidget_2)
        self.label_F5.setObjectName("label_F5")
        self.gridLayout_F.addWidget(self.label_F5, 3, 0, 1, 1)
        self.F5 = QtWidgets.QLabel(self.gridLayoutWidget_2)
        self.F5.setObjectName("F5")
        self.gridLayout_F.addWidget(self.F5, 3, 1, 1, 1)
        self.label_7 = QtWidgets.QLabel(self.gridLayoutWidget_2)
        self.label_7.setObjectName("label_7")
        self.gridLayout_F.addWidget(self.label_7, 4, 0, 1, 1)
        self.F7 = QtWidgets.QLabel(self.gridLayoutWidget_2)
        self.F7.setObjectName("F7")
        self.gridLayout_F.addWidget(self.F7, 4, 1, 1, 1)
        self.label_F9 = QtWidgets.QLabel(self.gridLayoutWidget_2)
        self.label_F9.setObjectName("label_F9")
        self.gridLayout_F.addWidget(self.label_F9, 5, 0, 1, 1)
        self.label_F0 = QtWidgets.QLabel(self.gridLayoutWidget_2)
        self.label_F0.setObjectName("label_F0")
        self.gridLayout_F.addWidget(self.label_F0, 1, 0, 1, 1)
        self.F0 = QtWidgets.QLabel(self.gridLayoutWidget_2)
        self.F0.setObjectName("F0")
        self.gridLayout_F.addWidget(self.F0, 1, 1, 1, 1)
        self.gridLayoutWidget_3 = QtWidgets.QWidget(self.centralwidget)
        self.gridLayoutWidget_3.setGeometry(QtCore.QRect(1120, 330, 121, 131))
        self.gridLayoutWidget_3.setObjectName("gridLayoutWidget_3")
        self.gridLayout_3 = QtWidgets.QGridLayout(self.gridLayoutWidget_3)
        self.gridLayout_3.setContentsMargins(0, 0, 0, 0)
        self.gridLayout_3.setObjectName("gridLayout_3")
        self.D7 = QtWidgets.QLabel(self.gridLayoutWidget_3)
        self.D7.setObjectName("D7")
        self.gridLayout_3.addWidget(self.D7, 4, 1, 1, 1)
        self.label_D7 = QtWidgets.QLabel(self.gridLayoutWidget_3)
        self.label_D7.setObjectName("label_D7")
        self.gridLayout_3.addWidget(self.label_D7, 4, 0, 1, 1)
        self.label_D0 = QtWidgets.QLabel(self.gridLayoutWidget_3)
        self.label_D0.setObjectName("label_D0")
        self.gridLayout_3.addWidget(self.label_D0, 1, 0, 1, 1)
        self.label_D_Title = QtWidgets.QLabel(self.gridLayoutWidget_3)
        self.label_D_Title.setObjectName("label_D_Title")
        self.gridLayout_3.addWidget(self.label_D_Title, 0, 0, 1, 2)
        self.D3 = QtWidgets.QLabel(self.gridLayoutWidget_3)
        self.D3.setObjectName("D3")
        self.gridLayout_3.addWidget(self.D3, 2, 1, 1, 1)
        self.label_D3 = QtWidgets.QLabel(self.gridLayoutWidget_3)
        self.label_D3.setObjectName("label_D3")
        self.gridLayout_3.addWidget(self.label_D3, 2, 0, 1, 1)
        self.D0 = QtWidgets.QLabel(self.gridLayoutWidget_3)
        self.D0.setObjectName("D0")
        self.gridLayout_3.addWidget(self.D0, 1, 1, 1, 1)
        self.D5 = QtWidgets.QLabel(self.gridLayoutWidget_3)
        self.D5.setObjectName("D5")
        self.gridLayout_3.addWidget(self.D5, 3, 1, 1, 1)
        self.label_D5 = QtWidgets.QLabel(self.gridLayoutWidget_3)
        self.label_D5.setObjectName("label_D5")
        self.gridLayout_3.addWidget(self.label_D5, 3, 0, 1, 1)
        self.label_D9 = QtWidgets.QLabel(self.gridLayoutWidget_3)
        self.label_D9.setObjectName("label_D9")
        self.gridLayout_3.addWidget(self.label_D9, 5, 0, 1, 1)
        self.D9 = QtWidgets.QLabel(self.gridLayoutWidget_3)
        self.D9.setObjectName("D9")
        self.gridLayout_3.addWidget(self.D9, 5, 1, 1, 1)
        self.radioBtn_F0 = QtWidgets.QRadioButton(self.centralwidget)
        self.radioBtn_F0.setGeometry(QtCore.QRect(1120, 540, 85, 17))
        self.radioBtn_F0.setAutoExclusive(False)
        self.radioBtn_F0.setObjectName("radioBtn_F0")
        self.radioBtn_F3 = QtWidgets.QRadioButton(self.centralwidget)
        self.radioBtn_F3.setGeometry(QtCore.QRect(1120, 560, 88, 17))
        self.radioBtn_F3.setAutoExclusive(False)
        self.radioBtn_F3.setObjectName("radioBtn_F3")
        self.radioBtn_F5 = QtWidgets.QRadioButton(self.centralwidget)
        self.radioBtn_F5.setGeometry(QtCore.QRect(1120, 580, 88, 17))
        self.radioBtn_F5.setAutoExclusive(False)
        self.radioBtn_F5.setObjectName("radioBtn_F5")
        self.radioBtn_F7 = QtWidgets.QRadioButton(self.centralwidget)
        self.radioBtn_F7.setGeometry(QtCore.QRect(1120, 600, 88, 17))
        self.radioBtn_F7.setAutoExclusive(False)
        self.radioBtn_F7.setObjectName("radioBtn_F7")
        self.radioBtn_F9 = QtWidgets.QRadioButton(self.centralwidget)
        self.radioBtn_F9.setGeometry(QtCore.QRect(1120, 620, 88, 17))
        self.radioBtn_F9.setAutoExclusive(False)
        self.radioBtn_F9.setObjectName("radioBtn_F9")
        MainWindow.setCentralWidget(self.centralwidget)
        self.menuBar = QtWidgets.QMenuBar(MainWindow)
        self.menuBar.setGeometry(QtCore.QRect(0, 0, 1297, 21))
        self.menuBar.setObjectName("menuBar")
        self.menuMenu_Bar = QtWidgets.QMenu(self.menuBar)
        self.menuMenu_Bar.setObjectName("menuMenu_Bar")
        MainWindow.setMenuBar(self.menuBar)
        self.menuBar.addAction(self.menuMenu_Bar.menuAction())

        self.retranslateUi(MainWindow)
        QtCore.QMetaObject.connectSlotsByName(MainWindow)

    def retranslateUi(self, MainWindow):
        _translate = QtCore.QCoreApplication.translate
        MainWindow.setWindowTitle(_translate("MainWindow", "openQCM NEXT - version 0.1.2 (BETA-TEST)"))
        self.pButton_Tswitch_ON.setText(_translate("MainWindow", "Temperature Ctrl ON"))
        self.label_Temperature.setText(_translate("MainWindow", "Temperature (° C)"))
        self.l1.setText(_translate("MainWindow", "Serial COM Port"))
        self.pButton_Temperature_Set.setText(_translate("MainWindow", "Temperature Set"))
        self.info11.setText(_translate("MainWindow", "Operation mode"))
        self.l2.setText(_translate("MainWindow", "Frequency - Quartz Sensors"))
        self.label_D_Share.setText(_translate("MainWindow", "D Share [(mA*s)/K]"))
        self.pButton_Tswitch_OFF.setText(_translate("MainWindow", "Temperature Crtl OFF"))
        self.label_P_Share.setText(_translate("MainWindow", "P Share [mA/K]"))
        self.label_Cycling_Time.setText(_translate("MainWindow", "Cycling Time [msec]"))
        self.pButton_PID_Set.setText(_translate("MainWindow", "PID Set"))
        self.label_I_Share.setText(_translate("MainWindow", "I Share [mA/(K+sec)]"))
        self.label_Temperature_state.setText(_translate("MainWindow", "Temperature Control"))
        self.indicator_temperature.setText(_translate("MainWindow", "0"))
        self.infobar.setText(_translate("MainWindow", "Infobar"))
        self.infostatus.setText(_translate("MainWindow", "Program status "))
        self.pButton_Start.setText(_translate("MainWindow", "Start"))
        self.pButton_Stop.setText(_translate("MainWindow", "Stop"))
        self.pButton_Reference.setText(_translate("MainWindow", "Set/Reset Reference"))
        self.pButton_Clear.setText(_translate("MainWindow", "Clear Plots"))
        self.label_F3.setText(_translate("MainWindow", "F 3"))
        self.label_F_Title.setText(_translate("MainWindow", "Frequency (Hz)"))
        self.F9.setText(_translate("MainWindow", "0"))
        self.F3.setText(_translate("MainWindow", "0"))
        self.label_F5.setText(_translate("MainWindow", "F 5"))
        self.F5.setText(_translate("MainWindow", "0"))
        self.label_7.setText(_translate("MainWindow", "F 7"))
        self.F7.setText(_translate("MainWindow", "0"))
        self.label_F9.setText(_translate("MainWindow", "F 9"))
        self.label_F0.setText(_translate("MainWindow", "F 0"))
        self.F0.setText(_translate("MainWindow", "0"))
        self.D7.setText(_translate("MainWindow", "0"))
        self.label_D7.setText(_translate("MainWindow", "D 7"))
        self.label_D0.setText(_translate("MainWindow", "D 0"))
        self.label_D_Title.setText(_translate("MainWindow", "Dissipation (ppm)"))
        self.D3.setText(_translate("MainWindow", "0"))
        self.label_D3.setText(_translate("MainWindow", "D 3"))
        self.D0.setText(_translate("MainWindow", "0"))
        self.D5.setText(_translate("MainWindow", "0"))
        self.label_D5.setText(_translate("MainWindow", "D 5"))
        self.label_D9.setText(_translate("MainWindow", "D 9"))
        self.D9.setText(_translate("MainWindow", "0"))
        self.radioBtn_F0.setText(_translate("MainWindow", "Fundamental"))
        self.radioBtn_F3.setText(_translate("MainWindow", "3rd Overtone"))
        self.radioBtn_F5.setText(_translate("MainWindow", "5th Overtone"))
        self.radioBtn_F7.setText(_translate("MainWindow", "7th Overtone"))
        self.radioBtn_F9.setText(_translate("MainWindow", "9th Overtone"))
        self.menuMenu_Bar.setTitle(_translate("MainWindow", "Menu Bar"))

from pyqtgraph import GraphicsLayoutWidget