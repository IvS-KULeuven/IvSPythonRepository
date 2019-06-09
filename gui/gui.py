# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'gui.ui'
#
# Created by: PyQt5 UI code generator 5.6
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets

class Ui_MainWindow(object):
    def setupUi(self, MainWindow):
        MainWindow.setObjectName("MainWindow")
        MainWindow.resize(1097, 800)
        self.centralwidget = QtWidgets.QWidget(MainWindow)
        self.centralwidget.setAutoFillBackground(False)
        self.centralwidget.setObjectName("centralwidget")
        self.main_layout = QtWidgets.QGridLayout(self.centralwidget)
        self.main_layout.setObjectName("main_layout")
        self.buttonframe = QtWidgets.QFrame(self.centralwidget)
        self.buttonframe.setFrameShape(QtWidgets.QFrame.NoFrame)
        self.buttonframe.setFrameShadow(QtWidgets.QFrame.Plain)
        self.buttonframe.setObjectName("buttonframe")
        self.button_layout = QtWidgets.QHBoxLayout(self.buttonframe)
        self.button_layout.setObjectName("button_layout")
        self.Button_filter = QtWidgets.QPushButton(self.buttonframe)
        self.Button_filter.setObjectName("Button_filter")
        self.button_layout.addWidget(self.Button_filter)
        self.Button_plot = QtWidgets.QPushButton(self.buttonframe)
        self.Button_plot.setObjectName("Button_plot")
        self.button_layout.addWidget(self.Button_plot)
        self.main_layout.addWidget(self.buttonframe, 0, 0, 1, 1)
        self.matplotlib = QtWidgets.QWidget(self.centralwidget)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.MinimumExpanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.matplotlib.sizePolicy().hasHeightForWidth())
        self.matplotlib.setSizePolicy(sizePolicy)
        self.matplotlib.setAutoFillBackground(False)
        self.matplotlib.setObjectName("matplotlib")
        self.mpl_layout = QtWidgets.QVBoxLayout(self.matplotlib)
        self.mpl_layout.setContentsMargins(0, 0, 0, 0)
        self.mpl_layout.setObjectName("mpl_layout")
        self.main_layout.addWidget(self.matplotlib, 1, 0, 1, 1)
        MainWindow.setCentralWidget(self.centralwidget)
        self.menubar = QtWidgets.QMenuBar(MainWindow)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 1097, 29))
        self.menubar.setObjectName("menubar")
        self.menuFile = QtWidgets.QMenu(self.menubar)
        self.menuFile.setObjectName("menuFile")
        MainWindow.setMenuBar(self.menubar)
        self.statusbar = QtWidgets.QStatusBar(MainWindow)
        self.statusbar.setObjectName("statusbar")
        MainWindow.setStatusBar(self.statusbar)
        self.actionExit = QtWidgets.QAction(MainWindow)
        self.actionExit.setMenuRole(QtWidgets.QAction.TextHeuristicRole)
        self.actionExit.setObjectName("actionExit")
        self.actionLoad_Data = QtWidgets.QAction(MainWindow)
        self.actionLoad_Data.setObjectName("actionLoad_Data")
        self.menuFile.addAction(self.actionExit)
        self.menubar.addAction(self.menuFile.menuAction())

        self.retranslateUi(MainWindow)
        QtCore.QMetaObject.connectSlotsByName(MainWindow)

    def retranslateUi(self, MainWindow):
        _translate = QtCore.QCoreApplication.translate
        MainWindow.setWindowTitle(_translate("MainWindow", "MainWindow"))
        self.Button_filter.setText(_translate("MainWindow", "Filter..."))
        self.Button_plot.setText(_translate("MainWindow", "Plot"))
        self.menuFile.setTitle(_translate("MainWindow", "Fi&le"))
        self.actionExit.setText(_translate("MainWindow", "&Exit"))
        self.actionLoad_Data.setText(_translate("MainWindow", "&Load Data"))

