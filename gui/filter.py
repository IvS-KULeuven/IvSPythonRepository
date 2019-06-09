# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'filter.ui'
#
# Created by: PyQt5 UI code generator 5.6
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets

class Ui_FilterDialog(object):
    def setupUi(self, FilterDialog):
        FilterDialog.setObjectName("FilterDialog")
        FilterDialog.resize(1130, 757)
        self.verticalLayout = QtWidgets.QVBoxLayout(FilterDialog)
        self.verticalLayout.setObjectName("verticalLayout")
        self.FiltertableView = QtWidgets.QTableView(FilterDialog)
        self.FiltertableView.setSizeAdjustPolicy(QtWidgets.QAbstractScrollArea.AdjustIgnored)
        self.FiltertableView.setEditTriggers(QtWidgets.QAbstractItemView.AnyKeyPressed|QtWidgets.QAbstractItemView.DoubleClicked|QtWidgets.QAbstractItemView.EditKeyPressed|QtWidgets.QAbstractItemView.SelectedClicked)
        self.FiltertableView.setAlternatingRowColors(False)
        self.FiltertableView.setSelectionMode(QtWidgets.QAbstractItemView.SingleSelection)
        self.FiltertableView.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectRows)
        self.FiltertableView.setObjectName("FiltertableView")
        self.FiltertableView.horizontalHeader().setSortIndicatorShown(False)
        self.FiltertableView.horizontalHeader().setStretchLastSection(True)
        self.verticalLayout.addWidget(self.FiltertableView)

        self.retranslateUi(FilterDialog)
        QtCore.QMetaObject.connectSlotsByName(FilterDialog)

    def retranslateUi(self, FilterDialog):
        _translate = QtCore.QCoreApplication.translate
        FilterDialog.setWindowTitle(_translate("FilterDialog", "Dialog"))

