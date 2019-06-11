#!/usr/bin/python3
# -*- coding: utf-8 -*-

import sys
# import numpy as np
import pandas as pd

from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure

from PyQt5.QtCore import QAbstractTableModel, Qt
from PyQt5.QtWidgets import QApplication, QMainWindow, QDialog
from gui import Ui_MainWindow
from filter import Ui_FilterDialog
from ivs.inout.fits import read_spectrum


class PandasModel(QAbstractTableModel):
    """
    Class to populate a table view with a pandas dataframe
    """
    def __init__(self, data, parent=None):
        QAbstractTableModel.__init__(self, parent)
        self._data = data

    def rowCount(self, parent=None):
        return self._data.shape[0]

    def columnCount(self, parent=None):
        return self._data.shape[1]

    def data(self, index, role=Qt.DisplayRole):
        if index.isValid():
            if role == Qt.DisplayRole:
                return str(self._data.iloc[index.row(), index.column()])
        return None

    def headerData(self, col, orientation, role):
        if orientation == Qt.Horizontal and role == Qt.DisplayRole:
            return self._data.columns[col]
        return None


class mainwindow_exec():

    def __init__(self):

        app = QApplication(sys.argv)
        self.window = QMainWindow()
        self.filter = QDialog()

        self.ui = Ui_MainWindow()
        self.ui.setupUi(self.window)
        self.filterui = Ui_FilterDialog()
        self.filterui.setupUi(self.filter)

        self.init_data(obs_only=True)
        self.add_mpl()
        self.window.show()

        self.ui.Button_filter.clicked.connect(lambda: self.filter_dialog())
        self.ui.Button_plot.clicked.connect(lambda: self.update_plot())

        sys.exit(app.exec_())

    def init_data(self, obs_only):

        # Load hermes tsv overview into pandas frame
        self.DATA = pd.read_csv('/STER/mercator/hermes/HermesFullDataOverview.tsv',
                                sep='\t', skiprows=2, header=None,
                                names=['unseq', 'prog_id', 'obsmode', 'bvcor',
                                       'observer', 'object', 'ra', 'dec',
                                       'bjd', 'exptime', 'pmtotal', 'date-avg',
                                       'airmass', 'filename'])

        # obsmode = HRF_OBJ
        if obs_only:
            self.DATA = self.DATA[self.DATA['obsmode'] == 'HRF_OBJ']

    def add_mpl(self):
        self.fig = Figure()
        self.ax = self.fig.add_subplot(111)
        self.canvas = FigureCanvas(self.fig)
        self.toolbar = NavigationToolbar(self.canvas, self.ui.matplotlib,
                                         coordinates=True)
        self.ui.mpl_layout.addWidget(self.toolbar)
        self.ui.mpl_layout.addWidget(self.canvas)
        self.canvas.draw()

    def update_plot(self):
        try:
            idx = self.filterui.FiltertableView.selectionModel().selectedRows()
            filename = self.model.data(self.model.index(idx[0].row(), 13))
            starname = self.model.data(self.model.index(idx[0].row(), 5))
            if 'raw' in filename:
                extname = 'HRF_OBJ_ext_CosmicsRemoved_log_merged_c'
                filename = filename.replace('raw', 'reduced')
                filename = filename.replace('HRF_OBJ', extname)
            wav, flux = read_spectrum(filename)
            self.ax.cla()
            self.ax.plot(wav, flux)
            self.ax.set_title(starname)
            self.canvas.draw()
        except(AttributeError, IndexError):
            pass
        except(FileNotFoundError):
            print('FileNotFoundError: \"' + filename + '\" not found...')

    def filter_dialog(self):
        self.model = PandasModel(self.DATA)
        self.filterui.FiltertableView.setModel(self.model)
        self.filter.show()


if __name__ == "__main__":
    mainwindow_exec()
