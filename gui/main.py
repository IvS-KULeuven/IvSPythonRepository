#!/usr/bin/python3
# -*- coding: utf-8 -*-

import sys
# import numpy as np
import pandas as pd

from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure
from PyQt5.QtWidgets import QApplication, QMainWindow
from PyQt5 import QtCore
from ivs.gui.gui import Ui_MainWindow
from ivs.inout.fits import read_spectrum
from ivs.gui.models import PandasModel, SelectionModel


class MainWindow(QMainWindow, Ui_MainWindow):
    """
    Main GUI window: used to interactively view HERMES spectra.

    Command line usage: $ python main.py

    Current functionality:
    - Use filter tab to select HERMES spectra
    - Use plot tab to display matplotlib plot of highlighted spectra

    TODO:
    - I/O functions for subselection
    - Implement more filter functions for overview table (RA,DEC,TIME)
    - Add FITS header/stellar/etc info into new panel below MPL canvas
    - Integrate with SIMBAD (for search, expanded stellar info, etc)
    - Integrate fitting/science functions into GUI
      (i.e. interactive lambda selection for profile fitting, etc)
    """
    def __init__(self, parent=None):
        # Initialise the Main window
        QMainWindow.__init__(self, parent=parent)
        self.ui = Ui_MainWindow()
        self.ui.setupUi(self)

        self.selected_data = []

        # Load HERMES overview TSV file
        self.load_overview(obs_only=True)
        # Initialise Filter/Data overview
        self.filtertab(self.overview)
        # Initialise MPL figure and toolbar
        self.add_mpl()

        # Connect Main window buttons to relevant functions
        self.ui.FilterAddButton.clicked.connect(lambda: self.filter_add())
        self.ui.PlotButton.clicked.connect(lambda: self.update_plot())
        self.ui.SearchPushButton.clicked.connect(lambda: self.Objectfilter())
        self.ui.Selectedpushbutton.clicked.connect(lambda: self.Remove_selected())

    def keyPressEvent(self, event):

        # Delete key used to remove highlighted spectra from selected sidebar
        if event.key() == QtCore.Qt.Key_Delete:
            self.Remove_selected()

    def Remove_selected(self):
        # Get highlighted rows
        idx = self.ui.SelectedDataTable.selectionModel().selectedRows()

        # Retrieve filename corresponding to highlighted rows
        # ( TODO probably a nicer way to do this?)
        filenames = []
        for idxxx in idx:
            filenames.append(self.selectionmodel.data(
                             self.selectionmodel.index(idxxx.row(), 0)))

        # Iterate over highlighted files, remove rom selection list
        for filename in filenames:
            if len(self.selected_data):
                if filename in self.selected_data:
                    self.selected_data.remove(filename)
                    continue

        # Regenerate selection model and tableview
        self.selectionmodel = SelectionModel(self.selected_data)
        self.ui.SelectedDataTable.setModel(self.selectionmodel)

    def load_overview(self, obs_only):

        # Load hermes tsv overview into pandas frame
        self.overview = pd.read_csv('/STER/mercator/hermes/HermesFullDataOverview.tsv',
                                sep='\t', skiprows=2, header=None,
                                names=['unseq', 'prog_id', 'obsmode', 'bvcor',
                                       'observer', 'object', 'ra', 'dec',
                                       'bjd', 'exptime', 'pmtotal', 'date-avg',
                                       'airmass', 'filename'])

        # obsmode = HRF_OBJ
        if obs_only:
            self.overview = self.overview[self.overview['obsmode'] == 'HRF_OBJ']

    def add_mpl(self):
        # Initial MPL figure and toolbar
        self.fig = Figure()
        self.ax = self.fig.add_subplot(111)
        self.canvas = FigureCanvas(self.fig)
        self.toolbar = NavigationToolbar(self.canvas, self.ui.matplotlib,
                                         coordinates=True)
        # Add MPL canvas and toolbar to Main window
        self.ui.mpl_layout.addWidget(self.toolbar)
        self.ui.mpl_layout.addWidget(self.canvas)
        self.canvas.draw()

    def update_plot(self):

        try:
            # Retrieve first highlighted row from tableview
            idx = self.ui.SelectedDataTable.selectionModel().selectedRows()[0]
            self.fig.clf()
            self.ax = self.fig.add_subplot(111)
            filename = self.selectionmodel.data(
                        self.selectionmodel.index(idx.row(), 0))
            # Defaulting to 'log_merged_c' due to incomplete '_cf'
            if 'raw' in filename:
                extname = 'HRF_OBJ_ext_CosmicsRemoved_log_merged_c'
                filename = filename.replace('raw', 'reduced')
                filename = filename.replace('HRF_OBJ', extname)
            # Load wav and flux of corresponding spectrum
            wav, flux = read_spectrum(filename)

            # Refresh and replot figure
            self.ax = self.fig.add_subplot(111)
            self.ax.plot(wav, flux)
            self.ax.set_title(filename)
            self.canvas.draw()
        except(AttributeError, IndexError):
            pass
        except(FileNotFoundError):
            error_text = 'FileNotFoundError: \"' + filename + '\" not found...'
            self.ui.statusBar.showMessage(error_text, 3000)

    def filtertab(self, input):
        # Initialise overview table from pandas dataframe, and show the dialog
        self.model = PandasModel(input)
        self.ui.FiltertableView.setModel(self.model)

    def filter_add(self):
            idx = self.ui.FiltertableView.selectionModel().selectedRows()
            skiptotal = 0
            for idxxx in idx:
                filename = self.model.data(self.model.index(idxxx.row(), 13))

                if len(self.selected_data):
                    if filename in self.selected_data:
                        skiptotal += 1
                        continue
                self.selected_data.append(filename)
            if len(self.selected_data):
                self.selectionmodel = SelectionModel(self.selected_data)
                self.ui.SelectedDataTable.setModel(self.selectionmodel)

            if skiptotal:
                if skiptotal == len(idx):
                    statustext = str(skiptotal) + ' spectra already added...'
                else:
                    statustext = (str(skiptotal) + ' spectra already added. ' +
                                  str(len(idx)-skiptotal) +
                                  ' new spectra added...')
            elif len(idx):
                statustext = str(len(idx))+' spectra addded...'
            else:
                statustext = ''
            self.ui.statusBar.showMessage(statustext, 3000)

    def Objectfilter(self):
        newdata = self.overview
        if len(self.ui.ObjLineEdit.text()):
            newdata = newdata[newdata['object'].str.contains(
                           self.ui.ObjLineEdit.text(), na=False, case=False)]
        if len(self.ui.ObserverLineEdit.text()):
            newdata = newdata[newdata['observer'].str.contains(
                           self.ui.ObserverLineEdit.text(), na=False, case=False)]
        if self.ui.ProgSpinBox.value():
            newdata = newdata[newdata['prog_id'] == self.ui.ProgSpinBox.value()]
        self.filtertab(newdata)

if __name__ == "__main__":
    app = QApplication(sys.argv)
    w = MainWindow()
    w.show()
    sys.exit(app.exec_())
