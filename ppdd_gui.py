import sys
from PyQt5.QtWidgets import (QMainWindow, QWidget, QTextEdit, QAction, QApplication,
                            QMessageBox, QGridLayout, QPushButton, QFileDialog, QErrorMessage,
                            QSizePolicy, QLabel, QLineEdit)
from PyQt5.QtGui import QIcon, QDoubleValidator
import matplotlib
# Make sure that we are using QT5
matplotlib.use('Qt5Agg')
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np

import ppdd

class PPDDWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.initUI()
        self.filename = ''
        self.ppdd = ppdd.PPDD()
        self.spectrum_vmax = 1e6        #used to display the amplitude spectrum
        
    def initUI(self): 
        self.layoutWidget = PPDDLayout()
        self.setCentralWidget(self.layoutWidget)

        openFile = QAction(QIcon.fromTheme('document-open'), 'Open', self)
        openFile.setShortcut('Ctrl+O')
        openFile.setStatusTip('Open File')
        openFile.triggered.connect(self.loadFileDialog)

        run = QAction(QIcon.fromTheme('media-playback-start'), 'Run', self)
        run.setShortcut('Ctrl+R')
        run.setStatusTip('Run ppdd')
        run.triggered.connect(self.runPPDD)

        exitAction = QAction(QIcon.fromTheme('application-exit'), 'Exit', self)
        exitAction.setShortcut('Ctrl+Q')
        exitAction.setStatusTip('Exit application')
        exitAction.triggered.connect(self.close)

        aboutAction = QAction(QIcon.fromTheme('help-about'), 'About', self)
        aboutAction.setStatusTip('About ppdd')
        aboutAction.triggered.connect(self.about)

        self.statusBar().showMessage('Ready')

        menubar = self.menuBar()
        fileMenu = menubar.addMenu('&File')
        fileMenu.addAction(openFile)
        fileMenu.addAction(run)
        fileMenu.addAction(exitAction)
        helpMenu = menubar.addMenu('&Help')
        helpMenu.addAction(aboutAction)

        toolbar = self.addToolBar('Open')
        toolbar.addAction(openFile)
        toolbar = self.addToolBar('Run')
        toolbar.addAction(run)
        toolbar = self.addToolBar('Exit')
        toolbar.addAction(exitAction)

        self.setWindowTitle('Python Plasma Density Diagnostics')    
        self.setWindowIcon(QIcon('ppdd.png'))   
        self.showMaximized()

    def closeEvent(self, event):
        reply = QMessageBox.question(self, 'Confirm to quit',
            "Are you sure to quit?", QMessageBox.Yes | QMessageBox.No, QMessageBox.No)
        if reply == QMessageBox.Yes:
            event.accept()
        else:
            event.ignore()        

    def loadFileDialog(self):
        filenames = QFileDialog.getOpenFileName(self, 'Open file', '', 'Data file (*.txt)', None, QFileDialog.DontUseNativeDialog)
        if filenames[0]:
            try :
                self.filename = filenames[0]
                self.ppdd.readfile(filenames[0])
                self.statusBar().showMessage('Scucessfully read file: {0}'.format(filenames[0]))
            except :
                error = QErrorMessage(self)
                error.showMessage('Error reading file: {0}'.format(filenames[0]))
                self.statusBar().showMessage('Error reading file: {0}'.format(filenames[0]))
                error.exec_() 

    def runPPDD(self):
        if not self.filename :      #check if a file is loaded
            self.statusBar().showMessage('No file is loaded, please load data file first!')
            return
        pyppd = self.ppdd
        if not pyppd.peak_fitted :  #if already fitted peaks, skip to speed up
            try :
                pyppd.find_peaks()
            except RuntimeError :
                self.statusBar().showMessage('Failed to find the secondary peak.')
                return

        self.layoutWidget.conf.fx.setText('{0:.6f}'.format(pyppd.fx))
        self.layoutWidget.conf.fy.setText('{0:.6f}'.format(pyppd.fy))
                
        #plot raw data and selected region
        ax_raw = self.layoutWidget.display.raw.axes
        ax_raw.set_title('Raw data')
        im1 = ax_raw.pcolormesh(pyppd.rawdata)
        rect1 = patches.Rectangle((pyppd.xmin, pyppd.ymin), pyppd.xmax-pyppd.xmin, pyppd.ymax-pyppd.ymin, linewidth=2, edgecolor='r', facecolor='none')
        ax_raw.set_xlim(0, pyppd.rawdata.shape[1])
        ax_raw.set_ylim(0, pyppd.rawdata.shape[0])
        ax_raw.add_patch(rect1)
        self.layoutWidget.display.raw.draw()

        #plot amplitude spectrum
        ax_spectrum = self.layoutWidget.display.spectrum.axes
        XYf2d_shifted = pyppd.XYf2d_shifted
        xfreq = np.fft.fftshift(np.fft.fftfreq(XYf2d_shifted.shape[1]))
        yfreq = np.fft.fftshift(np.fft.fftfreq(XYf2d_shifted.shape[0]))
        im3 = ax_spectrum.pcolormesh(xfreq, yfreq, XYf2d_shifted, vmax=self.spectrum_vmax)
        ax_spectrum.set_xlim(-0.2,0.2)
        ax_spectrum.set_ylim(-0.2,0.2)
        rect3 = patches.Rectangle((pyppd.fx-pyppd.xband,-np.abs(pyppd.fy)-pyppd.yband), 2*pyppd.xband, 2*(pyppd.yband+np.abs(pyppd.fy)), linewidth=2, edgecolor='r', facecolor='none')
        ax_spectrum.add_patch(rect3)
        self.layoutWidget.display.spectrum.draw()

        pyppd.filt_move()
        #find the center of phase spectrum
        try :
            pyppd.find_symmetry_axis()
        except RuntimeError :       #currently not possible because find_symmetry_axis always give a center in [ymin, ymax]
            self.statusBar().showMessage('Failed to find the symmetry axis.')
            return

        #plot phase spectrum
        ax_phase = self.layoutWidget.display.phase.axes
        ax_phase.set_title('Phase spectrum')
        im2 = ax_phase.pcolormesh(pyppd.phase)
        divider2 = make_axes_locatable(ax_phase)
        cax_phase = divider2.append_axes("right", size="5%", pad=0.05)
        plt.colorbar(im2, cax_phase)
        ax_phase.hlines(pyppd.ycenter, 0, pyppd.phase.shape[1], linewidth=3, colors='black')
        self.layoutWidget.display.phase.draw()

        #perform abel transform
        try :
            pyppd.abel()
        except ValueError :     #given invalid symmetry axis
            self.statusBar().showMessage('Failed to perform Abel transform.')
            return

        #plot desnity (relative refractivity)
        ax_density = self.layoutWidget.display.density.axes
        ax_density.set_title('Relative Refractivity')
        im4 = ax_density.pcolormesh(pyppd.AIM, vmax=0.1, vmin=0)
        ax_density.set_xlim(0, pyppd.AIM.shape[1])
        ax_density.set_ylim(0, pyppd.AIM.shape[0])
        divider4 = make_axes_locatable(ax_density)
        cax_density = divider4.append_axes("right", size="5%", pad=0.05)
        plt.colorbar(im4, cax_density)
        self.layoutWidget.display.density.draw()

    def about(self):
        QMessageBox.about(self, "About", """
        Python Plasma Density Diagnostics (PPDD) is an automatic and self-learning program to perform spectrum analysis and backward abel transform to get the plasma density distribution out of interferometric images, based on numpy, scipy and qt.

        Author: XU Zhiyi
        sgsdxzy@gmail.com
        State Key Laboratory of Nuclear Physics and Technology, and Key Laboratory of HEDP of the Ministry
of Education, CAPT, Peking University, Beijing 100871, China
        """)


class PPDDLayout(QWidget):
    def __init__(self):
        super().__init__()
        self.initUI()

    def initUI(self):
        self.conf = PPDDConf()
        self.display = PPDDDisplay()

        self.grid = QGridLayout()
        self.setLayout(self.grid)
        self.grid.setColumnStretch(0, 1)
        self.grid.setColumnStretch(1, 9)

        self.grid.addWidget(self.conf, 0, 0)
        self.grid.addWidget(self.display, 0, 1)

class PPDDConf(QWidget):
    def __init__(self):
        super().__init__()
        self.initUI()
        
    def initUI(self):
        self.grid = QGridLayout()
        self.grid.setColumnStretch(0, 3)
        self.grid.setColumnStretch(1, 7)
        self.setLayout(self.grid)

        self.grid.addWidget(QLabel(""), 0, 0)
        self.fx = QLineEdit()
        self.fx.setValidator(QDoubleValidator(0, 0.5, 6))
        self.grid.addWidget(self.fx, 0, 1)

        self.grid.addWidget(QLabel("fx"), 0, 0)
        self.fx = QLineEdit()
        self.fx.setValidator(QDoubleValidator(0, 0.5, 6))
        self.grid.addWidget(self.fx, 0, 1)

        self.grid.addWidget(QLabel("fx"), 0, 0)
        self.fx = QLineEdit()
        self.fx.setValidator(QDoubleValidator(0, 0.5, 6))
        self.grid.addWidget(self.fx, 0, 1)

        self.grid.addWidget(QLabel("fy"), 1, 0)
        self.fy = QLineEdit()
        self.fy.setValidator(QDoubleValidator(0, 0.5, 6))
        self.grid.addWidget(self.fy, 1, 1)



class PPDDDisplay(QWidget):
    def __init__(self):
        super().__init__()
        self.initUI()
        
    def initUI(self):
        grid = QGridLayout()
        grid.setColumnStretch(0, 45)
        grid.setColumnStretch(1, 55)
        self.setLayout(grid)

        self.raw = MplCanvas()
        self.phase = MplCanvas()
        self.spectrum = MplCanvas()
        self.density = MplCanvas()

        grid.addWidget(self.raw, 0, 0)
        grid.addWidget(self.phase, 0, 1)
        grid.addWidget(self.spectrum, 1, 0)
        grid.addWidget(self.density, 1, 1)
        

class MplCanvas(FigureCanvas):
    """Ultimately, this is a QWidget (as well as a FigureCanvasAgg, etc.)."""

    def __init__(self, parent=None, width=5, height=4, dpi=100):
        fig = Figure(figsize=(width, height), dpi=dpi)
        self.axes = fig.add_subplot(111)
        # We want the axes cleared every time plot() is called
        self.axes.hold(False)

        FigureCanvas.__init__(self, fig)
        self.setParent(parent)

        FigureCanvas.setSizePolicy(self, QSizePolicy.Expanding, QSizePolicy.Expanding)
        FigureCanvas.updateGeometry(self)


if __name__ == '__main__':
    app = QApplication(sys.argv)
    ex = PPDDWindow()
    sys.exit(app.exec_())
