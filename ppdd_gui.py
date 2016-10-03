import sys
from PyQt5.QtWidgets import (QMainWindow, QWidget, QTextEdit, QAction, QApplication,
                            QMessageBox, QGridLayout, QPushButton, QFileDialog, QErrorMessage,
                            QSizePolicy, QLabel, QLineEdit, QComboBox)
from PyQt5.QtGui import QIcon, QDoubleValidator, QIntValidator
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
        self.filename = ''
        self.ppdd = ppdd.PPDD()
        self.spectrum_vmax = 1e6        #used to display the amplitude spectrum
        self.method_dict = {
            "Hansen-Law": "hansenlaw",
            "Onion-bordas": "onion_bordas",
            "Basex": "basex"
        }
        self.initUI()
        
    def initUI(self): 
        #using a grid layout Qwidget as centeralwidget
        self.layout = QWidget()
        layout = self.layout
        self.setCentralWidget(layout)

        #divide layout to left and right two columns
        layout.grid = QGridLayout()
        layout.setLayout(layout.grid)
        layout.grid.setColumnStretch(0, 1)
        layout.grid.setColumnStretch(1, 9)
        layout.left = QWidget()
        layout.right = QWidget()
        layout.grid.addWidget(layout.left, 0, 0)
        layout.grid.addWidget(layout.right, 0, 1)

        #the left column is splited to up and down two rows
        left = layout.left
        left.grid = QGridLayout()
        left.setLayout(left.grid)
        #left.grid.setRowStretch(0, 3)
        #left.grid.setRowStretch(1, 7)
        left.up = QWidget()
        left.down = QWidget()
        left.grid.addWidget(left.up, 0, 0)
        left.grid.addWidget(left.down, 1, 0)

        #the left up part
        up = left.up
        up.grid = QGridLayout()
        up.setLayout(up.grid)
        up.grid.setColumnStretch(0, 5)
        up.grid.setColumnStretch(1, 5)

        up.grid.addWidget(QLabel("X range"), 0, 0, 1, 2)
        up.xmin = QLineEdit()
        up.xmin.setValidator(QIntValidator())
        up.grid.addWidget(up.xmin, 1, 0)
        up.xmax = QLineEdit()
        up.xmax.setValidator(QIntValidator())
        up.grid.addWidget(up.xmax, 1, 1)

        up.grid.addWidget(QLabel("Y range"), 2, 0, 1, 2)
        up.ymin = QLineEdit()
        up.ymin.setValidator(QIntValidator())
        up.grid.addWidget(up.ymin, 3, 0)
        up.ymax = QLineEdit()
        up.ymax.setValidator(QIntValidator())
        up.grid.addWidget(up.ymax, 3, 1)

        up.grid.addWidget(QLabel("Symmetry range"), 4, 0, 1, 2)
        up.symin = QLineEdit()
        up.symin.setValidator(QIntValidator())
        up.grid.addWidget(up.symin, 5, 0)
        up.symax = QLineEdit()
        up.symax.setValidator(QIntValidator())
        up.grid.addWidget(up.symax, 5, 1)

        #the left down part
        down = left.down
        down.grid = QGridLayout()
        down.setLayout(down.grid)
        down.grid.setColumnStretch(0, 3)
        down.grid.setColumnStretch(1, 7)

        down.grid.addWidget(QLabel("fx"), 0, 0)
        down.fx = QLineEdit()
        down.fx.setValidator(QDoubleValidator(0, 0.5, 6))
        down.grid.addWidget(down.fx, 0, 1)

        down.grid.addWidget(QLabel("fy"), 1, 0)
        down.fy = QLineEdit()
        down.fy.setValidator(QDoubleValidator(0, 0.5, 6))
        down.grid.addWidget(down.fy, 1, 1)

        down.grid.addWidget(QLabel("xband"), 2, 0)
        down.xband = QLineEdit()
        down.xband.setValidator(QDoubleValidator(0, 0.5, 6))
        down.grid.addWidget(down.xband, 2, 1)

        down.grid.addWidget(QLabel("yband"), 3, 0)
        down.yband = QLineEdit()
        down.yband.setValidator(QDoubleValidator(0, 0.5, 6))
        down.grid.addWidget(down.yband, 3, 1)

        down.grid.addWidget(QLabel("method"), 4, 0)
        down.method = QComboBox()
        down.method.addItem("Hansen-Law")
        down.method.addItem("Onion-bordas")
        down.method.addItem("Basex")
        down.grid.addWidget(down.method, 4, 1)

        #the right column is used to display 4 plots in (2, 2) style
        right = layout.right
        right.grid = QGridLayout()
        right.setLayout(right.grid)
        right.grid.setColumnStretch(0, 45)
        right.grid.setColumnStretch(1, 55)

        right.raw = MplCanvas()
        right.phase = MplCanvas()
        right.phase.cax = make_axes_locatable(right.phase.axes).append_axes("right", size="5%", pad=0.05)       #create a cax for plotting colorbar
        right.phase.cax.hold(False)
        right.phase.cax.tick_params(axis='both', which='both', bottom='off', labelbottom='off')
        right.spectrum = MplCanvas()
        right.density = MplCanvas()
        right.density.cax = make_axes_locatable(right.density.axes).append_axes("right", size="5%", pad=0.05)   #create a cax for plotting colorbar
        right.density.cax.hold(False)
        right.density.cax.tick_params(axis='both', which='both', bottom='off', labelbottom='off')

        right.grid.addWidget(right.raw, 0, 0)
        right.grid.addWidget(right.phase, 0, 1)
        right.grid.addWidget(right.spectrum, 1, 0)
        right.grid.addWidget(right.density, 1, 1)

        #write inital values
        pypdd = self.ppdd
        layout.left.up.xmin.setText(str(pypdd.xmin))
        layout.left.up.xmax.setText(str(pypdd.xmax))
        layout.left.up.ymin.setText(str(pypdd.ymin))
        layout.left.up.ymax.setText(str(pypdd.ymax))
        layout.left.up.symin.setText(str(pypdd.symin))
        layout.left.up.symax.setText(str(pypdd.symax))
        layout.left.down.fx.setText('{0:.4f}'.format(pypdd.guess.fx))
        layout.left.down.fy.setText('{0:.4f}'.format(pypdd.guess.fy))
        layout.left.down.xband.setText('{0:.3f}'.format(pypdd.xband))
        layout.left.down.yband.setText('{0:.3f}'.format(pypdd.yband))


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

    def update_conf(self):
        #TODO validate inputs

        pypdd = self.ppdd
        newxmin = int(self.layout.left.up.xmin.text())
        newxmax = int(self.layout.left.up.xmax.text())
        newymin = int(self.layout.left.up.ymin.text())
        newymax = int(self.layout.left.up.ymax.text())

        if not ((newxmin == pypdd.xmin) and (newymin == pypdd.ymin) and (newxmax == pypdd.xmax) and (newymax == pypdd.ymax)) :
            #input has changed
            pypdd.xmin = newxmin
            pypdd.xmax = newxmax
            pypdd.ymin = newymin
            pypdd.ymax = newymax
            pypdd.peak_fitted = False #bcause input has changed

        pypdd.symin = int(self.layout.left.up.symin.text())
        pypdd.symax = int(self.layout.left.up.symax.text())
        pypdd.xband = float(self.layout.left.down.xband.text())
        pypdd.yband = float(self.layout.left.down.yband.text())

        pypdd.method = self.method_dict[str(self.layout.left.down.method.currentText())]


    def runPPDD(self):
        if not self.filename :      #check if a file is loaded
            self.statusBar().showMessage('No file is loaded, please load data file first!')
            return

        self.update_conf()
        pypdd = self.ppdd

        #plot raw data and selected region
        ax_raw = self.layout.right.raw.axes
        ax_raw.set_title('Raw data')
        im1 = ax_raw.pcolormesh(pypdd.rawdata)
        rect1 = patches.Rectangle((pypdd.xmin, pypdd.ymin), pypdd.xmax-pypdd.xmin, pypdd.ymax-pypdd.ymin, linewidth=2, edgecolor='r', facecolor='none')
        ax_raw.set_xlim(0, pypdd.rawdata.shape[1])
        ax_raw.set_ylim(0, pypdd.rawdata.shape[0])
        ax_raw.add_patch(rect1)
        self.layout.right.raw.draw()

        if not pypdd.peak_fitted :  #if already fitted peaks, skip to speed up
            try :
                pypdd.guess.fx = float(self.layout.left.down.fx.text())
                pypdd.guess.fy = float(self.layout.left.down.fy.text())
                pypdd.find_peaks()
            except RuntimeError :
                self.statusBar().showMessage('Failed to find the secondary peak. Please input fx and fy and run again.')
                #plot amplitude spectrum without passbands
                ax_spectrum = self.layout.right.spectrum.axes
                XYf2d = np.fft.fftn(pypdd.xy2d)
                XYf2d_shifted = np.abs(np.fft.fftshift(XYf2d))     
                xfreq = np.fft.fftshift(np.fft.fftfreq(XYf2d_shifted.shape[1]))
                yfreq = np.fft.fftshift(np.fft.fftfreq(XYf2d_shifted.shape[0]))
                im3 = ax_spectrum.pcolormesh(xfreq, yfreq, XYf2d_shifted, vmax=self.spectrum_vmax)
                ax_spectrum.set_xlim(-0.2,0.2)
                ax_spectrum.set_ylim(-0.2,0.2)
                self.layout.right.spectrum.draw()
                return

        self.layout.left.down.fx.setText('{0:.4f}'.format(pypdd.fx))
        self.layout.left.down.fy.setText('{0:.4f}'.format(pypdd.fy))
                
        #plot amplitude spectrum
        ax_spectrum = self.layout.right.spectrum.axes
        XYf2d_shifted = pypdd.XYf2d_shifted
        xfreq = np.fft.fftshift(np.fft.fftfreq(XYf2d_shifted.shape[1]))
        yfreq = np.fft.fftshift(np.fft.fftfreq(XYf2d_shifted.shape[0]))
        im3 = ax_spectrum.pcolormesh(xfreq, yfreq, XYf2d_shifted, vmax=self.spectrum_vmax)
        ax_spectrum.set_xlim(-0.2,0.2)
        ax_spectrum.set_ylim(-0.2,0.2)
        rect3 = patches.Rectangle((pypdd.fx-pypdd.xband,-np.abs(pypdd.fy)-pypdd.yband), 2*pypdd.xband, 2*(pypdd.yband+np.abs(pypdd.fy)), linewidth=2, edgecolor='r', facecolor='none')
        ax_spectrum.add_patch(rect3)
        self.layout.right.spectrum.draw()

        pypdd.filt_move()

        #find the center of phase spectrum
        try :
            pypdd.find_symmetry_axis()
        except RuntimeError :       #currently not possible because find_symmetry_axis always give a center in [ymin, ymax]
            self.statusBar().showMessage('Failed to find the symmetry axis.')
            #plot phase spectrum without center line
            ax_phase = self.layout.right.phase.axes
            ax_phase.set_title('Phase spectrum')
            im2 = ax_phase.pcolormesh(pypdd.phase)
            plt.colorbar(im2, self.layout.right.phase.cax)
            self.layout.right.phase.draw()
            return

        #plot phase spectrum
        ax_phase = self.layout.right.phase.axes
        ax_phase.set_title('Phase spectrum')
        im2 = ax_phase.pcolormesh(pypdd.phase)
        plt.colorbar(im2, self.layout.right.phase.cax)
        ax_phase.hlines(pypdd.ycenter, 0, pypdd.phase.shape[1], linewidth=3, colors='black')
        self.layout.right.phase.draw()

        #perform abel transform
        try :
            pypdd.abel()
        except ValueError :     #given invalid symmetry axis
            self.statusBar().showMessage('Failed to perform Abel transform.')
            return

        #plot desnity (relative refractivity)
        ax_density = self.layout.right.density.axes
        ax_density.set_title('Relative Refractivity')
        im4 = ax_density.pcolormesh(pypdd.AIM, vmax=0.1, vmin=0)
        ax_density.set_xlim(0, pypdd.AIM.shape[1])
        ax_density.set_ylim(0, pypdd.AIM.shape[0])
        plt.colorbar(im2, self.layout.right.density.cax)
        self.layout.right.density.draw()

    def about(self):
        QMessageBox.about(self, "About", """
        Python Plasma Density Diagnostics (PPDD) is an automatic and self-learning program to perform spectrum analysis and backward abel transform to get the plasma density distribution out of interferometric images, based on numpy, scipy and qt.

        Author: XU Zhiyi
        sgsdxzy@gmail.com
        State Key Laboratory of Nuclear Physics and Technology, and Key Laboratory of HEDP of the Ministry
of Education, CAPT, Peking University, Beijing 100871, China
        """)


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
