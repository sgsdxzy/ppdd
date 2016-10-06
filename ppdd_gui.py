import sys
import os
from PyQt5.QtWidgets import *
from PyQt5.QtGui import *

import matplotlib
# Make sure that we are using QT5
matplotlib.use('Qt5Agg')
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from matplotlib import colorbar
import numpy as np

import ppdd

class PPDDWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.filename = ''
        self.ppdd = ppdd.PPDD()
        self.batch = batchRunWindow(self)
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
        left.grid.setRowStretch(0, 3)
        left.grid.setRowStretch(1, 4)
        left.grid.setRowStretch(2, 3)
        left.up = QWidget()
        left.down = QWidget()
        left.grid.addWidget(left.up, 0, 0)
        left.grid.addWidget(left.down, 2, 0)

        #the left up part
        up = left.up
        up.grid = QGridLayout()
        up.setLayout(up.grid)
        up.grid.setColumnStretch(0, 5)
        up.grid.setColumnStretch(1, 5)

        up.grid.addWidget(QLabel("X range"), 0, 0, 1, 2)
        up.xmin = QLineEdit()
        up.xmin.setValidator(QIntValidator())
        up.xmin.textEdited.connect(self.ppdd.reset)
        up.grid.addWidget(up.xmin, 1, 0)
        up.xmax = QLineEdit()
        up.xmax.setValidator(QIntValidator())
        up.xmax.textEdited.connect(self.ppdd.reset)
        up.grid.addWidget(up.xmax, 1, 1)

        up.grid.addWidget(QLabel("Y range"), 2, 0, 1, 2)
        up.ymin = QLineEdit()
        up.ymin.setValidator(QIntValidator())
        up.ymin.textEdited.connect(self.ppdd.reset)
        up.grid.addWidget(up.ymin, 3, 0)
        up.ymax = QLineEdit()
        up.ymax.setValidator(QIntValidator())
        up.ymax.textEdited.connect(self.ppdd.reset)
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
        for i in range(7):
            down.grid.setRowStretch(i, 1)

        down.grid.addWidget(QLabel("fx"), 0, 0)
        down.fx = QLineEdit()
        down.fx.setValidator(QDoubleValidator(0, 0.5, 6))
        down.fx.textEdited.connect(self.ppdd.reset)
        down.grid.addWidget(down.fx, 0, 1)

        down.grid.addWidget(QLabel("fy"), 1, 0)
        down.fy = QLineEdit()
        down.fy.setValidator(QDoubleValidator(0, 0.5, 6))
        down.fy.textEdited.connect(self.ppdd.reset)
        down.grid.addWidget(down.fy, 1, 1)

        down.grid.addWidget(QLabel("xband"), 2, 0)
        down.xband = QLineEdit()
        down.xband.setValidator(QDoubleValidator(0, 0.5, 6))
        down.grid.addWidget(down.xband, 2, 1)

        down.grid.addWidget(QLabel("yband"), 3, 0)
        down.yband = QLineEdit()
        down.yband.setValidator(QDoubleValidator(0, 0.5, 6))
        down.grid.addWidget(down.yband, 3, 1)

        down.learning = QCheckBox("Learning")
        down.learning.stateChanged.connect(self.ppdd.reset)
        down.grid.addWidget(down.learning, 4, 0, 1, 2)


        down.grid.addWidget(QLabel("method"), 5, 0, 1, 2)
        down.method = QComboBox()
        down.method.addItem("Hansen-Law")
        down.method.addItem("Onion-bordas")
        down.method.addItem("Basex")
        down.grid.addWidget(down.method, 6, 0, 1, 2)

        #the right column is used to display 4 plots in (2, 2) style
        right = layout.right
        right.grid = QGridLayout()
        right.setLayout(right.grid)
        right.grid.setColumnStretch(0, 45)
        right.grid.setColumnStretch(1, 55)
        right.grid.setRowStretch(0, 5)
        right.grid.setRowStretch(1, 5)

        right.raw = MplCanvas()
        right.phase = MplCanvas()
        right.phase.cax, _ = colorbar.make_axes(right.phase.axes, fraction = 0.05,  pad = 0.01, aspect = 10)
        right.phase.cax.hold(False)
        right.phase.cax.tick_params(axis='both', which='both', bottom='off', labelbottom='off')
        right.phase.detached.cax, _ = colorbar.make_axes(right.phase.detached.axes, fraction = 0.05,  pad = 0.01, aspect = 10)
        right.phase.detached.cax.hold(False)
        right.phase.detached.cax.tick_params(axis='both', which='both', bottom='off', labelbottom='off')
        right.spectrum = MplCanvas()
        right.density = MplCanvas()
        right.density.cax, _ = colorbar.make_axes(right.density.axes, fraction = 0.05,  pad = 0.01, aspect = 10)
        right.density.cax.hold(False)
        right.density.cax.tick_params(axis='both', which='both', bottom='off', labelbottom='off')
        right.density.detached.cax, _ = colorbar.make_axes(right.density.detached.axes, fraction = 0.05,  pad = 0.01, aspect = 10)
        right.density.detached.cax.hold(False)
        right.density.detached.cax.tick_params(axis='both', which='both', bottom='off', labelbottom='off')

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
        layout.left.down.learning.setChecked(pypdd.learning)

        openFile = QAction(QIcon.fromTheme('document-open'), 'Open', self)
        openFile.setShortcut('Ctrl+O')
        openFile.setStatusTip('Open file')
        openFile.triggered.connect(self.loadFileDialog)

        run = QAction(QIcon.fromTheme('media-playback-start'), 'Run', self)
        run.setShortcut('Ctrl+R')
        run.setStatusTip('Run ppdd')
        run.triggered.connect(self.runPPDD)

        saveFile = QAction(QIcon.fromTheme('document-save'), 'Save', self)
        saveFile.setShortcut('Ctrl+S')
        saveFile.setStatusTip('Save file')
        saveFile.triggered.connect(self.saveFile)

        saveFileAs = QAction(QIcon.fromTheme('document-save-as'), 'Save as', self)
        saveFileAs.setShortcut('Ctrl+D')
        saveFileAs.setStatusTip('Save file As')
        saveFileAs.triggered.connect(self.saveFileAs)

        batchRun = QAction(QIcon.fromTheme('media-seek-forward'), 'Batch run', self)
        batchRun.setShortcut('Ctrl+B')
        batchRun.setStatusTip('Batch run')
        batchRun.triggered.connect(self.batchRun)

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
        fileMenu.addAction(saveFile)
        fileMenu.addAction(saveFileAs)
        fileMenu.addAction(batchRun)
        fileMenu.addAction(exitAction)
        helpMenu = menubar.addMenu('&Help')
        helpMenu.addAction(aboutAction)

        toolbar = self.addToolBar('Open')
        toolbar.addAction(openFile)
        toolbar = self.addToolBar('Run')
        toolbar.addAction(run)
        toolbar = self.addToolBar('Save')
        toolbar.addAction(saveFile)
        toolbar = self.addToolBar('Save as')
        toolbar.addAction(saveFileAs)
        toolbar = self.addToolBar('Batch run')
        toolbar.addAction(batchRun)
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
        filenames = QFileDialog.getOpenFileName(self, 'Open file', '', 'Data file (*.txt);;Any file (*)', None, QFileDialog.DontUseNativeDialog)
        if filenames[0]:
            self.update_conf()
            pypdd = self.ppdd
            try :
                pypdd.readfile(filenames[0])
            except :
                error = QErrorMessage(self)
                error.showMessage('Error reading file: {0}'.format(filenames[0]))
                self.statusBar().showMessage('Error reading file: {0}'.format(filenames[0]))
                error.exec_() 
                return

            self.filename = filenames[0]
            self.statusBar().showMessage('Scucessfully read file: {0}'.format(filenames[0]))

            #plot raw data and selected region
            figure = self.layout.right.raw
            pypdd.plot_raw(figure.axes, region = (pypdd.xmin, pypdd.xmax, pypdd.ymin, pypdd.ymax))
            pypdd.plot_raw(figure.detached.axes, region = (pypdd.xmin, pypdd.xmax, pypdd.ymin, pypdd.ymax))
            figure.draw()
            figure.detached.draw()

    def saveFile(self):
        if not self.filename :      #check if a file is loaded
            self.statusBar().showMessage('No file is loaded, please load data file first!')
            return
        
        #make sure the result is what user wants
        self.runPPDD()

        outputpath = os.path.join(os.path.dirname(os.path.realpath(sys.argv[0])), 'output')
        os.makedirs(outputpath, exist_ok=True)
        outputfile = os.path.join(outputpath, os.path.basename(self.filename).rsplit('.', 1)[0]+'.txt')
        np.savetxt(outputfile, self.ppdd.AIM, fmt = '%1.6f', newline = os.linesep)
        self.statusBar().showMessage('Successfully saved file: {0}'.format(outputfile))

    def saveFileAs(self):
        if not self.filename :      #check if a file is loaded
            self.statusBar().showMessage('No file is loaded, please load data file first!')
            return

        outputpath = os.path.join(os.path.dirname(os.path.realpath(sys.argv[0])), 'output')
        os.makedirs(outputpath, exist_ok=True)
        outputfile = os.path.join(outputpath, os.path.basename(self.filename).rsplit('.', 1)[0]+'.txt')
        filenames = QFileDialog.getSaveFileName(self, 'Save File', outputfile, 'Data file (*.txt);;Any file (*)', None, QFileDialog.DontUseNativeDialog) 
        if filenames[0]:
            #only if user actually selected a file
            #make sure the result is what user wants
            self.runPPDD()

            np.savetxt(filenames[0], self.ppdd.AIM, fmt = '%1.6f', newline = os.linesep)
            self.statusBar().showMessage('Successfully saved file: {0}'.format(filenames[0]))

    def batchRun(self):
        #toggle show/hide
        if self.batch.isVisible() :
            self.batch.hide()
        else :
            self.batch.show()

    def update_conf(self):
        #TODO validate inputs

        pypdd = self.ppdd
        pypdd.xmin = int(self.layout.left.up.xmin.text())
        pypdd.xmax = int(self.layout.left.up.xmax.text())
        pypdd.ymin = int(self.layout.left.up.ymin.text())
        pypdd.ymax = int(self.layout.left.up.ymax.text())

        pypdd.symin = int(self.layout.left.up.symin.text())
        pypdd.symax = int(self.layout.left.up.symax.text())
        pypdd.xband = float(self.layout.left.down.xband.text())
        pypdd.yband = float(self.layout.left.down.yband.text())

        pypdd.learning = self.layout.left.down.learning.isChecked()
        pypdd.method = self.method_dict[str(self.layout.left.down.method.currentText())]


    def runPPDD(self):
        try :
            if not self.filename :      #check if a file is loaded
                self.statusBar().showMessage('No file is loaded, please load data file first!')
                return

            self.update_conf()
            pypdd = self.ppdd
            #plot raw data and selected region
            figure = self.layout.right.raw
            pypdd.plot_raw(figure.axes, region = (pypdd.xmin, pypdd.xmax, pypdd.ymin, pypdd.ymax))
            pypdd.plot_raw(figure.detached.axes, region = (pypdd.xmin, pypdd.xmax, pypdd.ymin, pypdd.ymax))
            figure.draw()
            figure.detached.draw()

            #get fx and fy from input 
            pypdd.guess.fx = float(self.layout.left.down.fx.text())
            pypdd.guess.fy = float(self.layout.left.down.fy.text())
            try :
                pypdd.find_peaks()
            except RuntimeError :
                self.statusBar().showMessage('Failed to find the secondary peak. Please enable learning, input fx and fy and run again.')
                #plot amplitude spectrum without passbands
                figure = self.layout.right.spectrum
                pypdd.plot_amplitude(figure.axes)
                pypdd.plot_amplitude(figure.detached.axes)
                figure.draw()
                figure.detached.draw()
                return

            #display fx and fy
            self.layout.left.down.fx.setText('{0:.4f}'.format(pypdd.fx))
            self.layout.left.down.fy.setText('{0:.4f}'.format(pypdd.fy))
                    
            #plot amplitude spectrum
            figure = self.layout.right.spectrum
            pypdd.plot_amplitude(figure.axes, bands=(pypdd.xband, pypdd.yband))
            pypdd.plot_amplitude(figure.detached.axes, bands=(pypdd.xband, pypdd.yband))
            figure.draw()
            figure.detached.draw()

            #get phase spectrum
            try :
                pypdd.filt_move()
            except ValueError as e:
                self.statusBar().showMessage(str(e))
                return
                

            #find the center of phase spectrum
            try :
                pypdd.find_symmetry_axis()
            except RuntimeError :       #currently not possible because find_symmetry_axis always give a center in [ymin, ymax]
                self.statusBar().showMessage('Failed to find the symmetry axis.')
                #plot phase spectrum without center line
                figure = self.layout.right.phase
                pypdd.plot_phase(figure.axes, figure.cax, limits = (pypdd.symin, pypdd.symax))
                pypdd.plot_phase(figure.detached.axes, figure.detached.cax, limits = (pypdd.symin, pypdd.symax))
                figure.draw()
                figure.detached.draw()
                return

            #plot phase spectrum
            figure = self.layout.right.phase
            pypdd.plot_phase(figure.axes, figure.cax, limits = (pypdd.symin, pypdd.symax), symmetry = pypdd.ycenter)
            pypdd.plot_phase(figure.detached.axes, figure.detached.cax, limits = (pypdd.symin, pypdd.symax), symmetry = pypdd.ycenter)
            figure.draw()
            figure.detached.draw()

            #perform abel transform
            try :
                pypdd.abel()
            except ValueError :     #given invalid symmetry axis
                self.statusBar().showMessage('Failed to perform Abel transform.')
                return

            #plot desnity (relative refractivity)
            figure = self.layout.right.density
            pypdd.plot_density(figure.axes, figure.cax)
            pypdd.plot_density(figure.detached.axes, figure.detached.cax)
            figure.draw()
            figure.detached.draw()

        except Exception as e:
            self.statusBar().showMessage(str(e))
            raise

    def about(self):
        QMessageBox.about(self, "About", """
        Python Plasma Density Diagnostics (PPDD) is an automatic and self-learning program to perform spectrum analysis and backward abel transform to get the plasma density distribution out of interferometric images, based on numpy, scipy and qt.

        Author: XU Zhiyi
        sgsdxzy@gmail.com
        State Key Laboratory of Nuclear Physics and Technology, and Key Laboratory of HEDP of the Ministry
of Education, CAPT, Peking University, Beijing 100871, China
        """)


class batchRunWindow(QWidget):
    """A new windows used to dsiplay batch run status"""

    def __init__(self, mainWindow):
        super().__init__()
        self.mainWindow = mainWindow
        self.initUI()
        
    def initUI(self): 
        self.setGeometry(300, 300, 800, 600)
        self.setWindowTitle('Batch run') 

        self.grid = QGridLayout()
        self.setLayout(self.grid)

        self.toolbar = QToolBar("Batch run", self)

        openFiles = QAction(QIcon.fromTheme('document-open'), 'Open files', self)
        openFiles.setShortcut('Ctrl+O')
        openFiles.setStatusTip('Open files')
        openFiles.triggered.connect(self.loadFilesDialog)

        run = QAction(QIcon.fromTheme('media-seek-forward'), 'Batch run', self)
        run.setShortcut('Ctrl+R')
        run.setStatusTip('Batch run')
        run.triggered.connect(self.run)

        self.toolbar.addAction(openFiles)
        self.toolbar.addAction(run)
        self.grid.addWidget(self.toolbar, 0, 0)

        self.progress = QProgressBar(self)
        self.grid.addWidget(self.progress, 1, 0)

        #create a canvas for showing a snapshot of the result of current file
        self.fig = Figure()
        self.axes = self.fig.add_subplot(111)
        self.axes.hold(False)
        self.cax, _ = colorbar.make_axes(self.axes, fraction = 0.05,  pad = 0.01, aspect = 10)
        self.cax.hold(False)
        self.cax.tick_params(axis='both', which='both', bottom='off', labelbottom='off')
        self.canvas = FigureCanvas(self.fig)
        self.canvas.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        self.canvas.updateGeometry()
        self.grid.addWidget(self.canvas, 2, 0)

        self.status = QStatusBar(self)
        self.grid.addWidget(self.status, 3, 0)

        self.grid.setRowStretch(0, 0)
        self.grid.setRowStretch(1, 0)
        self.grid.setRowStretch(3, 0)

        self.status.showMessage('Ready.')

    def loadFilesDialog(self):
        filenames = QFileDialog.getOpenFileNames(self, 'Open files', '', 'Data file (*.txt);;Any file (*)', None, QFileDialog.DontUseNativeDialog)
        if filenames[0]:
            self.status.showMessage('Selected {0} file(s): {1}...'.format(len(filenames[0]), filenames[0][0]))
            self.filenames = filenames[0]

    def run(self):
        if not self.filenames :
            self.status.showMessage('No file is selected, please select data files first!')
            return

        self.mainWindow.update_conf()
        outputpath = os.path.join(os.path.dirname(os.path.realpath(sys.argv[0])), 'output')
        os.makedirs(outputpath, exist_ok=True)
        total = len(self.filenames)
        success = 0
        fail = 0
        self.progress.setMinimum(0)
        self.progress.setMaximum(total)
        pypdd = self.mainWindow.ppdd

        #batch run
        failed_files = []
        for f in self.filenames:
            self.progress.setValue(success + fail + 1)
            try :
                pypdd.readfile(f)
                pypdd.find_peaks()
                pypdd.filt_move()
                pypdd.find_symmetry_axis()
                pypdd.abel()
                pypdd.plot_density(self.axes, self.cax)

                outputfile = os.path.join(outputpath, os.path.basename(f).rsplit('.', 1)[0]+'.txt')
                np.savetxt(outputfile, pypdd.AIM, fmt = '%1.6f', newline = os.linesep)
                self.status.showMessage('Successfully saved file: {0}'.format(outputfile))
                success += 1
            except :
                failed_files.append(f)
                fail += 1
                continue

        self.status.showMessage('Total number of files: {0}, success number of files: {1}.'.format(total, success))

        if failed_files : 
            error = QErrorMessage(self)
            failed_str = os.linesep.join(failed_files)
            error.showMessage('Failed to analyze the following {0} file(s): {1}{2}'.format(fail, os.linesep, failed_str))
            error.exec_() 


class MplCanvas(FigureCanvas):
    """Ultimately, this is a QWidget (as well as a FigureCanvasAgg, etc.)."""

    def __init__(self, parent=None, width=5, height=4, dpi=100):
        self.fig = Figure(figsize=(width, height), dpi=dpi)
        self.axes = self.fig.add_subplot(111)
        # We want the axes cleared every time plot() is called
        self.axes.hold(False)

        FigureCanvas.__init__(self, self.fig)
        self.setParent(parent)

        FigureCanvas.setSizePolicy(self, QSizePolicy.Expanding, QSizePolicy.Expanding)
        FigureCanvas.updateGeometry(self)

        self.detached = DetachedCanvas()

    def mousePressEvent(self, event):
        #toggle show/hide
        if self.detached.isVisible() :
            self.detached.hide()
        else :
            self.detached.show()

class DetachedCanvas(FigureCanvas):
    """A detached canvas used to show enlarged plots, and save individual plots."""

    def __init__(self, parent=None):
        self.fig = Figure()
        self.axes = self.fig.add_subplot(111)
        # We want the axes cleared every time plot() is called
        self.axes.hold(False)

        FigureCanvas.__init__(self, self.fig)
        self.setParent(parent)

        FigureCanvas.setSizePolicy(self, QSizePolicy.Expanding, QSizePolicy.Expanding)
        FigureCanvas.updateGeometry(self)




if __name__ == '__main__':
    app = QApplication(sys.argv)
    ex = PPDDWindow()
    sys.exit(app.exec_())
