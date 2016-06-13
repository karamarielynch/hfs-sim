import sys
from PyQt4 import QtGui, QtCore
import pyqtgraph as pg
import numpy as np
import time
from scipy import signal 
import matplotlib.pyplot as plt
import ast
from hfs_creation import *
import os.path

x = np.linspace(-20*10**3, 60*10**3, 50*10**3, endpoint = False)

c = 299792458.0
		
class Example(QtGui.QMainWindow):
	
	def __init__(self):
		super(Example, self).__init__()

		self.initUI()
		
	def initUI(self):
			
		widget = QtGui.QWidget()
		self.setCentralWidget(widget)

		self.grid = QtGui.QGridLayout(widget)
		self.grid.setSpacing(10)
		
		self.setWindowTitle('Hyperfine spectrum simulator')

		self.IsotopeFilename = sys.argv[1]
		HFS_params = np.loadtxt(self.IsotopeFilename, delimiter='\t', dtype='S')
		self.IsotopeNames = []
		names = HFS_params.T[0] #.decode("utf-8")
		for i in range(len(names)):
			newnames = names[i].decode("utf-8")
			self.IsotopeNames.append(newnames)
		Isotope_value = HFS_params[0]
		self.I 			= float(Isotope_value[1])
		self.J_lower 	= float(Isotope_value[2])
		self.J_upper 	= float(Isotope_value[3])
		self.CF 		= float(Isotope_value[4])
		self.A_lower 	= float(Isotope_value[5])
		self.A_upper 	= float(Isotope_value[6])
		self.B_lower 	= float(Isotope_value[7])
		self.B_upper 	= float(Isotope_value[8])
		self.FWHM 		= float(Isotope_value[9])
		self.Mass 		= float(Isotope_value[10])
		self.RISLine	= float(Isotope_value[11])
		self.Mass_ref 	= float(Isotope_value[12])
		self.ARatio		= float(Isotope_value[13])
		self.Line		= float(Isotope_value[14])
		self.Harmonic	= float(Isotope_value[15])

		Isotope	= QtGui.QLabel('Isotope')
		self.grid.addWidget(Isotope, 2, 0)
		IsotopeEdit = QtGui.QLabel()
		self.IsotopeList = QtGui.QComboBox(self)
		self.IsotopeList.addItems(self.IsotopeNames)
		self.grid.addWidget(self.IsotopeList, 2, 1)
		self.IsotopeList.activated[str].connect(self.updateIsotope)
		self.IsotopeList.activated[str].connect(self.updateTopPlot)
		self.IsotopeList.activated[str].connect(self.updateBottomPlot)
		self.IsotopeList.activated[str].connect(self.updateHFSPeaks)

		Spin	= QtGui.QLabel('Spin')
		self.grid.addWidget(Spin, 3, 0)
		self.SpinEdit = pg.SpinBox(value=self.I, dec=True, minStep=0.5)
		self.SpinEdit.setRange(0, 12)
		self.grid.addWidget(self.SpinEdit, 3, 1)
		self.SpinEdit.valueChanged.connect(self.updateSpin)
		self.SpinEdit.valueChanged.connect(self.updateTopPlot)
		self.SpinEdit.valueChanged.connect(self.updateBottomPlot)
		self.SpinEdit.valueChanged.connect(self.updateHFSPeaks)

		Jl	= QtGui.QLabel('Jl')
		self.grid.addWidget(Jl, 3, 2)
		self.JlEdit = pg.SpinBox(value=self.J_lower, dec=True, minStep=0.5)
		self.JlEdit.setRange(0, 2.5)
		self.grid.addWidget(self.JlEdit, 3, 3)
		self.JlEdit.valueChanged.connect(self.updateJl)
		self.JlEdit.valueChanged.connect(self.updateTopPlot)
		self.JlEdit.valueChanged.connect(self.updateBottomPlot)
		self.JlEdit.valueChanged.connect(self.updateHFSPeaks)

		Ju	= QtGui.QLabel('Ju')
		self.grid.addWidget(Ju, 3, 4)
		self.JuEdit = pg.SpinBox(value=self.J_upper, dec=True, minStep=0.5)
		self.JuEdit.setRange(0, 2.5)
		self.grid.addWidget(self.JuEdit, 3, 5)
		self.JuEdit.valueChanged.connect(self.updateJu)
		self.JuEdit.valueChanged.connect(self.updateTopPlot)
		self.JuEdit.valueChanged.connect(self.updateBottomPlot)
		self.JuEdit.valueChanged.connect(self.updateHFSPeaks)
		

		Al	= QtGui.QLabel('Al (MHz)')
		self.grid.addWidget(Al, 4, 0)
		self.AlEdit = pg.SpinBox(value=self.A_lower, dec=True, minStep=0.1)
		self.AlEdit.setRange(-20000, 20000)
		self.grid.addWidget(self.AlEdit, 4, 1)
		self.AlEdit.valueChanged.connect(self.updateAl)
		self.AlEdit.valueChanged.connect(self.updateTopPlot)
		self.AlEdit.valueChanged.connect(self.updateBottomPlot)
		self.AlEdit.valueChanged.connect(self.updateHFSPeaks)

		Au	= QtGui.QLabel('Au (MHz)')
		self.grid.addWidget(Au, 4, 2)
		self.AuEdit = pg.SpinBox(value=self.A_upper, dec=True, minStep=0.1)
		self.AuEdit.setRange(-20000, 20000)
		self.grid.addWidget(self.AuEdit, 4, 3)
		self.AuEdit.valueChanged.connect(self.updateAu)
		self.AuEdit.valueChanged.connect(self.updateTopPlot)
		self.AuEdit.valueChanged.connect(self.updateBottomPlot)
		self.AuEdit.valueChanged.connect(self.updateHFSPeaks)

		Bl	= QtGui.QLabel('Bl (MHz)')
		self.grid.addWidget(Bl, 4, 4)
		self.BlEdit = pg.SpinBox(value=self.B_lower, dec=True, minStep=0.1)
		self.BlEdit.setRange(-20000, 20000)
		self.grid.addWidget(self.BlEdit, 4, 5)
		self.BlEdit.valueChanged.connect(self.updateBl)
		self.BlEdit.valueChanged.connect(self.updateTopPlot)
		self.BlEdit.valueChanged.connect(self.updateBottomPlot)
		self.BlEdit.valueChanged.connect(self.updateHFSPeaks)

		Bu	= QtGui.QLabel('Bu (MHz)')
		self.grid.addWidget(Bu, 4, 6)
		self.BuEdit = pg.SpinBox(value=self.B_upper, dec=True, minStep=0.1)
		self.BuEdit.setRange(-20000, 20000)
		self.grid.addWidget(self.BuEdit, 4, 7)
		self.BuEdit.valueChanged.connect(self.updateBu)
		self.BuEdit.valueChanged.connect(self.updateTopPlot)
		self.BuEdit.valueChanged.connect(self.updateBottomPlot)
		self.BuEdit.valueChanged.connect(self.updateHFSPeaks)

		CF	= QtGui.QLabel('CF (MHz)')
		self.grid.addWidget(CF, 5, 0)
		self.CFEdit = pg.SpinBox(value=self.CF, dec=True, minStep=1)
		self.CFEdit.setRange(-50000, 50000)
		self.grid.addWidget(self.CFEdit, 5, 1)
		self.CFEdit.valueChanged.connect(self.updateCF)
		self.CFEdit.valueChanged.connect(self.updateTopPlot)
		self.CFEdit.valueChanged.connect(self.updateBottomPlot)
		self.CFEdit.valueChanged.connect(self.updateHFSPeaks)

		FWHM	= QtGui.QLabel('FWHM (MHz)')
		self.grid.addWidget(FWHM, 5, 2)
		self.FWHMEdit = pg.SpinBox(value=self.FWHM, dec=True, minStep=1)
		self.FWHMEdit.setRange(0, 3000)
		self.grid.addWidget(self.FWHMEdit, 5, 3)
		self.FWHMEdit.valueChanged.connect(self.updateFWHM)
		self.FWHMEdit.valueChanged.connect(self.updateTopPlot)
		self.FWHMEdit.valueChanged.connect(self.updateBottomPlot)
		self.FWHMEdit.valueChanged.connect(self.updateHFSPeaks)

		ARatio	= QtGui.QLabel('Au/Al Ratio')
		self.grid.addWidget(ARatio, 5, 4)
		self.ARatioEnable = QtGui.QCheckBox('', self)
		self.grid.addWidget(self.ARatioEnable, 5, 5)
		self.ARatioEnable.stateChanged.connect(self.enableARatio)
		self.ARatioEnable.stateChanged.connect(self.updateTopPlot)
		self.ARatioEnable.stateChanged.connect(self.updateBottomPlot)
		self.ARatioEnable.stateChanged.connect(self.updateHFSPeaks)


		ISCOOL	= QtGui.QLabel('ISCOOL (V)')
		self.grid.addWidget(ISCOOL, 8, 0)
		self.ISCOOLEdit = pg.SpinBox(value=30000., dec=True, minStep=0.01)
		self.ISCOOLEdit.setRange(0, 60000)
		self.grid.addWidget(self.ISCOOLEdit, 8, 1)
		self.ISCOOLEdit.valueChanged.connect(self.updateISCOOL)
		self.ISCOOLEdit.valueChanged.connect(self.updateBottomPlot)
		self.ISCOOLEdit.valueChanged.connect(self.updateHFSPeaks)


		centroid = self.RISLine/self.Harmonic
		alpha	= self.ISCOOLEdit.value()/(self.Mass_ref*931.494061*10**6)
		centroid_doppler = centroid/( 1 + alpha - sqrt(2*alpha + alpha*alpha))		
		offset = centroid_doppler

		minrange = round(centroid_doppler,1) - 2.0 # cm-1
		maxrange = round(centroid_doppler,1) + 2.0 # cm-1

		calculateButton = QtGui.QPushButton('Reset Range', self)
		self.grid.addWidget(calculateButton, 9, 0)
		calculateButton.clicked.connect(self.calculateRanges)
		calculateButton.clicked.connect(self.updateRange)
		calculateButton.clicked.connect(self.updateBottomPlot)
		#calculateButton.clicked.connect(self.updateHFSPeaks)
		calculateButton.setToolTip('Recalculate wavenumber range')

		MinRange	= QtGui.QLabel('From (cm-1)')
		self.grid.addWidget(MinRange, 8, 2)
		self.MinRangeEdit = QtGui.QDoubleSpinBox(decimals=3)
		self.MinRangeEdit.setSingleStep(0.001)
		self.MinRangeEdit.setRange(0, 100000)
		self.MinRangeEdit.setValue(minrange)
		self.grid.addWidget(self.MinRangeEdit, 8, 3)
		self.MinRangeEdit.valueChanged.connect(self.updateRange)
		self.MinRangeEdit.valueChanged.connect(self.updateBottomPlot)

		MaxRange	= QtGui.QLabel('To (cm-1)')
		self.grid.addWidget(MaxRange, 8, 4)
		self.MaxRangeEdit = QtGui.QDoubleSpinBox(decimals=3)
		self.MaxRangeEdit.setSingleStep(0.001)
		self.MaxRangeEdit.setRange(0, 100000)
		self.MaxRangeEdit.setValue(maxrange)
		self.grid.addWidget(self.MaxRangeEdit, 8, 5)
		self.MaxRangeEdit.valueChanged.connect(self.updateRange)
		self.MaxRangeEdit.valueChanged.connect(self.updateBottomPlot)

		Offset	= QtGui.QLabel('Offset (cm-1)')
		self.grid.addWidget(Offset, 8, 6)
		self.OffsetEdit = QtGui.QDoubleSpinBox(decimals=3)
		self.OffsetEdit.setSingleStep(0.001)
		self.OffsetEdit.setRange(0, 100000)
		self.OffsetEdit.setValue(offset)
		self.grid.addWidget(self.OffsetEdit, 8, 7)
		self.OffsetEdit.valueChanged.connect(self.updateOffset)
		self.OffsetEdit.valueChanged.connect(self.updateBottomPlot)

		Peaks 	= QtGui.QLabel('Peaks (cm-1)')
		Intensity 	= QtGui.QLabel('Intensities')
		self.grid.addWidget(Peaks, 9, 1)
		self.grid.addWidget(Intensity, 10, 1)
		self.PeaksEdits = []
		self.IntensityEdits = []
		for i in range(len(HF_function(self.I, self.J_lower, self.J_upper, self.CF, self.A_lower, self.A_upper, self.B_lower, self.B_upper)[0])):
			PeakValue_freq = HF_function(self.I, self.J_lower, self.J_upper, self.CF, self.A_lower, self.A_upper, self.B_lower, self.B_upper)[0][i]
			alpha	= self.ISCOOLEdit.value()/(self.Mass*931.494061*10**6)
			PeakValue_dopp = PeakValue_freq/( 1 + alpha - sqrt(2*alpha + alpha*alpha))
			PeakValue_wave = PeakValue_dopp*10**4/(self.Harmonic*c) + offset
			self.PeaksEdit = QtGui.QLabel(str(round(PeakValue_wave, 3)))
			self.grid.addWidget(self.PeaksEdit, 9, i+2)
			self.PeaksEdits.append(self.PeaksEdit)
			IntensityValue = HF_function(self.I, self.J_lower, self.J_upper, self.CF, self.A_lower, self.A_upper, self.B_lower, self.B_upper)[1][i]
			self.IntensityEdit = QtGui.QLabel(str(round(IntensityValue, 1)))
			self.grid.addWidget(self.IntensityEdit, 10, i+2)
			self.IntensityEdits.append(self.IntensityEdit)

		self.plotWidget = pg.PlotWidget()
		self.grid.addWidget(self.plotWidget, 0, 0, 2, 11)
		self.plotWidget.setLabel('bottom', "Frequency relative to reference isotope", units='Hz', color='#0092ff')
		self.plotWidget.setTitle(str(self.Line)+" nm transition", color='#0092ff')
		self.plotWidget.showAxis('left', False)
		self.plotWidget.plot(x*10**6, HFS(self.I, self.J_lower, self.J_upper, self.CF, self.A_lower, self.A_upper, self.B_lower, self.B_upper, self.FWHM, 10, 0, x), pen=tuple([0,146,255]))


		self.wavenumber_range = np.linspace(minrange, maxrange, 10**5)
		self.centroid = (self.RISLine*c)/10**4 # RISLine taken from NIST and corrected for reference isotope (units of cm-1)
		self.freq_range_lab = (self.Harmonic*self.wavenumber_range*c)/10**4
		self.freq_range = Doppler_correction(self.freq_range_lab, self.Mass, self.ISCOOLEdit.value())


		self.plotWidget2 = pg.PlotWidget()
		self.grid.addWidget(self.plotWidget2, 6, 0, 2, 11)
		self.plotWidget2.setLabel('bottom', "Wavenumber offset by "+str(self.OffsetEdit.value()), units='cm-1', color='#00ffff')
		self.plotWidget2.setTitle(str(self.Line)+" nm transition", color='#00ffff')
		self.plotWidget2.showAxis('left', False)
		self.plotWidget2.plot(self.wavenumber_range-self.OffsetEdit.value(), HFS(self.I, self.J_lower, self.J_upper, self.CF, self.A_lower, self.A_upper, self.B_lower, self.B_upper, self.FWHM, 10, 0, (self.freq_range-self.centroid) ), pen=tuple([0,255,255]))

		self.old_isotope = self.IsotopeList.currentText()
		
		self.show()
		


	def updateTopPlot(self):
		self.plotWidget.plot(clear=True)
		self.plotWidget.setTitle(str(self.Line)+" nm transition", color='#0092ff')
		try:
			self.plotWidget.plot(x*10**6, HFS(self.I, self.J_lower, self.J_upper, self.CF, self.A_lower, self.A_upper, self.B_lower, self.B_upper, self.FWHM, 10, 0, x ), pen=tuple([0,146,255]))
		except ValueError:
			self.showInfo()

	def showInfo(self):
		msg = QtGui.QMessageBox()
		msg.setIcon(QtGui.QMessageBox.Warning)
		msg.setWindowTitle('Error for Jl -> Ju')
		msg.setText("Please choose a valid transition")
		msg.exec_()

	def updateBottomPlot(self):
		self.centroid = (self.RISLine*c)/10**4
		self.freq_range_lab = (self.Harmonic*self.wavenumber_range*c)/10**4
		freq_range = Doppler_correction(self.freq_range_lab, self.Mass, self.ISCOOLEdit.value())
		offset = self.OffsetEdit.value()
		
		self.plotWidget2.plot(clear=True)
		self.plotWidget2.setLabel('bottom', "Wavenumber offset by "+str(self.OffsetEdit.value()), units='cm-1', color='#00ffff')
		self.plotWidget2.setTitle(str(self.Line)+" nm transition", color='#00ffff')
		try:
			self.plotWidget2.plot(self.wavenumber_range-offset, HFS(self.I, self.J_lower, self.J_upper, self.CF, self.A_lower, self.A_upper, self.B_lower, self.B_upper, self.FWHM, 10, 0, (freq_range-self.centroid) ), pen=tuple([0,255,255]))
		except ValueError:
			pass
	
	def updateIsotope(self):
		isotope = self.sender().currentText()
		HFS_params = np.loadtxt(self.IsotopeFilename, delimiter='\t', dtype='S')
		for i in range(len(HFS_params)):
			if HFS_params[i][0].decode("utf-8") == isotope:
				Isotope_value = HFS_params[i]
		self.I 			= float(Isotope_value[1])
		self.J_lower 	= float(Isotope_value[2])
		self.J_upper 	= float(Isotope_value[3])
		self.CF 		= float(Isotope_value[4])
		self.A_lower 	= float(Isotope_value[5])
		self.A_upper 	= float(Isotope_value[6])
		self.B_lower 	= float(Isotope_value[7])
		self.B_upper 	= float(Isotope_value[8])
		self.FWHM 		= float(Isotope_value[9])
		self.Mass 		= float(Isotope_value[10])
		self.RISLine 	= float(Isotope_value[11])
		self.Mass_ref 	= float(Isotope_value[12])
		self.ARatio		= float(Isotope_value[13])
		self.Line		= float(Isotope_value[14])
		self.Harmonic	= float(Isotope_value[15])

		self.SpinEdit.setValue(self.I)
		self.JlEdit.setValue(self.J_lower)
		self.JuEdit.setValue(self.J_upper)
		self.AlEdit.setValue(self.A_lower)
		if self.ARatioEnable.isChecked():
			self.AuEdit.setValue(self.A_lower*self.ARatio)
		else:
			self.AuEdit.setValue(self.A_upper)
		self.BlEdit.setValue(self.B_lower)
		self.BuEdit.setValue(self.B_upper)
		self.FWHMEdit.setValue(self.FWHM)
		self.CFEdit.setValue(self.CF)

		if self.old_isotope.split('-')[1] != isotope.split('-')[1]:
			centroid = self.RISLine/self.Harmonic
			alpha	= self.ISCOOLEdit.value()/(self.Mass_ref*931.494061*10**6)
			centroid_doppler = centroid/( 1 + alpha - sqrt(2*alpha + alpha*alpha))		
			offset = centroid_doppler

			minrange = round(centroid_doppler,1) - 2.0 # cm-1
			maxrange = round(centroid_doppler,1) + 2.0 # cm-1

			self.wavenumber_range = np.linspace(minrange, maxrange, 10**5)
			self.freq_range_lab = (self.Harmonic*self.wavenumber_range*c)/10**4
			self.freq_range = Doppler_correction(self.freq_range_lab, self.Mass, self.ISCOOLEdit.value())

			self.MinRangeEdit.setValue(minrange)
			self.MaxRangeEdit.setValue(maxrange)
			self.OffsetEdit.setValue(offset)
		else:
			pass

		self.old_isotope = self.IsotopeList.currentText()

	def updateHFSPeaks(self):
		for i in range(len(self.PeaksEdits)):
			self.PeaksEdits[i].setText('')
			self.IntensityEdits[i].setText('')
		PeakValues = []
		IntensityValues = []
		for i in range(len(HF_function(self.I, self.J_lower, self.J_upper, self.CF, self.A_lower, self.A_upper, self.B_lower, self.B_upper)[0])):
			PeakValue_freq = HF_function(self.I, self.J_lower, self.J_upper, self.CF, self.A_lower, self.A_upper, self.B_lower, self.B_upper)[0][i]
			IntensityValue = HF_function(self.I, self.J_lower, self.J_upper, self.CF, self.A_lower, self.A_upper, self.B_lower, self.B_upper)[1][i]

			if IntensityValue > 0:
				centroid = self.RISLine/self.Harmonic
				alpha	= self.ISCOOLEdit.value()/(self.Mass*931.494061*10**6)
				centroid_doppler = centroid/( 1 + alpha - sqrt(2*alpha + alpha*alpha))		
				offset = centroid_doppler

				PeakValue_wave = PeakValue_freq*10**4/(self.Harmonic*c) + offset
				PeakValues.append(PeakValue_wave)
				IntensityValues.append(IntensityValue)

		PeakValues, IntensityValues = np.array(PeakValues), np.array(IntensityValues)
		PeakValues_sort = sort(PeakValues)
		inds = PeakValues.argsort()
		IntensityValues_sort = IntensityValues[inds]

		self.PeaksEdits = []
		self.IntensityEdits = []
		for i in range(len(PeakValues_sort)):
			self.PeaksEdit = QtGui.QLabel(str(round(PeakValues_sort[i], 3)))
			self.grid.addWidget(self.PeaksEdit, 9, i+2)
			self.PeaksEdits.append(self.PeaksEdit)

			self.IntensityEdit = QtGui.QLabel(str(round(IntensityValues_sort[i], 1)))
			self.grid.addWidget(self.IntensityEdit, 10, i+2)
			self.IntensityEdits.append(self.IntensityEdit)

		with open("hfs_peaks.txt", 'w') as f:
			f.write(str(self.IsotopeList.currentText().split('-')[0])+'\n')
			f.write(str(self.ISCOOLEdit.value())+'\n')
			for p in PeakValues_sort:
				f.write(str(p)+"\n")	

	def updateSpin(self):
		self.I = self.sender().value()

	def updateJl(self):
		self.J_lower = self.sender().value()
		
	def updateJu(self):
		self.J_upper = self.sender().value()

	def updateAl(self):
		self.A_lower = self.sender().value()
		if self.ARatioEnable.isChecked():
			self.A_upper = self.A_lower*self.ARatio
			self.AuEdit.setValue(self.A_upper)
		
	def updateAu(self):
		self.A_upper = self.sender().value()

	def updateBl(self):
		self.B_lower = self.sender().value()
		
	def updateBu(self):
		self.B_upper = self.sender().value()

	def updateFWHM(self):
		self.FWHM = self.sender().value()

	def updateCF(self):
		self.CF = self.sender().value()
	
	def updateISCOOL(self):
		self.ISCOOL = self.sender().value()

	def updateOffset(self):
		self.Offset = self.sender().value()

	def updateRange(self):
		MinRange = self.MinRangeEdit.value()
		MaxRange = self.MaxRangeEdit.value()
		self.wavenumber_range = np.linspace(MinRange, MaxRange, 10**5)

	def calculateRanges(self):
		centroid = self.RISLine/self.Harmonic
		alpha	= self.ISCOOLEdit.value()/(self.Mass_ref*931.494061*10**6)
		centroid_doppler = centroid/( 1 + alpha - sqrt(2*alpha + alpha*alpha))		
		offset = centroid_doppler

		minrange = round(centroid_doppler,1) - 2.0 # cm-1
		maxrange = round(centroid_doppler,1) + 2.0 # cm-1

		self.wavenumber_range = np.linspace(minrange, maxrange, 10**5)
		self.freq_range_lab = (self.Harmonic*self.wavenumber_range*c)/10**4
		self.freq_range = Doppler_correction(self.freq_range_lab, self.Mass, self.ISCOOLEdit.value())

		self.MinRangeEdit.setValue(minrange)
		self.MaxRangeEdit.setValue(maxrange)
		self.OffsetEdit.setValue(offset)

	def enableARatio(self, state):
		if state == QtCore.Qt.Checked:
			self.A_upper = self.A_lower*self.ARatio
			self.AuEdit.setValue(self.A_upper)

def main():
	
	app = QtGui.QApplication(sys.argv)
	ex = Example()
	sys.exit(app.exec_())


if __name__ == '__main__':
	main()
	
	
	
	
	
	