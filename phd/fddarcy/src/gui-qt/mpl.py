# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'mpl.ui'
#
# Created: Sun Oct 21 20:14:38 2012
#      by: pyside-uic 0.2.14 running on PySide 1.1.1
#
# WARNING! All changes made in this file will be lost!

from PySide import QtCore, QtGui

class Ui_MainWindow(object):
    def setupUi(self, MainWindow):
        MainWindow.setObjectName("MainWindow")
        MainWindow.resize(749, 585)
        self.centralwidget = QtGui.QWidget(MainWindow)
        self.centralwidget.setObjectName("centralwidget")
        self.RUNButton = QtGui.QPushButton(self.centralwidget)
        self.RUNButton.setGeometry(QtCore.QRect(500, 470, 74, 20))
        self.RUNButton.setObjectName("RUNButton")
        self.verticalLayoutWidget = QtGui.QWidget(self.centralwidget)
        self.verticalLayoutWidget.setGeometry(QtCore.QRect(10, 20, 118, 81))
        self.verticalLayoutWidget.setObjectName("verticalLayoutWidget")
        self.verticalLayout = QtGui.QVBoxLayout(self.verticalLayoutWidget)
        self.verticalLayout.setContentsMargins(0, 0, 0, 0)
        self.verticalLayout.setObjectName("verticalLayout")
        self.splitter_2 = QtGui.QSplitter(self.verticalLayoutWidget)
        self.splitter_2.setOrientation(QtCore.Qt.Horizontal)
        self.splitter_2.setObjectName("splitter_2")
        self.label_2 = QtGui.QLabel(self.splitter_2)
        self.label_2.setObjectName("label_2")
        self.dim_x = QtGui.QLineEdit(self.splitter_2)
        self.dim_x.setObjectName("dim_x")
        self.verticalLayout.addWidget(self.splitter_2)
        self.splitter = QtGui.QSplitter(self.verticalLayoutWidget)
        self.splitter.setOrientation(QtCore.Qt.Horizontal)
        self.splitter.setObjectName("splitter")
        self.label = QtGui.QLabel(self.splitter)
        self.label.setObjectName("label")
        self.dim_y = QtGui.QLineEdit(self.splitter)
        self.dim_y.setObjectName("dim_y")
        self.verticalLayout.addWidget(self.splitter)
        self.splitter_3 = QtGui.QSplitter(self.verticalLayoutWidget)
        self.splitter_3.setOrientation(QtCore.Qt.Horizontal)
        self.splitter_3.setObjectName("splitter_3")
        self.splitter_4 = QtGui.QSplitter(self.splitter_3)
        self.splitter_4.setOrientation(QtCore.Qt.Horizontal)
        self.splitter_4.setObjectName("splitter_4")
        self.label_3 = QtGui.QLabel(self.splitter_4)
        self.label_3.setObjectName("label_3")
        self.spacing = QtGui.QLineEdit(self.splitter_4)
        self.spacing.setObjectName("spacing")
        self.verticalLayout.addWidget(self.splitter_3)
        self.verticalLayoutWidget_2 = QtGui.QWidget(self.centralwidget)
        self.verticalLayoutWidget_2.setGeometry(QtCore.QRect(10, 200, 121, 147))
        self.verticalLayoutWidget_2.setObjectName("verticalLayoutWidget_2")
        self.verticalLayout_2 = QtGui.QVBoxLayout(self.verticalLayoutWidget_2)
        self.verticalLayout_2.setContentsMargins(0, 0, 0, 0)
        self.verticalLayout_2.setObjectName("verticalLayout_2")
        self.label_5 = QtGui.QLabel(self.verticalLayoutWidget_2)
        self.label_5.setObjectName("label_5")
        self.verticalLayout_2.addWidget(self.label_5)
        self.splitter_6 = QtGui.QSplitter(self.verticalLayoutWidget_2)
        self.splitter_6.setOrientation(QtCore.Qt.Horizontal)
        self.splitter_6.setObjectName("splitter_6")
        self.label_6 = QtGui.QLabel(self.splitter_6)
        self.label_6.setObjectName("label_6")
        self.head_up = QtGui.QLineEdit(self.splitter_6)
        self.head_up.setObjectName("head_up")
        self.verticalLayout_2.addWidget(self.splitter_6)
        self.splitter_7 = QtGui.QSplitter(self.verticalLayoutWidget_2)
        self.splitter_7.setOrientation(QtCore.Qt.Horizontal)
        self.splitter_7.setObjectName("splitter_7")
        self.label_7 = QtGui.QLabel(self.splitter_7)
        self.label_7.setObjectName("label_7")
        self.head_down = QtGui.QLineEdit(self.splitter_7)
        self.head_down.setObjectName("head_down")
        self.verticalLayout_2.addWidget(self.splitter_7)
        self.splitter_9 = QtGui.QSplitter(self.verticalLayoutWidget_2)
        self.splitter_9.setOrientation(QtCore.Qt.Horizontal)
        self.splitter_9.setObjectName("splitter_9")
        self.label_9 = QtGui.QLabel(self.splitter_9)
        self.label_9.setObjectName("label_9")
        self.k = QtGui.QLineEdit(self.splitter_9)
        self.k.setObjectName("k")
        self.verticalLayout_2.addWidget(self.splitter_9)
        self.splitter_8 = QtGui.QSplitter(self.verticalLayoutWidget_2)
        self.splitter_8.setOrientation(QtCore.Qt.Horizontal)
        self.splitter_8.setObjectName("splitter_8")
        self.label_8 = QtGui.QLabel(self.splitter_8)
        self.label_8.setObjectName("label_8")
        self.n = QtGui.QLineEdit(self.splitter_8)
        self.n.setObjectName("n")
        self.verticalLayout_2.addWidget(self.splitter_8)
        self.horizontalLayout = QtGui.QHBoxLayout()
        self.horizontalLayout.setObjectName("horizontalLayout")
        self.line_bc = QtGui.QRadioButton(self.verticalLayoutWidget_2)
        self.line_bc.setEnabled(True)
        self.line_bc.setObjectName("line_bc")
        self.horizontalLayout.addWidget(self.line_bc)
        self.constant_bc = QtGui.QRadioButton(self.verticalLayoutWidget_2)
        self.constant_bc.setObjectName("constant_bc")
        self.horizontalLayout.addWidget(self.constant_bc)
        self.verticalLayout_2.addLayout(self.horizontalLayout)
        self.verticalLayoutWidget_3 = QtGui.QWidget(self.centralwidget)
        self.verticalLayoutWidget_3.setGeometry(QtCore.QRect(10, 100, 136, 94))
        self.verticalLayoutWidget_3.setObjectName("verticalLayoutWidget_3")
        self.verticalLayout_3 = QtGui.QVBoxLayout(self.verticalLayoutWidget_3)
        self.verticalLayout_3.setContentsMargins(0, 0, 0, 0)
        self.verticalLayout_3.setObjectName("verticalLayout_3")
        self.label_13 = QtGui.QLabel(self.verticalLayoutWidget_3)
        self.label_13.setObjectName("label_13")
        self.verticalLayout_3.addWidget(self.label_13)
        self.splitter_10 = QtGui.QSplitter(self.verticalLayoutWidget_3)
        self.splitter_10.setOrientation(QtCore.Qt.Horizontal)
        self.splitter_10.setObjectName("splitter_10")
        self.label_10 = QtGui.QLabel(self.splitter_10)
        self.label_10.setObjectName("label_10")
        self.max_iter_flow = QtGui.QLineEdit(self.splitter_10)
        self.max_iter_flow.setObjectName("max_iter_flow")
        self.verticalLayout_3.addWidget(self.splitter_10)
        self.splitter_5 = QtGui.QSplitter(self.verticalLayoutWidget_3)
        self.splitter_5.setOrientation(QtCore.Qt.Horizontal)
        self.splitter_5.setObjectName("splitter_5")
        self.label_4 = QtGui.QLabel(self.splitter_5)
        self.label_4.setObjectName("label_4")
        self.w_SOR = QtGui.QLineEdit(self.splitter_5)
        self.w_SOR.setObjectName("w_SOR")
        self.verticalLayout_3.addWidget(self.splitter_5)
        self.splitter_11 = QtGui.QSplitter(self.verticalLayoutWidget_3)
        self.splitter_11.setOrientation(QtCore.Qt.Horizontal)
        self.splitter_11.setObjectName("splitter_11")
        self.label_12 = QtGui.QLabel(self.splitter_11)
        self.label_12.setObjectName("label_12")
        self.convergence_limit = QtGui.QLineEdit(self.splitter_11)
        self.convergence_limit.setObjectName("convergence_limit")
        self.verticalLayout_3.addWidget(self.splitter_11)
        self.verticalLayoutWidget_4 = QtGui.QWidget(self.centralwidget)
        self.verticalLayoutWidget_4.setGeometry(QtCore.QRect(10, 350, 123, 179))
        self.verticalLayoutWidget_4.setObjectName("verticalLayoutWidget_4")
        self.verticalLayout_4 = QtGui.QVBoxLayout(self.verticalLayoutWidget_4)
        self.verticalLayout_4.setContentsMargins(0, 0, 0, 0)
        self.verticalLayout_4.setObjectName("verticalLayout_4")
        self.do_random_walk = QtGui.QCheckBox(self.verticalLayoutWidget_4)
        self.do_random_walk.setObjectName("do_random_walk")
        self.verticalLayout_4.addWidget(self.do_random_walk)
        self.splitter_12 = QtGui.QSplitter(self.verticalLayoutWidget_4)
        self.splitter_12.setOrientation(QtCore.Qt.Horizontal)
        self.splitter_12.setObjectName("splitter_12")
        self.label_14 = QtGui.QLabel(self.splitter_12)
        self.label_14.setObjectName("label_14")
        self.deltaT_rw = QtGui.QLineEdit(self.splitter_12)
        self.deltaT_rw.setObjectName("deltaT_rw")
        self.verticalLayout_4.addWidget(self.splitter_12)
        self.splitter_13 = QtGui.QSplitter(self.verticalLayoutWidget_4)
        self.splitter_13.setOrientation(QtCore.Qt.Horizontal)
        self.splitter_13.setObjectName("splitter_13")
        self.label_15 = QtGui.QLabel(self.splitter_13)
        self.label_15.setObjectName("label_15")
        self.max_iter_rw = QtGui.QLineEdit(self.splitter_13)
        self.max_iter_rw.setObjectName("max_iter_rw")
        self.verticalLayout_4.addWidget(self.splitter_13)
        self.splitter_15 = QtGui.QSplitter(self.verticalLayoutWidget_4)
        self.splitter_15.setOrientation(QtCore.Qt.Horizontal)
        self.splitter_15.setObjectName("splitter_15")
        self.label_17 = QtGui.QLabel(self.splitter_15)
        self.label_17.setObjectName("label_17")
        self.DL = QtGui.QLineEdit(self.splitter_15)
        self.DL.setObjectName("DL")
        self.verticalLayout_4.addWidget(self.splitter_15)
        self.splitter_16 = QtGui.QSplitter(self.verticalLayoutWidget_4)
        self.splitter_16.setOrientation(QtCore.Qt.Horizontal)
        self.splitter_16.setObjectName("splitter_16")
        self.label_18 = QtGui.QLabel(self.splitter_16)
        self.label_18.setObjectName("label_18")
        self.DT = QtGui.QLineEdit(self.splitter_16)
        self.DT.setObjectName("DT")
        self.verticalLayout_4.addWidget(self.splitter_16)
        self.splitter_14 = QtGui.QSplitter(self.verticalLayoutWidget_4)
        self.splitter_14.setOrientation(QtCore.Qt.Horizontal)
        self.splitter_14.setObjectName("splitter_14")
        self.label_16 = QtGui.QLabel(self.splitter_14)
        self.label_16.setObjectName("label_16")
        self.particle_num = QtGui.QLineEdit(self.splitter_14)
        self.particle_num.setObjectName("particle_num")
        self.verticalLayout_4.addWidget(self.splitter_14)
        self.splitter_18 = QtGui.QSplitter(self.verticalLayoutWidget_4)
        self.splitter_18.setOrientation(QtCore.Qt.Horizontal)
        self.splitter_18.setObjectName("splitter_18")
        self.splitter_17 = QtGui.QSplitter(self.splitter_18)
        self.splitter_17.setOrientation(QtCore.Qt.Horizontal)
        self.splitter_17.setObjectName("splitter_17")
        self.label_19 = QtGui.QLabel(self.splitter_17)
        self.label_19.setObjectName("label_19")
        self.x_rw = QtGui.QLineEdit(self.splitter_17)
        self.x_rw.setObjectName("x_rw")
        self.y_rw = QtGui.QLineEdit(self.splitter_18)
        self.y_rw.setObjectName("y_rw")
        self.verticalLayout_4.addWidget(self.splitter_18)
        self.verticalLayoutWidget_5 = QtGui.QWidget(self.centralwidget)
        self.verticalLayoutWidget_5.setGeometry(QtCore.QRect(150, 400, 160, 128))
        self.verticalLayoutWidget_5.setObjectName("verticalLayoutWidget_5")
        self.verticalLayout_5 = QtGui.QVBoxLayout(self.verticalLayoutWidget_5)
        self.verticalLayout_5.setContentsMargins(0, 0, 0, 0)
        self.verticalLayout_5.setObjectName("verticalLayout_5")
        self.checkBox_2 = QtGui.QCheckBox(self.verticalLayoutWidget_5)
        self.checkBox_2.setObjectName("checkBox_2")
        self.verticalLayout_5.addWidget(self.checkBox_2)
        self.splitter_19 = QtGui.QSplitter(self.verticalLayoutWidget_5)
        self.splitter_19.setOrientation(QtCore.Qt.Horizontal)
        self.splitter_19.setObjectName("splitter_19")
        self.label_20 = QtGui.QLabel(self.splitter_19)
        self.label_20.setObjectName("label_20")
        self.background_C = QtGui.QLineEdit(self.splitter_19)
        self.background_C.setObjectName("background_C")
        self.verticalLayout_5.addWidget(self.splitter_19)
        self.splitter_21 = QtGui.QSplitter(self.verticalLayoutWidget_5)
        self.splitter_21.setOrientation(QtCore.Qt.Horizontal)
        self.splitter_21.setObjectName("splitter_21")
        self.label_22 = QtGui.QLabel(self.splitter_21)
        self.label_22.setObjectName("label_22")
        self.initial_C = QtGui.QLineEdit(self.splitter_21)
        self.initial_C.setObjectName("initial_C")
        self.verticalLayout_5.addWidget(self.splitter_21)
        self.splitter_20 = QtGui.QSplitter(self.verticalLayoutWidget_5)
        self.splitter_20.setOrientation(QtCore.Qt.Horizontal)
        self.splitter_20.setObjectName("splitter_20")
        self.label_21 = QtGui.QLabel(self.splitter_20)
        self.label_21.setObjectName("label_21")
        self.max_iter_adv = QtGui.QLineEdit(self.splitter_20)
        self.max_iter_adv.setObjectName("max_iter_adv")
        self.verticalLayout_5.addWidget(self.splitter_20)
        self.horizontalLayout_2 = QtGui.QHBoxLayout()
        self.horizontalLayout_2.setObjectName("horizontalLayout_2")
        self.line_bc_adv = QtGui.QRadioButton(self.verticalLayoutWidget_5)
        self.line_bc_adv.setEnabled(True)
        self.line_bc_adv.setObjectName("line_bc_adv")
        self.horizontalLayout_2.addWidget(self.line_bc_adv)
        self.constant_bc_adv = QtGui.QRadioButton(self.verticalLayoutWidget_5)
        self.constant_bc_adv.setObjectName("constant_bc_adv")
        self.horizontalLayout_2.addWidget(self.constant_bc_adv)
        self.verticalLayout_5.addLayout(self.horizontalLayout_2)
        self.widget = MatplotlibWidget(self.centralwidget)
        self.widget.setGeometry(QtCore.QRect(200, 20, 501, 371))
        self.widget.setObjectName("widget")
        MainWindow.setCentralWidget(self.centralwidget)
        self.menubar = QtGui.QMenuBar(MainWindow)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 749, 20))
        self.menubar.setObjectName("menubar")
        MainWindow.setMenuBar(self.menubar)
        self.statusbar = QtGui.QStatusBar(MainWindow)
        self.statusbar.setObjectName("statusbar")
        MainWindow.setStatusBar(self.statusbar)

        self.retranslateUi(MainWindow)
        QtCore.QMetaObject.connectSlotsByName(MainWindow)

    def retranslateUi(self, MainWindow):
        MainWindow.setWindowTitle(QtGui.QApplication.translate("MainWindow", "MainWindow", None, QtGui.QApplication.UnicodeUTF8))
        self.RUNButton.setText(QtGui.QApplication.translate("MainWindow", "Run", None, QtGui.QApplication.UnicodeUTF8))
        self.label_2.setText(QtGui.QApplication.translate("MainWindow", "Dim Y", None, QtGui.QApplication.UnicodeUTF8))
        self.label.setText(QtGui.QApplication.translate("MainWindow", "Dim X", None, QtGui.QApplication.UnicodeUTF8))
        self.label_3.setText(QtGui.QApplication.translate("MainWindow", "Spacing", None, QtGui.QApplication.UnicodeUTF8))
        self.label_5.setText(QtGui.QApplication.translate("MainWindow", "Initial and Boundary C", None, QtGui.QApplication.UnicodeUTF8))
        self.label_6.setText(QtGui.QApplication.translate("MainWindow", "Head up", None, QtGui.QApplication.UnicodeUTF8))
        self.label_7.setText(QtGui.QApplication.translate("MainWindow", "Head Down", None, QtGui.QApplication.UnicodeUTF8))
        self.label_9.setText(QtGui.QApplication.translate("MainWindow", "Hid. Conduc.", None, QtGui.QApplication.UnicodeUTF8))
        self.label_8.setText(QtGui.QApplication.translate("MainWindow", "Porosity", None, QtGui.QApplication.UnicodeUTF8))
        self.line_bc.setText(QtGui.QApplication.translate("MainWindow", "Line", None, QtGui.QApplication.UnicodeUTF8))
        self.constant_bc.setText(QtGui.QApplication.translate("MainWindow", "Constant", None, QtGui.QApplication.UnicodeUTF8))
        self.label_13.setText(QtGui.QApplication.translate("MainWindow", "Flow Simulation Parameters", None, QtGui.QApplication.UnicodeUTF8))
        self.label_10.setText(QtGui.QApplication.translate("MainWindow", "Max Iter", None, QtGui.QApplication.UnicodeUTF8))
        self.label_4.setText(QtGui.QApplication.translate("MainWindow", "w SOR", None, QtGui.QApplication.UnicodeUTF8))
        self.label_12.setText(QtGui.QApplication.translate("MainWindow", "Conver. Limit", None, QtGui.QApplication.UnicodeUTF8))
        self.do_random_walk.setText(QtGui.QApplication.translate("MainWindow", "RW AD", None, QtGui.QApplication.UnicodeUTF8))
        self.label_14.setText(QtGui.QApplication.translate("MainWindow", "DeltaT", None, QtGui.QApplication.UnicodeUTF8))
        self.label_15.setText(QtGui.QApplication.translate("MainWindow", "Max Iter", None, QtGui.QApplication.UnicodeUTF8))
        self.label_17.setText(QtGui.QApplication.translate("MainWindow", "Longi. Diffusion", None, QtGui.QApplication.UnicodeUTF8))
        self.label_18.setText(QtGui.QApplication.translate("MainWindow", "Transv. Diffusion", None, QtGui.QApplication.UnicodeUTF8))
        self.label_16.setText(QtGui.QApplication.translate("MainWindow", "Particle Num", None, QtGui.QApplication.UnicodeUTF8))
        self.label_19.setText(QtGui.QApplication.translate("MainWindow", "Init Pos.", None, QtGui.QApplication.UnicodeUTF8))
        self.checkBox_2.setText(QtGui.QApplication.translate("MainWindow", "FD Advection", None, QtGui.QApplication.UnicodeUTF8))
        self.label_20.setText(QtGui.QApplication.translate("MainWindow", "Background C", None, QtGui.QApplication.UnicodeUTF8))
        self.label_22.setText(QtGui.QApplication.translate("MainWindow", "Initial C", None, QtGui.QApplication.UnicodeUTF8))
        self.label_21.setText(QtGui.QApplication.translate("MainWindow", "Max Iter", None, QtGui.QApplication.UnicodeUTF8))
        self.line_bc_adv.setText(QtGui.QApplication.translate("MainWindow", "Line", None, QtGui.QApplication.UnicodeUTF8))
        self.constant_bc_adv.setText(QtGui.QApplication.translate("MainWindow", "Constant", None, QtGui.QApplication.UnicodeUTF8))

from mplwidget import MatplotlibWidget
