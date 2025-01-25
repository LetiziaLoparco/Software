# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'GUI_setup.ui'
#
# Created by: PyQt5 UI code generator 5.15.11
#
# WARNING: Any manual changes made to this file will be lost when pyuic5 is
# run again.  Do not edit this file unless you know what you are doing.


from PyQt5 import QtCore, QtGui, QtWidgets


class Ui_MainWindow(object):
    def setupUi(self, MainWindow):
        MainWindow.setObjectName("MainWindow")
        MainWindow.resize(1217, 658)
        self.centralwidget = QtWidgets.QWidget(MainWindow)
        self.centralwidget.setObjectName("centralwidget")
        self.Parameters_label = QtWidgets.QLineEdit(self.centralwidget)
        self.Parameters_label.setGeometry(QtCore.QRect(670, 30, 81, 22))
        self.Parameters_label.setObjectName("Parameters_label")
        self.U_h_label = QtWidgets.QLineEdit(self.centralwidget)
        self.U_h_label.setGeometry(QtCore.QRect(670, 80, 101, 22))
        self.U_h_label.setObjectName("U_h_label")
        self.Omega_max_label = QtWidgets.QLineEdit(self.centralwidget)
        self.Omega_max_label.setGeometry(QtCore.QRect(670, 120, 101, 22))
        self.Omega_max_label.setObjectName("Omega_max_label")
        self.Delta_final_label = QtWidgets.QLineEdit(self.centralwidget)
        self.Delta_final_label.setGeometry(QtCore.QRect(670, 200, 101, 22))
        self.Delta_final_label.setObjectName("Delta_final_label")
        self.Delta_initial_label = QtWidgets.QLineEdit(self.centralwidget)
        self.Delta_initial_label.setGeometry(QtCore.QRect(670, 160, 101, 22))
        self.Delta_initial_label.setObjectName("Delta_initial_label")
        self.t_rise_label = QtWidgets.QLineEdit(self.centralwidget)
        self.t_rise_label.setGeometry(QtCore.QRect(890, 80, 101, 22))
        self.t_rise_label.setObjectName("t_rise_label")
        self.t_sweep_label = QtWidgets.QLineEdit(self.centralwidget)
        self.t_sweep_label.setGeometry(QtCore.QRect(890, 120, 101, 22))
        self.t_sweep_label.setObjectName("t_sweep_label")
        self.t_fall_label = QtWidgets.QLineEdit(self.centralwidget)
        self.t_fall_label.setGeometry(QtCore.QRect(890, 160, 101, 22))
        self.t_fall_label.setObjectName("t_fall_label")
        self.U_h_value = QtWidgets.QDoubleSpinBox(self.centralwidget)
        self.U_h_value.setGeometry(QtCore.QRect(790, 80, 61, 22))
        self.U_h_value.setLocale(QtCore.QLocale(QtCore.QLocale.English, QtCore.QLocale.UnitedStates))
        self.U_h_value.setDecimals(1)
        self.U_h_value.setSingleStep(0.1)
        self.U_h_value.setProperty("value", 2.7)
        self.U_h_value.setObjectName("U_h_value")
        self.Omega_max_value = QtWidgets.QDoubleSpinBox(self.centralwidget)
        self.Omega_max_value.setGeometry(QtCore.QRect(790, 120, 61, 22))
        self.Omega_max_value.setLocale(QtCore.QLocale(QtCore.QLocale.English, QtCore.QLocale.UnitedStates))
        self.Omega_max_value.setDecimals(1)
        self.Omega_max_value.setSingleStep(0.1)
        self.Omega_max_value.setProperty("value", 1.8)
        self.Omega_max_value.setObjectName("Omega_max_value")
        self.Delta_initial_value = QtWidgets.QDoubleSpinBox(self.centralwidget)
        self.Delta_initial_value.setGeometry(QtCore.QRect(790, 160, 61, 22))
        self.Delta_initial_value.setLocale(QtCore.QLocale(QtCore.QLocale.English, QtCore.QLocale.UnitedStates))
        self.Delta_initial_value.setDecimals(1)
        self.Delta_initial_value.setMinimum(-20.0)
        self.Delta_initial_value.setSingleStep(0.1)
        self.Delta_initial_value.setProperty("value", -6.0)
        self.Delta_initial_value.setObjectName("Delta_initial_value")
        self.Delta_final_value = QtWidgets.QDoubleSpinBox(self.centralwidget)
        self.Delta_final_value.setGeometry(QtCore.QRect(790, 200, 61, 22))
        self.Delta_final_value.setLocale(QtCore.QLocale(QtCore.QLocale.English, QtCore.QLocale.UnitedStates))
        self.Delta_final_value.setDecimals(1)
        self.Delta_final_value.setMinimum(-20.0)
        self.Delta_final_value.setSingleStep(0.1)
        self.Delta_final_value.setProperty("value", 4.5)
        self.Delta_final_value.setObjectName("Delta_final_value")
        self.t_rise_value = QtWidgets.QDoubleSpinBox(self.centralwidget)
        self.t_rise_value.setGeometry(QtCore.QRect(1010, 80, 61, 22))
        self.t_rise_value.setLocale(QtCore.QLocale(QtCore.QLocale.English, QtCore.QLocale.UnitedStates))
        self.t_rise_value.setDecimals(1)
        self.t_rise_value.setMaximum(10000.0)
        self.t_rise_value.setProperty("value", 250.0)
        self.t_rise_value.setObjectName("t_rise_value")
        self.t_sweep_value_1 = QtWidgets.QDoubleSpinBox(self.centralwidget)
        self.t_sweep_value_1.setGeometry(QtCore.QRect(1010, 120, 61, 22))
        self.t_sweep_value_1.setLocale(QtCore.QLocale(QtCore.QLocale.English, QtCore.QLocale.UnitedStates))
        self.t_sweep_value_1.setDecimals(1)
        self.t_sweep_value_1.setMaximum(10000.0)
        self.t_sweep_value_1.setProperty("value", 500.0)
        self.t_sweep_value_1.setObjectName("t_sweep_value_1")
        self.t_sweep_value_2 = QtWidgets.QDoubleSpinBox(self.centralwidget)
        self.t_sweep_value_2.setGeometry(QtCore.QRect(1130, 120, 61, 22))
        self.t_sweep_value_2.setLocale(QtCore.QLocale(QtCore.QLocale.English, QtCore.QLocale.UnitedStates))
        self.t_sweep_value_2.setDecimals(1)
        self.t_sweep_value_2.setMaximum(10000.0)
        self.t_sweep_value_2.setProperty("value", 800.0)
        self.t_sweep_value_2.setObjectName("t_sweep_value_2")
        self.t_fall_value = QtWidgets.QDoubleSpinBox(self.centralwidget)
        self.t_fall_value.setGeometry(QtCore.QRect(1010, 160, 61, 22))
        self.t_fall_value.setLocale(QtCore.QLocale(QtCore.QLocale.English, QtCore.QLocale.UnitedStates))
        self.t_fall_value.setDecimals(1)
        self.t_fall_value.setMaximum(10000.0)
        self.t_fall_value.setProperty("value", 250.0)
        self.t_fall_value.setObjectName("t_fall_value")
        self.to_label = QtWidgets.QLineEdit(self.centralwidget)
        self.to_label.setGeometry(QtCore.QRect(1090, 120, 21, 22))
        self.to_label.setObjectName("to_label")
        self.verticalLayoutWidget = QtWidgets.QWidget(self.centralwidget)
        self.verticalLayoutWidget.setGeometry(QtCore.QRect(20, 40, 621, 521))
        self.verticalLayoutWidget.setObjectName("verticalLayoutWidget")
        self.Plot = QtWidgets.QVBoxLayout(self.verticalLayoutWidget)
        self.Plot.setContentsMargins(0, 0, 0, 0)
        self.Plot.setObjectName("Plot")
        self.Correlation_plot_button = QtWidgets.QPushButton(self.centralwidget)
        self.Correlation_plot_button.setGeometry(QtCore.QRect(670, 390, 321, 28))
        self.Correlation_plot_button.setObjectName("Correlation_plot_button")
        self.Neel_structure_plot_button = QtWidgets.QPushButton(self.centralwidget)
        self.Neel_structure_plot_button.setGeometry(QtCore.QRect(670, 490, 321, 28))
        self.Neel_structure_plot_button.setObjectName("Neel_structure_plot_button")
        self.Parameters_label_2 = QtWidgets.QLineEdit(self.centralwidget)
        self.Parameters_label_2.setGeometry(QtCore.QRect(670, 330, 91, 31))
        self.Parameters_label_2.setObjectName("Parameters_label_2")
        self.Show_next_command = QtWidgets.QCommandLinkButton(self.centralwidget)
        self.Show_next_command.setGeometry(QtCore.QRect(510, 560, 131, 41))
        self.Show_next_command.setObjectName("Show_next_command")
        self.Reg_seq_plot_button = QtWidgets.QPushButton(self.centralwidget)
        self.Reg_seq_plot_button.setGeometry(QtCore.QRect(670, 250, 191, 28))
        self.Reg_seq_plot_button.setObjectName("Reg_seq_plot_button")
        MainWindow.setCentralWidget(self.centralwidget)
        self.menubar = QtWidgets.QMenuBar(MainWindow)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 1217, 26))
        self.menubar.setObjectName("menubar")
        MainWindow.setMenuBar(self.menubar)
        self.statusbar = QtWidgets.QStatusBar(MainWindow)
        self.statusbar.setObjectName("statusbar")
        MainWindow.setStatusBar(self.statusbar)

        self.retranslateUi(MainWindow)
        QtCore.QMetaObject.connectSlotsByName(MainWindow)

    def retranslateUi(self, MainWindow):
        _translate = QtCore.QCoreApplication.translate
        MainWindow.setWindowTitle(_translate("MainWindow", "MainWindow"))
        self.Parameters_label.setText(_translate("MainWindow", "Parameters:"))
        self.U_h_label.setText(_translate("MainWindow", "U/h (MHz):"))
        self.Omega_max_label.setText(_translate("MainWindow", "Ωmax/2π (MHz):"))
        self.Delta_final_label.setText(_translate("MainWindow", "δfinal/2π (MHz):"))
        self.Delta_initial_label.setText(_translate("MainWindow", "δinit/2π (MHz):"))
        self.t_rise_label.setText(_translate("MainWindow", "t rise (ns):"))
        self.t_sweep_label.setText(_translate("MainWindow", "t sweep (ns):"))
        self.t_fall_label.setText(_translate("MainWindow", "t fall (ns):"))
        self.to_label.setText(_translate("MainWindow", "to"))
        self.Correlation_plot_button.setText(_translate("MainWindow", "Show 2d antiferromagnet correlation"))
        self.Neel_structure_plot_button.setText(_translate("MainWindow", "Show Néel structure factor depending on t"))
        self.Parameters_label_2.setText(_translate("MainWindow", "Show results:"))
        self.Show_next_command.setText(_translate("MainWindow", "Show next"))
        self.Reg_seq_plot_button.setText(_translate("MainWindow", "Show Register"))


if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    MainWindow = QtWidgets.QMainWindow()
    ui = Ui_MainWindow()
    ui.setupUi(MainWindow)
    MainWindow.show()
    sys.exit(app.exec_())
