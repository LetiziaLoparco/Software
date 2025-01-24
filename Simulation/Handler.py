import numpy as np
import sys
from PyQt5 import QtWidgets
from PyQt5.QtWidgets import QApplication, QMainWindow
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from GUI_setup import Ui_MainWindow
from Functions import prepare_register  
from Functions import run_simulation
from Functions import plot_correlation_vs_t
from Functions import plot_Neel_structure_factor_different_t
#from Functions import show_next_correlation_plot
from Functions import Create_figure_correlation_matrix
from Config import *





class GUIHandler(QMainWindow):
    def __init__(self):
        super().__init__()
        self.ui = Ui_MainWindow()
        self.ui.setupUi(self)
        self.setup_connections()

        # Create a Matplotlib figure and canvas
        self.figure = Figure()
        self.canvas = FigureCanvas(self.figure)
        self.ui.Plot.addWidget(self.canvas)  # Add canvas to QVBoxLayout

    def setup_connections(self):
        self.ui.Reg_seq_plot_button.clicked.connect(self.initialize_parameters)
        self.ui.Reg_seq_plot_button.clicked.connect(self.show_register_sequence)
        self.ui.Correlation_plot_button.clicked.connect(self.initialize_parameters)
        self.ui.Correlation_plot_button.clicked.connect(self.plot_correlation)
        self.ui.Neel_structure_plot_button.clicked.connect(self.initialize_parameters)
        self.ui.Neel_structure_plot_button.clicked.connect(self.show_neel_structure)
        self.ui.Different_neighbors_plot_button.clicked.connect(self.initialize_parameters)
        self.ui.Different_neighbors_plot_button.clicked.connect(self.show_different_neighbors)
        self.ui.Clear_command.clicked.connect(self.clear_plot)
        self.ui.Show_next_command.clicked.connect(self.show_next_correlation_plot)

    def initialize_parameters(self):
        # Parameters in rad/µs and ns
        self.U = self.ui.U_h_value.value() * 2 * np.pi
        self.Omega_max = self.ui.Omega_max_value.value() * 2 * np.pi
        self.delta_0 = self.ui.Delta_initial_value.value() * 2 * np.pi
        self.delta_f = self.ui.Delta_final_value.value() * 2 * np.pi
        self.t_rise = self.ui.t_rise_value.value()
        self.t_fall = self.ui.t_fall_value.value()
        self.t_sweep_values = np.arange(self.ui.t_sweep_value_1.value(), self.ui.t_sweep_value_2.value(), 100)
        self.t_sweep_values = np.ceil(self.t_sweep_values / 4) * 4  # Ensure values are multiples of 4
        return self.U, self.Omega_max, self.delta_0, self.delta_f, self.t_rise, self.t_fall, self.t_sweep_values
    
    def clear_plot(self):
        self.figure.clear()

    def show_register_sequence(self): 
        """
        Display the register sequence plot in the GUI.
        """
        U, _, _, _, _, _, _ = self.initialize_parameters()
        reg, _ = prepare_register(U)
        reg.draw()  


    def plot_correlation(self): # Per fare next.button crea una nuova istanza self.results_correlation
        U, Omega_max, delta_0, delta_f, t_rise, t_fall, t_sweep_range =self.initialize_parameters()
        reg, R_interatomic = prepare_register(U)
        results_correlation = run_simulation(reg, R_interatomic, Omega_max, delta_0, delta_f, t_rise, t_fall, t_sweep_range)
        # Clear the figure and redraw the plot
        self.figure.clear()
        ax = self.figure.add_subplot(111)  # Add a subplot
        plot_correlation_vs_t(results_correlation, ax=ax)  # Pass the subplot to the plot function
        #IMAGE_LIST_INDEX += 1
        self.canvas.draw()  # Refresh the canvas

    def show_next_correlation_plot(self):
        global IMAGE_LIST_INDEX
        
        self.figure.clear()
        ax = self.figure.add_subplot(111)
        if IMAGE_LIST_INDEX < len(IMAGE_LIST):
            IMAGE_LIST_INDEX += 1
            print(IMAGE_LIST_INDEX)
            print(N_SIDE)
            Create_figure_correlation_matrix(IMAGE_LIST[IMAGE_LIST_INDEX][0], IMAGE_LIST[IMAGE_LIST_INDEX][1], N_SIDE, ax = ax)
            self.canvas.draw() 

    def show_neel_structure(self): # mi dà problemi ax, got multiple values for ax
        U, Omega_max, delta_0, delta_f, t_rise, t_fall, t_sweep_range =self.initialize_parameters()
        reg, R_interatomic = prepare_register(U, N_SIDE)
        results_correlation = run_simulation(reg, R_interatomic, Omega_max, delta_0, delta_f, t_rise, t_fall, t_sweep_range)

        self.figure.clear()
        ax = self.figure.add_subplot(111)  # Add a subplot
        plot_Neel_structure_factor_different_t(reg, R_interatomic, results_correlation, ax=ax)    
        self.canvas.draw() 

    def show_different_neighbors(self):
        U, Omega_max, delta_0, delta_f, t_rise, t_fall, t_sweep_range =self.initialize_parameters()
        reg, R_interatomic = prepare_register(U, N_SIDE)
        results_correlation = run_simulation(reg, R_interatomic, Omega_max, delta_0, delta_f, t_rise, t_fall, t_sweep_range)

        self.figure.clear()
        ax = self.figure.add_subplot(111)  # Add a subplot
        plot_correlation_vs_t(results_correlation, N_SIDE, ax=ax)
        self.canvas.draw() 


        



