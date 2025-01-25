import numpy as np
from PyQt5.QtWidgets import QMainWindow
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from GUI_setup import Ui_MainWindow
from Functions import prepare_register  
from Functions import run_simulation
from Functions import prepare_and_show_first_figure_correlation_matrix
from Functions import prepare_and_show_Neel_structure_factor
from Utilities import create_figure_correlation_matrix
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
        self.ui.Plot.addWidget(self.canvas)  

    def setup_connections(self):
        self.ui.Reg_seq_plot_button.clicked.connect(self.initialize_parameters)
        self.ui.Reg_seq_plot_button.clicked.connect(self.show_register_sequence)
        self.ui.Correlation_plot_button.clicked.connect(self.initialize_parameters)
        self.ui.Correlation_plot_button.clicked.connect(self.plot_correlation)
        self.ui.Neel_structure_plot_button.clicked.connect(self.initialize_parameters)
        self.ui.Neel_structure_plot_button.clicked.connect(self.show_neel_structure)
        self.ui.Show_next_command.clicked.connect(self.show_next_correlation_plot)

    def initialize_parameters(self):
        # Parameters in rad/µs and ns
        self.U = self.ui.U_h_value.value() * 2 * np.pi
        self.Omega_max = self.ui.Omega_max_value.value() * 2 * np.pi
        self.delta_0 = self.ui.Delta_initial_value.value() * 2 * np.pi
        self.delta_f = self.ui.Delta_final_value.value() * 2 * np.pi
        self.t_rise = self.ui.t_rise_value.value()
        self.t_fall = self.ui.t_fall_value.value()
        self.t_sweep_values = np.arange(self.ui.t_sweep_value_1.value(), self.ui.t_sweep_value_2.value(), NUMBER_STEPS)
        self.t_sweep_values = np.ceil(self.t_sweep_values / 4) * 4  # Ensure values are multiples of 4
        return self.U, self.Omega_max, self.delta_0, self.delta_f, self.t_rise, self.t_fall, self.t_sweep_values
    

    def show_register_sequence(self): 
        """
        Display the register sequence plot in the GUI.
        """
        U, _, _, _, _, _, _ = self.initialize_parameters()
        reg, _ = prepare_register(U)
        reg.draw()  


    
    def plot_correlation(self): 
        # Ti assicura un indice valido, percchè a ciclo finito show_next_correlation_plot() manda a un indice non valido/Managment of global variable
        global IMAGE_LIST_INDEX 
        IMAGE_LIST_INDEX = 0 

        # Simulation 
        U, Omega_max, delta_0, delta_f, t_rise, t_fall, t_sweep_range =self.initialize_parameters()
        reg, R_interatomic = prepare_register(U)
        results_sim = run_simulation(reg, R_interatomic, Omega_max, delta_0, delta_f, t_rise, t_fall, t_sweep_range)

        # Plot figure
        self.figure.clear()
        ax = self.figure.add_subplot(111)  
        prepare_and_show_first_figure_correlation_matrix(results_sim, ax=ax)  
        self.canvas.draw()  

    def show_next_correlation_plot(self):
        global IMAGE_LIST_INDEX
        IMAGE_LIST_INDEX += 1
        
        self.figure.clear()
        ax = self.figure.add_subplot(111)
        if IMAGE_LIST_INDEX < len(IMAGE_LIST):
            print(IMAGE_LIST_INDEX)
            create_figure_correlation_matrix(IMAGE_LIST[IMAGE_LIST_INDEX][0], IMAGE_LIST[IMAGE_LIST_INDEX][1], ax = ax)
            self.canvas.draw() 
        else:
            IMAGE_LIST_INDEX = -1
            print("No more plots to show, if you want to see them again, press the button again")
            self.figure.clear()
            self.canvas.draw()


    # Commento
    def show_neel_structure(self): 
        U, Omega_max, delta_0, delta_f, t_rise, t_fall, t_sweep_range =self.initialize_parameters()
        reg, R_interatomic = prepare_register(U)
        results_sim = run_simulation(reg, R_interatomic, Omega_max, delta_0, delta_f, t_rise, t_fall, t_sweep_range)

        self.figure.clear()
        ax = self.figure.add_subplot(111)  
        prepare_and_show_Neel_structure_factor(reg, R_interatomic, results_sim, ax=ax)    
        self.canvas.draw() 



