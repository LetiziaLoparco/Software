import numpy as np
from PyQt5.QtWidgets import QMainWindow, QFileDialog
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from GUI_setup import Ui_MainWindow
from Functions import prepare_register  
from Functions import run_simulation
from Functions import prepare_and_show_first_figure_correlation_matrix
from Functions import prepare_and_show_Neel_structure_factor
from Functions import load_config
from Functions import save_config
from Utilities import create_figure_correlation_matrix
from SharedVariables import IMAGE_MATRIX_LIST, IMAGE_MATRIX_LIST_INDEX



class GUIHandler(QMainWindow): 
    def __init__(self):
        super().__init__() 
        self.ui = Ui_MainWindow()
        self.ui.setupUi(self)
        self.setup_connections()

        # Set up Matplotlib figure and canvas for plotting
        self.figure = Figure()
        self.canvas = FigureCanvas(self.figure) 
        self.ui.Plot.addWidget(self.canvas)  

        # Loading the configuration
        self.config = load_config() # If exist configuration file, it will be loaded when the GUI is opened
        self.sync_gui_with_config() # If the configuration file is not found, the default values will be used



    def setup_connections(self):
        """
        Set up connections between UI buttons and their respective handler functions.
        """

        self.ui.Reg_seq_plot_button.clicked.connect(self.initialize_parameters)
        self.ui.Reg_seq_plot_button.clicked.connect(self.show_register_sequence)
        self.ui.Correlation_plot_button.clicked.connect(self.initialize_parameters)
        self.ui.Correlation_plot_button.clicked.connect(self.plot_correlation)
        self.ui.Neel_structure_plot_button.clicked.connect(self.initialize_parameters)
        self.ui.Neel_structure_plot_button.clicked.connect(self.show_neel_structure)
        self.ui.Show_next_command.clicked.connect(self.show_next_correlation_plot)
        self.ui.Config_button.clicked.connect(self.load_config_from_button) # If the user wants to load a configuration file, the button will be clicked



    def sync_gui_with_config(self):
        """
        Sync the values in the GUI with the current configuration.
        """
        self.N_SIDE = self.config.get("N_SIDE", 3)  # Default to 3 if missing
        self.NUMBER_STEPS = self.config.get("NUMBER_STEPS", 100)  # Default to 100 if missing
        save_config(self.config)  # Save after loading



    def load_config_from_button(self):
        """
        Open a file dialog to load a configuration file and sync the GUI.
        """
        options = QFileDialog.Options()
        file_path, _ = QFileDialog.getOpenFileName(
            self, "Select Configuration File", "", "JSON Files (*.json);;All Files (*)", options=options
        )
        if file_path:
            try:
                self.config = load_config(file_path)
                self.sync_gui_with_config()
                print(f"Configuration loaded successfully from {file_path}")
            except Exception as e:
                print(f"Error loading configuration: {e}")



    def initialize_parameters(self):
        """
        Retrieve user-defined parameters from the GUI inputs and initialize them.

        Converts user inputs into physical units (e.g., radians/µs) and ensures time steps
        are adjusted to be multiples of 4, as required by the simulation.

        Returns
        -------
        tuple
            Initialized parameters including interaction strength, pulse settings,
            time sweep values and Number of atoms per side.
        """

        self.U = self.ui.U_h_value.value() * 2 * np.pi
        self.Omega_max = self.ui.Omega_max_value.value() * 2 * np.pi
        self.delta_0 = self.ui.Delta_initial_value.value() * 2 * np.pi
        self.delta_f = self.ui.Delta_final_value.value() * 2 * np.pi
        self.t_rise = self.ui.t_rise_value.value()
        self.t_fall = self.ui.t_fall_value.value()
        self.t_sweep_values = np.arange(self.ui.t_sweep_value_1.value(), self.ui.t_sweep_value_2.value(),  self.NUMBER_STEPS)
        self.t_sweep_values = np.ceil(self.t_sweep_values / 4) * 4  # Ensure values are multiples of 4
        return self.U, self.Omega_max, self.delta_0, self.delta_f, self.t_rise, self.t_fall, self.t_sweep_values, self.N_SIDE



    def show_register_sequence(self): 
        """
        Display the quantum register sequence in the GUI.

        Uses the `prepare_register` function to generate the register and draw its configuration.
        """
        U, _, _, _, _, _, _, N_SIDE= self.initialize_parameters()
        reg, _ = prepare_register(U, N_SIDE)
        reg.draw()  


    
    def plot_correlation(self): 
        """
        Plot the correlation function for the first time sweep and reset the index for navigating plots.

        Retrieves simulation results, processes them, and displays the first correlation plot in the GUI.
        """
        global IMAGE_MATRIX_LIST_INDEX 
        IMAGE_MATRIX_LIST_INDEX = 0 
        # Ensures the index is valid by resetting it, as the show_next_correlation_plot() 
        # function increments the index and may lead to an invalid value at the end of the cycle.
        # This handles the global variable management for navigation between plots.

        # Simulation 
        U, Omega_max, delta_0, delta_f, t_rise, t_fall, t_sweep_range, N_SIDE  =self.initialize_parameters()
        reg, R_interatomic = prepare_register(U, N_SIDE)
        results_sim = run_simulation(reg, R_interatomic, Omega_max, delta_0, delta_f, t_rise, t_fall, t_sweep_range, N_SIDE)

        # Plot figure
        self.figure.clear()
        ax = self.figure.add_subplot(111)  
        prepare_and_show_first_figure_correlation_matrix(results_sim, N_SIDE, ax=ax)  
        self.canvas.draw()  



    def show_next_correlation_plot(self):
        """
        Navigate and display the next correlation plot.

        If the index exceeds the available plots, reset the index and show a message.
        """
        global IMAGE_MATRIX_LIST_INDEX
        IMAGE_MATRIX_LIST_INDEX += 1
        N_SIDE = self.N_SIDE

        self.figure.clear()
        ax = self.figure.add_subplot(111)
        if IMAGE_MATRIX_LIST_INDEX < len(IMAGE_MATRIX_LIST):

            create_figure_correlation_matrix(IMAGE_MATRIX_LIST[IMAGE_MATRIX_LIST_INDEX][0], IMAGE_MATRIX_LIST[IMAGE_MATRIX_LIST_INDEX][1], N_SIDE, ax = ax)
            self.canvas.draw() 
        else:
            # Reset the index if no more plots are available
            IMAGE_MATRIX_LIST_INDEX = -1

            print("No more plots to show, if you want to see them again, press the button again")
            self.figure.clear()
            self.canvas.draw()


    
    def show_neel_structure(self): 
        """
        Plot the Néel structure factor for the given time sweeps.

        Runs the simulation, calculates the Néel structure factor, and displays it in the GUI.
        """
        # Simulation
        U, Omega_max, delta_0, delta_f, t_rise, t_fall, t_sweep_range, N_SIDE =self.initialize_parameters()
        reg, R_interatomic = prepare_register(U, N_SIDE)
        results_sim = run_simulation(reg, R_interatomic, Omega_max, delta_0, delta_f, t_rise, t_fall, t_sweep_range, N_SIDE)

        # Plot figure
        self.figure.clear()
        ax = self.figure.add_subplot(111)  
        prepare_and_show_Neel_structure_factor(reg, R_interatomic, results_sim, N_SIDE, ax=ax)    
        self.canvas.draw() 



