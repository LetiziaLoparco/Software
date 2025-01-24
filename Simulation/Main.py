from Handler import GUIHandler
import sys
from PyQt5 import QtWidgets


if __name__ == "__main__":
    app = QtWidgets.QApplication(sys.argv)
    main_window = GUIHandler() # Initialize the GUI handler
    main_window.show()   # Show the GUI window
    sys.exit(app.exec_()) # Start the Qt event loop