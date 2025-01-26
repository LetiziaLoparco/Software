from Handler import GUIHandler
import sys
from PyQt5 import QtWidgets


if __name__ == "__main__":
    app = QtWidgets.QApplication(sys.argv)
    main_window = GUIHandler() 
    main_window.show()   
    sys.exit(app.exec_()) 