from multiprocessing import freeze_support
import sys
import os #add

# VER 0.1.2
# importing from PyQt5 or PySide2 
try:
    from PyQt5 import QtGui
except:  
    from PySide2 import QtGui

from openQCM.common.architecture import Architecture,OSType
from openQCM.common.arguments import Arguments
from openQCM.common.logger import Logger as Log
from openQCM.core.constants import MinimalPython, Constants
from openQCM.ui import mainWindow

TAG = ""#"[Application]"


###############################################################################
# Main Application
###############################################################################
class OPENQCM:
    
    ###########################################################################
    # Initializing values for application
    ###########################################################################
    def __init__(self, argv=sys.argv):
          
        freeze_support()
        self._args = self._init_logger()
        self._app = QtGui.QApplication(argv)
        ##
        #if Architecture.get_os() is OSType.windows:
        '''
            Python console position and dimensions
          import win32gui
          xpos = 10
          ypos = 10
          width = 980
          length = 510
          def enumHandler(hwnd, lParam):
             if win32gui.IsWindowVisible(hwnd):
                win32gui.MoveWindow(hwnd, xpos, ypos, width, length, True)   
          win32gui.EnumWindows(enumHandler, None)
          # Set python console title 
          import ctypes
          ctypes.windll.kernel32.SetConsoleTitleW("Real-Time openQCM GUI - command line")
          # icon do the magin in windos OS solution here https://stackoverflow.com/a/12522799/4030282
          myappid = "openqcm.arbitrary.string"
          ctypes.windll.shell32.SetCurrentProcessExplicitAppUserModelID(myappid)
        '''
    
    ###########################################################################
    # Runs the application
    ###########################################################################
    def run(self):
        if Architecture.is_python_version(MinimalPython.major, minor=MinimalPython.minor):
            print(TAG,"Path:",os.path.dirname(__file__)) #add
            print('')
            print(TAG,"Application started")
            Log.i(TAG, "Application started")
            win = mainWindow.MainWindow(samples=self._args.get_user_samples())
            #win.setWindowTitle("{} - {}".format(Constants.app_title, Constants.app_version))
            #win.move(500, 20) #GUI position (x,y) on the screen 
            #win.show()
            
            # set application icon
            self._app.setWindowIcon(QtGui.QIcon('\\icon\\favicon.ico'))
            
            # VER 0.1.2
            # set style 
            # remove set style 
            # self._app.setStyle('Fusion')
            
            # the show must go on
            win.show()
            
            # VER 0.1.2
            # execute the application in PyQt or PySide style 
            if 'PyQt5' in sys.modules:
                self._app.exec()
            else: 
                self._app.exec_()
                
            print(TAG, "Finishing Application...")
            print(TAG, "Application closed")
            Log.i(TAG, "Finishing Application...\n")
            Log.i(TAG, "Application closed\n")
            win.close()
        else:
            self._fail()
        self.close()

    ###########################################################################
    # Closes application
    ###########################################################################
    def close(self):
        self._app.exit()
        Log.close()
        sys.exit()
               
    ###########################################################################
    # Initializing logger
    ###########################################################################
    @staticmethod
    def _init_logger():
        args = Arguments()
        args.create()
        args.set_user_log_level()
        return args
    
    ###########################################################################
    # Specifies the minimal Python version required
    ###########################################################################
    @staticmethod
    def _fail():
        txt = str("Application requires Python {}.{} to run".format(MinimalPython.major, MinimalPython.minor))
        print(TAG, txt)
        Log.e(TAG, txt)


if __name__ == '__main__':
    OPENQCM().run()
