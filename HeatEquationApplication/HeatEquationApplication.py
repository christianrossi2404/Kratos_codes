from KratosHeatEquationApplication import *
application = KratosHeatEquationApplication()
application_name = "KratosHeatEquationApplication"
application_folder = "HeatEquationApplication"

# The following lines are common for all applications
from . import application_importer
import inspect
caller = inspect.stack()[1] # Information about the file that imported this, to check for unexpected imports
application_importer.ImportApplication(application,application_name,application_folder,caller)
