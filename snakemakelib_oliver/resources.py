import os
import re
import shutil
from . import OLIVER_PATH, OLIVER_TEMPLATES_PATH

OLIVER_TEMPLATES_BASE = os.path.join(OLIVER_TEMPLATES_PATH, 'base.html')

def copy_bootstrap(path):
    """ Copy bootstrap files to project path. """
    os.makedirs(path, exist_ok=True)

    shutil.copyfile(os.path.join(OLIVER_PATH, 'static/css/bootstrap-theme.min.css'), os.path.join(path, 'css/bootstrap-theme.min.css')) 
    shutil.copyfile(os.path.join(OLIVER_PATH, 'static/css/bootstrap.min.css'), os.path.join(path, 'css/bootstrap.min.css')) 
    shutil.copyfile(os.path.join(OLIVER_PATH, 'static/js/bootstrap.min.js'), os.path.join(path, 'js/bootstrap.min.js')) 

    shutil.copyfile(os.path.join(OLIVER_PATH, 'static/css/dashboard.css'), os.path.join(path, 'css/dashboard.css')) 

    shutil.copyfile(os.path.join(OLIVER_PATH, 'static/js/bokeh-0.11.1.min.js'), os.path.join(path, 'js/bokeh-0.11.1.min.js')) 
    shutil.copyfile(os.path.join(OLIVER_PATH, 'static/css/bokeh-0.11.1.min.css'), os.path.join(path, 'js/bokeh-0.11.1.min.css')) 
