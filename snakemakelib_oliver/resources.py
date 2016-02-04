import os
import re
import shutil
from . import OLIVER_TEMPLATES_PATH

OLIVER_TEMPLATES_BASE = os.path.join(OLIVER_TEMPLATES_PATH, 'base.html')

def copy_bootstrap(path):
    """ Copy bootstrap files to project path. """
    os.makedirs(path, exist_ok=True)

    shutil.copyfile('static/css/bootstrap-theme.min.css', os.path.join(path, 'css/bootstrap-theme.min.css')) 
    shutil.copyfile('static/css/bootstrap.min.css', os.path.join(path, 'css/bootstrap.min.css')) 
    shutil.copyfile('static/js/bootstrap.min.js', os.path.join(path, 'js/bootstrap.min.js')) 
