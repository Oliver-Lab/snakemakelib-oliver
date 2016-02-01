import os
from jinja2 import Environment, PackageLoader

__import__('pkg_resources').declare_namespace(__name__)

OLIVER_PATH = os.path.dirname(__file__)
OLIVER_TEMPLATES_PATH = os.path.join(OLIVER_PATH, "_templates")

# Template path and templates
OliverTemplateEnv = Environment(loader = PackageLoader("snakemakelib_oliver", "_templates"))
OliverTemplateEnv.globals.update(zip=zip)
