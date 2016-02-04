# Generating formated results
import os
import re
import shutil
from glob import glob
import jinja2
from snakemake.report import data_uri
from bokeh.plotting import gridplot
from snakemakelib_oliver import OLIVER_TEMPLATES_PATH

# Template environment
loader = jinja2.FileSystemLoader(OLIVER_TEMPLATES_PATH)
ENV = jinja2.Environment(loader=loader)

def aggregate_results(res):
    res.aggregate().save_aggregate_data()

def qc_summary_fastqc_copy_images(FASTQC, _samples):
    # Crawl FASTQC figures and move to static/images
    for grp in FASTQC.targets['images']:
        sample = None
        for s in _samples:
            if re.search(s['SM'], grp) is not None:
                sample = s['SM']
                break

        for img in glob(grp):
            shutil.copyfile(img, os.path.join(odir, 'static/images', sample + '_' + os.path.basename(img)))

def qc_summary(config, input, output, FASTQC):
    # Pull out regex and samples
    _samples = config['_samples']
    

    # Prepare directory
    odir = os.path.dirname(output.html)
    os.makedirs(os.path.join(odir, 'static/css'), exist_ok=True)
    os.makedirs(os.path.join(odir, 'static/js'), exist_ok=True)
    os.makedirs(os.path.join(odir, 'static/images'), exist_ok=True)

    # Dictionary of figure
    d = {}

    # Copy FASTQC images to static
    qc_summary_fastqc_copy_images(FASTQC, _samples)


    #FASTQC.read_aggregate_data()
    #d.update({'fastqc': {'fig': Cutadapt.plot('cutadapt')[0]}})

    # Write the resulting html
    #tp = Env.get_template('qc_workflow_report.html')
    #with open(output, 'w') as fh:
    #    fh.write(static_html(tp, tempalte_variables=d, css_raw=css_files))
