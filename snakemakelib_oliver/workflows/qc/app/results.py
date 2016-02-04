# Generating formated results
import os
import re
import shutil
from glob import glob
import jinja2
from snakemake.report import data_uri
from bokeh.plotting import gridplot
from snakemakelib_oliver import OLIVER_TEMPLATES_PATH
from snakemakelib_oliver.resources import OLIVER_TEMPLATES_BASE, copy_bootstrap

__all__ = ['aggregate_results', 'qc_summary_fastqc_copy_images']

# Template environment
loader = jinja2.FileSystemLoader(OLIVER_TEMPLATES_PATH)
ENV = jinja2.Environment(loader=loader)

def aggregate_results(res):
    res.aggregate().save_aggregate_data()

def qc_summary_fastqc_copy_images(fastqc, name, _samples):
    # Crawl FASTQC figures and move to static/images
    fig = []
    for tgt in fastqc.targets['images']:
        sample = None
        for s in _samples:
            if re.search(s['SM'], tgt) is not None:
                sample = s['SM']
                break

        for img in glob(tgt):
            fname = os.path.join(odir, 'static/images', sample + '_' + name + '_' + os.path.basename(img))
            fig.append(fname)
            shutil.copyfile(img, fname)

    return {name: {'fig': fig}}

def qc_fastqc_summary(Fastqc, Fastqc_trimmed):
    pass
    # Write the resulting html
    #tp = Env.get_template('qc_workflow_report.html')
    #with open(output, 'w') as fh:
    #    fh.write(static_html(tp, tempalte_variables=d, css_raw=css_files))

def qc_summary(config, input, output, Fastqc, Fastqc_trimmed, Cutadapt, Align, Align_trimmed):
    # Get sample list
    _samples = config['_samples']

    # Prepare directory
    odir = os.path.dirname(output.html)
    os.makedirs(os.path.join(odir, 'static/css'), exist_ok=True)
    os.makedirs(os.path.join(odir, 'static/js'), exist_ok=True)
    os.makedirs(os.path.join(odir, 'static/images'), exist_ok=True)

    # Add bootstrap
    copy_bootstrap(os.path.join(odir, 'static'))

    # Generate FASTQC results page
    qc_fastqc_summary(Fastqc, Fastqc_trimmed)

    # Dictionary of figure
    d = {}

    # Copy FASTQC images to static
    d.update(qc_summary_fastqc_copy_images(Fastqc, 'fastqc_fq', _samples))
    d.update(qc_summary_fastqc_copy_images(Fastqc_trimmed, 'fastqc_trimmed', _samples))


    #Fastqc.read_aggregate_data()
    #d.update({'fastqc': {'fig': Cutadapt.plot('cutadapt')[0]}})

    # Write the resulting html
    #tp = Env.get_template('qc_workflow_report.html')
    #with open(output, 'w') as fh:
    #    fh.write(static_html(tp, tempalte_variables=d, css_raw=css_files))
