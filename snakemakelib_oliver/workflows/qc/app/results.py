# Generating formated results
import os
import re
import shutil
from glob import glob
from collections import defaultdict, OrderedDict

from bokeh.plotting import gridplot
import jinja2

from snakemake.report import data_uri

from snakemakelib_oliver import OLIVER_TEMPLATES_PATH
from snakemakelib_oliver.resources import OLIVER_TEMPLATES_BASE, copy_bootstrap


__all__ = ['aggregate_results', 'qc_summary', 'FastqcResults']

# Template environment
loader = jinja2.FileSystemLoader(OLIVER_TEMPLATES_PATH)
ENV = jinja2.Environment(loader=loader)

def aggregate_results(res):
    res.aggregate().save_aggregate_data()

class FastqcResults(object):
    """ Class to aggregate FASTQC results across different types 
    
    {section: {
        sample: {
            type: {
                caption: '',
                file: ''
    """


    def __init__(self, reportDir):
        self.odir = reportDir
        self._fastqcSections = ['Adapter Content',
                           'Duplication Levels',
                           'Kmer Profiles',
                           'Per Base N Content',
                           'Per Base Quality',
                           'Per Base Sequence Content',
                           'Per Sequence GC Content',
                           'Per Sequence Quality',
                           'Per Tile Quality',
                           'Sequence Length Distribution']

        self.fastqcSections = OrderedDict()
        self.fastqcSummary = OrderedDict()
        for section in self._fastqcSections:
            self.fastqcSections[section] = OrderedDict()

    def add_appliction(self, app):
        name = app.name
        _samples = app.units

        # Grab Summary tables
        app.read_aggregate_data('summary', index_col="SM")
        self.fastqcSummary[name] = app.aggregate_data['summary'].to_html(classes='table table-striped')

        # Store FASTQC images in a dictionary for easy templating
        for section in self.fastqcSections:
            tgt_re = app.iotargets[section][0]
            for sample in _samples:
                tgt = tgt_re.format(**sample)
                png = os.path.basename(tgt)
                new = 'static/images/{SM}_{PU}_{TYPE}_{FILE}'.format(FILE=png, TYPE=name, **sample)
                shutil.copyfile(tgt, os.path.join(self.odir, new))
                sampleKey = '{SM}_{PU}'.format(**sample)
                if sampleKey not in self.fastqcSections[section]:
                    self.fastqcSections[section][sampleKey] = OrderedDict()

                self.fastqcSections[section][sampleKey][name] = {
                        'caption': name,
                        'file': new
                        }

    def render(self):
        # Write the resulting html
        tp = ENV.get_template('fastqc.html')
        with open(os.path.join(self.odir, 'fastqc.html'), 'w') as fh:
            fh.write(tp.render(vars=self.fastqcSections, summaries=self.fastqcSummary))



def qc_summary(config, input, output, Fastqc, Fastqc_trimmed, Fastqc_bam, Fastqc_trimmed_bam, Cutadapt, Align, Align_trimmed):


    # Prepare directory
    odir = os.path.dirname(output.html)
    os.makedirs(os.path.join(odir, 'static/css'), exist_ok=True)
    os.makedirs(os.path.join(odir, 'static/js'), exist_ok=True)
    os.makedirs(os.path.join(odir, 'static/images'), exist_ok=True)

    # Add bootstrap
    copy_bootstrap(os.path.join(odir, 'static'))

    # Generate FASTQC results page
    fastqcRes = FastqcResults(odir)
    fastqcRes.add_appliction(Fastqc)
    fastqcRes.add_appliction(Fastqc_trimmed)
    fastqcRes.add_appliction(Fastqc_bam)
    fastqcRes.add_appliction(Fastqc_trimmed_bam)
    fastqcRes.render()

    # Dictionary of figure
    #d = {}

    # Copy FASTQC images to static
    #d.update(qc_summary_fastqc_copy_images(Fastqc, 'fastqc_fq', _samples))
    #d.update(qc_summary_fastqc_copy_images(Fastqc_trimmed, 'fastqc_trimmed', _samples))


    #Fastqc.read_aggregate_data()
    #d.update({'fastqc': {'fig': Cutadapt.plot('cutadapt')[0]}})

    # Write the resulting html
    #tp = Env.get_template('qc_workflow_report.html')
    #with open(output, 'w') as fh:
    #    fh.write(static_html(tp, tempalte_variables=d, css_raw=css_files))
