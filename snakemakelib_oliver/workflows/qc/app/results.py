# Generating formated results
import os
import re
import shutil
from glob import glob
from collections import defaultdict

from bokeh.plotting import gridplot
import jinja2

from snakemake.report import data_uri

from snakemakelib_oliver import OLIVER_TEMPLATES_PATH
from snakemakelib_oliver.resources import OLIVER_TEMPLATES_BASE, copy_bootstrap


__all__ = ['aggregate_results', 'qc_summary_fastqc_copy_images', 'qc_fastqc_summary', 'qc_summary']

# Template environment
loader = jinja2.FileSystemLoader(OLIVER_TEMPLATES_PATH)
ENV = jinja2.Environment(loader=loader)

def aggregate_results(res):
    res.aggregate().save_aggregate_data()

def qc_summary_fastqc_copy_images(fastqc, name, _samples, d):
    images = ['adapter_content',
              'duplication_levels',
              'kmer_profiles',
              'per_base_n_content',
              'per_base_quality',
              'per_base_sequence_content',
              'per_sequence_gc_content',
              'per_sequence_quality',
              'per_tile_quality',
              'sequence_length_distribution'
              ]

    # Crawl FASTQC figures and move to static/images
    for i in images:
        for tgt in fastqc.targets[i]:
            sample = None
            for s in _samples:
                if re.search(s['SM'], tgt) is not None:
                    sample = s['SM']
                    fname = os.path.join(odir, 'static/images', sample + '_' + name + '_' + i + '.png')
                    shutil.copyfile(tgt, fname)
                    oname = os.path.join('static/images', sample + '_' + name + '_' + i + '.png')
                    d['fig'][i][sample][name] = oname


def qc_fastqc_summary(Fastqc, Fastqc_trimmed, Fastqc_bam, Fastqc_trimmed_bam):
    # Get Plots
    d = {}
    d['samples'] = set([x['SM'] for x in _samples])

    d['fig'] = defaultdict(lambda: defaultdict(dict))
    qc_summary_fastqc_copy_images(Fastqc, 'fastqc', _samples, d)
    qc_summary_fastqc_copy_images(Fastqc_trimmed, 'fastqc_trimmed', _samples, d)
    qc_summary_fastqc_copy_images(Fastqc_bam, 'fastqc_bam', _samples, d)
    qc_summary_fastqc_copy_images(Fastqc_trimmed_bam, 'fastqc_trimmed_bam', _samples, d)

    # Get data table
    #Fastqc.read_aggregate_data()
    #d['fastqc']['table'] = Fastqc.aggregate_data['summary'].to_html()

    #Fastqc_trimmed.read_aggregate_data()
    #d['fastqc_trimmed']['table'] = Fastqc_trimmed.aggregate_data['summary'].to_html()

    #Fastqc_bam.read_aggregate_data()
    #d['fastqc_bam']['table'] = Fastqc_bam.aggregate_data['summary'].to_html()

    #Fastqc_trimmed.read_aggregate_data()
    #d['fastqc_trimmed_bam']['table'] = Fastqc_trimmed_bam.aggregate_data['summary'].to_html()

    # Write the resulting html
    tp = ENV.get_template('fastqc.html')
    with open(os.path.join(odir, 'fastqc.html'), 'w') as fh:
        fh.write(tp.render(vars=d))


def qc_summary(config, input, output, Fastqc, Fastqc_trimmed, Fastqc_bam, Fastqc_trimmed_bam, Cutadapt, Align, Align_trimmed):

    # Get sample list
    global _samples
    _samples = config['_samples']

    global sample_sort
    if 'sample_sort' in config:
        sample_sort = config['sample_sort']
    else:
        sample_sort = set([x['SM'] for x in _samples]).tolist.sort()

    # Prepare directory
    global odir
    odir = os.path.dirname(output.html)
    os.makedirs(os.path.join(odir, 'static/css'), exist_ok=True)
    os.makedirs(os.path.join(odir, 'static/js'), exist_ok=True)
    os.makedirs(os.path.join(odir, 'static/images'), exist_ok=True)

    # Add bootstrap
    copy_bootstrap(os.path.join(odir, 'static'))

    # Generate FASTQC results page
    qc_fastqc_summary(Fastqc, Fastqc_trimmed, Fastqc_bam, Fastqc_trimmed_bam)

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
