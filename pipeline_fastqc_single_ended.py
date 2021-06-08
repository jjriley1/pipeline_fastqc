"""====================
ReadQc pipeline
====================


This fastqc pipeline imports unmapped reads in :term:`fastq` format
from SRA accesion numbers and mapped reads in :term:`bam` format 
from TCGA, via GDC accesion numbers. Both SRA and GDC accession numbers 
are stored in remote input files.
The pipeline then performs basic quality control steps, but no processors
for quality trimming or adapter removal is used. The poor and good quality
samples are filtered/selected manually for subsequent mapping steps.

Quality metrics are based on the FastQC tools, see
http://www.bioinformatics.bbsrc.ac.uk/projects/fastqc/ for further
details.

Usage
=====

See :ref:`PipelineSettingUp` and :ref:`PipelineRunning`
on general information how to use CGAT pipelines.

When pre-processing reads before mapping, the workflow of the pipeline
is as follows:

1. Run the ``full`` target to perform initial QC on the raw data. Then
   build the report (``build_report`` target).

2. Inspect the output to decide the quality thresholds and separate poor
   and good samples for subsequent mapping steps.

3. Compare quality score results with previous stats (if available) or with
   attribute/score list from GTEx (for SRA-downloaded samples).

Configuration
-------------

See :file:`pipeline.yml` for setting configuration values affecting
the workflow (pre-processing or no pre-processing) and options for
various pre-processing tools.

Input
-----

Reads are imported by placing files or linking to files in the :term:
`working directory`.

The default file format assumes the following convention:

   <sample>.<suffix>

The ``suffix`` determines the file type. The following suffixes/file
types are possible:

remote
   Text files containing SRA or GDC accession numbers required for downloading
   unmapped reads in :term:`fastq` format and aligned data in :term:`bam` format
   files, respectively.

fastq.gz
   Single-end reads in fastq format.

fastq.1.gz, fastq.2.gz
   Paired-end reads in fastq format.
   The two fastq files must be sorted by read-pair.

bam
   Aligned/mapped reads

.. note::

   Quality scores need to be of the same scale for all input files.
   Thus it might be difficult to mix different formats.

Pipeline output
----------------

The major output is a set of HTML pages and plots reporting on the quality of
the sequence archive.

Example
=======

Example data is available at

https://www.cgat.org/downloads/public/cgatpipelines/pipeline_test_data/test_readqc.tgz

To run the example, simply unpack and untar::

   wget -qO- https://www.cgat.org/downloads/public/cgatpipelines/
             pipeline_test_data/test_readqc.tgz | tar -xvz
   cd test_readqc
   python <srcdir>/pipeline_readqc.py make full

Code
====

"""

# import ruffus
from ruffus import *

# import useful standard python modules
import sys
import os
import re
import shutil
import sqlite3
import glob

# import modules from the CGAT code collection
import CGATCore.Experiment as E
import CGATPipelines.PipelineMapping as PipelineMapping
from CGATCore import Pipeline as P
import CGATPipelines.PipelineReadqc as PipelineReadqc
import CGATCore.IOTools as IOTools
from CGATPipelines.Report import run_report
import CGAT.Sra as Sra

# load options from the config file
PARAMS = P.get_parameters(
    ["%s/pipeline.yml" % os.path.splitext(__file__)[0],
     "../pipeline.yml",
     "pipeline.yml"])

#Add parameter values from associated pipeline (pipeline_annotations.py)
PARAMS.update(P.peek_parameters(
    PARAMS["annotations_dir"],
    "genesets",
    prefix="annotations_",
    update_interface=True,
    restrict_interface=True))

# define input files and preprocessing steps
# list of acceptable input formats
INPUT_FORMATS = ["*.fastq.1.gz", "*.fastq.gz",
                 "*.sra", "*.csfasta.gz", "*.bam"]

# Regular expression to extract a track from an input file. Does not preserve
# a directory as part of the track.
REGEX_TRACK = r"(?P<track>[^/]+).(?P<suffix>fastq.1.gz|fastq.gz|sra|csfasta.gz|bam)"

SEQUENCEFILES_REGEX = r"(\S+).(?P<suffix>fastq.1.gz|fastq.gz|sra|csfasta.gz|bam)"

######################################################################################

def connect():
    '''connect to database.
    This method also attaches to helper databases.
    '''

    dbh = sqlite3.connect(PARAMS["database_url"])

    if not os.path.exists(PARAMS["annotations_database"]):
        raise ValueError(
            "can't find database '%s'" %
            PARAMS["annotations_database"])

    statement = '''ATTACH DATABASE '%s' as annotations''' % \
                (PARAMS["annotations_database"])

    cc = dbh.cursor()
    cc.execute(statement)
    cc.close()
    return dbh

################################################################################
# Download fastq (SRA) / bam (GDC) files using accession numbers in remote files
################################################################################

@follows(mkdir("temp_bams"))
@transform("input_files.dir/*.remote", 
           formatter(),
           "fastq.dir/{basename[0]}.fastq")
def downloadFiles(infiles, outfile):
    
    infile = infiles
    basefile = os.path.basename(infile)
    filename = "temp_bams/%s" % basefile
    baseoutfile = os.path.basename(outfile)
    outdir = os.path.dirname(outfile)

    if infile.endswith(".remote"):
        for line in IOTools.open_file(infile):
            repo, acc = line.strip().split("\t")[:2]
            if repo == "SRA":
                if not os.path.isfile(outfile + ".1.gz"):
                    statement = "; ".join([Sra.prefetch(acc),
                                           Sra.extract(acc, outdir)])
                    P.run(statement)
                else:
                    pass  
    
            elif repo == "GDC":
                base = os.path.splitext(basefile)               
                outfile = "bam.dir/" + base[0] + ".bam"

                token = glob.glob("gdc-user-token*")
                if len(token) > 0:
                    token = token[0]
                else:
                    token = None

                s, infile = Sra.process_remote_BAM(
                    infile, token, filename,
                    filter_bed=os.path.join(
                        PARAMS["annotations_dir"],
                        PARAMS["annotations_interface_contigs_bed"]))

                infile = " ".join(infile)
                if not os.path.isfile(outfile):
                    statement = "; ".join(
                        ["mkdir -p %(filename)s",
                         s,
                         '''cp %(infile)s %(outfile)s;
                            rm -r %(filename)s'''])
                    P.run(statement)
                else:
                    pass
          
            else:
                raise ValueError("Unknown repository: %s" % repo)
    else:
        pass
   

#######################################################################
# Run quality control
#######################################################################

@follows(downloadFiles, mkdir("fastqc.dir"))
@transform(["fastq.dir/*", "bam.dir/*", "input_files.dir/*"],
           formatter(REGEX_TRACK),
           "fastqc.dir/{basename[0]}.fastqc")
def runFastQC(infiles, outfile):
    '''run FastQC on each input file.

    check mapping qualities are in solexa format for downloaded .fastq 
    and .bam files.  Perform quality control checks on reads from
    .fastq and .bam files.

    '''

    infile = infiles
    outdir = os.path.dirname(outfile)

    if infile.endswith(".bam"):
        statement = '''fastqc --extract --outdir=%(outdir)s %(infile)s >& %(outfile)s'''

    else:
#        outfile = os.path.join(outdir, os.path.basename(outfile.split(os.extsep, 2)[1] + ".fastqc"))
        m = PipelineMapping.FastQC(nogroup=PARAMS["readqc_no_group"],
                                   outdir=outdir,
                                   qual_format=PARAMS['readqc_qual_format'])
        statement = m.build((infile,), outfile)
    
    if not os.path.isfile(outfile):
        P.run(statement)
    else:
        pass


@follows(runFastQC)
@split(runFastQC, ["fastqc_basic_statistics.tsv.gz", "fastqc_*.tsv.gz"])
def summarizeFastQC(infiles, outfiles):
    all_files = []
    for infile in infiles:
        track = P.snip(infile, ".fastqc")
        all_files.extend(glob.glob(
            os.path.join(track + "*_fastqc",
                         "fastqc_data.txt")))

    dfs = PipelineReadqc.read_fastqc(
        all_files)

    for key, df in dfs.items():
        fn = re.sub("basic_statistics", key, outfiles[0])
        E.info("writing to {}".format(fn))
        with IOTools.open_file(fn, "w") as outf:
            df.to_csv(outf, sep="\t", index=True)


@follows(summarizeFastQC)
@merge(runFastQC, "fastqc_status_summary.tsv.gz")
def buildFastQCSummaryStatus(infiles, outfile):
    '''load FastQC status summaries into a single table.'''
    PipelineReadqc.buildFastQCSummaryStatus(
        infiles,
        outfile,
        "fastqc.dir")


@follows(summarizeFastQC, buildFastQCSummaryStatus)
@jobs_limit(PARAMS.get("jobs_limit_db", 1), "db")
@transform((summarizeFastQC, buildFastQCSummaryStatus),
           suffix(".tsv.gz"), ".load")
def loadFastQC(infile, outfile):
    '''load FASTQC stats into database.'''
    P.load(infile, outfile, options="--add-index=track")


############################################################################
# Build report
############################################################################

@follows(loadFastQC, mkdir("MultiQC_report.dir"))
@originate("MultiQC_report.dir/multiqc_report.html")
def renderMultiqc(infile):
    '''build mulitqc report'''

    statement = (
        "export LANG=en_GB.UTF-8 && "
        "export LC_ALL=en_GB.UTF-8 && "
        "multiqc . -f && "
        "mv multiqc_report.html MultiQC_report.dir/")

    P.run(statement)


@follows(renderMultiqc)
@follows(mkdir("report"))
def build_report():
    '''build report from scratch.'''
    if len(os.listdir("multiqc_data") ) == 0:
        E.info("starting documentation build process from scratch")
        run_report(clean=True)
        os.rmdir("report")
    else:
        pass

@follows(build_report)
def filterFastQC():
    if not os.path.exists("fastqc.dir/poor_quality_samples/"):
        statement = ''' bash /shared/sudlab1/General/projects/UTRONs/MyFiles/scripts/filterby_qc_dups.sh '''
    else:
        pass

    P.run(statement)


#######################################################################################

@follows(runFastQC, loadFastQC,summarizeFastQC,buildFastQCSummaryStatus, 
         loadFastQC, renderMultiqc, build_report, filterFastQC)
def full():
    pass

@follows()
def filterBamsOrFastqs():
    statement = ''' bash /shared/sudlab1/General/projects/UTRONs/MyFiles/scripts/filterBamsOrFastqc_script.sh '''
    P.run(statement)

def main(argv=None):
    if argv is None:
        argv = sys.argv
    P.main(argv)

if __name__ == "__main__":
    sys.exit(P.main(sys.argv))
