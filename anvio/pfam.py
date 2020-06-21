#!/usr/bin/env python
# -*- coding: utf-8
"""
    This file contains PfamSetup and Pfam classes.

"""
import os
import gzip
import shutil
import requests
from io import BytesIO
import glob

import anvio
import anvio.dbops as dbops
import anvio.utils as utils
import anvio.terminal as terminal
import anvio.filesnpaths as filesnpaths

from anvio.drivers.hmmer import HMMer
from anvio.parsers import parser_modules
from anvio.errors import ConfigError
from anvio.tables.genefunctions import TableForGeneFunctions


__author__ = "Developers of anvi'o (see AUTHORS.txt)"
__copyright__ = "Copyleft 2015-2018, the Meren Lab (http://merenlab.org/)"
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "Özcan Esen"
__email__ = "ozcanesen@gmail.com"


run = terminal.Run()
progress = terminal.Progress()
pp = terminal.pretty_print


def read_remote_file(url, is_gzip=True):
    remote_file = requests.get(url)

    if remote_file.status_code == 404:
        raise Exception("'%s' returned 404 Not Found. " % url)

    if is_gzip:
        buf = BytesIO(remote_file.content)
        fg = gzip.GzipFile(fileobj=buf)
        return fg.read().decode('utf-8')

    return remote_file.content.decode('utf-8')


class PfamSetup(object):
    def __init__(self, args, run=run, progress=progress):
        """Setup a Pfam database for anvi'o

        Parameters
        ==========
        args : argparse.Namespace
            See `bin/anvi-setup-pfams` for available arguments
        """

        self.args = args
        self.run = run
        self.progress = progress
        self.pfam_data_dir = args.pfam_data_dir

        filesnpaths.is_program_exists('hmmpress')

        if self.pfam_data_dir and args.reset:
            raise ConfigError("You are attempting to run Pfam setup on a non-default data directory (%s) using the --reset flag. "
                              "To avoid automatically deleting a directory that may be important to you, anvi'o refuses to reset "
                              "directories that have been specified with --pfam-data-dir. If you really want to get rid of this "
                              "directory and regenerate it with Pfam data inside, then please remove the directory yourself using "
                              "a command like `rm -r %s`. We are sorry to make you go through this extra trouble, but it really is "
                              "the safest way to handle things." % (self.pfam_data_dir, self.pfam_data_dir))

        if not self.pfam_data_dir:
            self.pfam_data_dir = os.path.join(os.path.dirname(anvio.__file__), 'data/misc/Pfam')

        filesnpaths.is_output_dir_writable(os.path.dirname(os.path.abspath(self.pfam_data_dir)))

        if not args.reset and not anvio.DEBUG:
            self.is_database_exists()

        filesnpaths.gen_output_directory(self.pfam_data_dir, delete_if_exists=args.reset)

        self.resolve_database_url()
        self.files = ['Pfam-A.hmm.gz', 'Pfam.version.gz', 'Pfam-A.clans.tsv.gz']


    def resolve_database_url(self):
        page_index = 'releases/Pfam%s' % self.args.pfam_version if self.args.pfam_version else 'current_release'
        self.database_url = "http://ftp.ebi.ac.uk/pub/databases/Pfam/%s" % page_index


    def is_database_exists(self):
        if os.path.exists(os.path.join(self.pfam_data_dir, 'Pfam-A.hmm') or os.path.exists(os.path.join(self.pfam_data_dir, 'Pfam-A.hmm.gz'))):
            raise ConfigError("It seems you already have Pfam database installed in '%s', please use --reset flag if you want to re-download it." % self.pfam_data_dir)


    def get_remote_version(self):
        content = read_remote_file(self.database_url + '/Pfam.version.gz')

        # below we are parsing this, not so elegant.
        # Pfam release       : 31.0
        # Pfam-A families    : 16712
        # Date               : 2017-02
        # Based on UniProtKB : 2016_10

        version = content.strip().split('\n')[0].split(':')[1].strip()
        release_date = content.strip().split('\n')[2].split(':')[1].strip()

        self.run.info("Current Pfam version on EBI", "%s (%s)" % (version, release_date))


    def download(self):
        self.run.info("Database URL", self.database_url)

        for file_name in self.files:
            utils.download_file(self.database_url + '/' + file_name,
                os.path.join(self.pfam_data_dir, file_name), progress=self.progress, run=self.run)

        self.confirm_downloaded_files()
        self.decompress_files()


    def confirm_downloaded_files(self):
        try:
            checksums_file = read_remote_file(self.database_url + '/md5_checksums', is_gzip=False).strip()
            checksums = {}
        except:
            self.run.warning("Checksum file '%s' is not available in FTP, Anvi'o won't be able to verify downloaded files." % (self.database_url + '/md5_checksums'))
            return

        for line in checksums_file.split('\n'):
            checksum, file_name = [item.strip() for item in line.strip().split()]
            checksums[file_name] = checksum

        for file_name in self.files:
            if not filesnpaths.is_file_exists(os.path.join(self.pfam_data_dir, file_name), dont_raise=True):
                 # TO DO: Fix messages :(
                raise ConfigError("Have missing file %s, please run --reset" % file_name)

            hash_on_disk = utils.get_file_md5(os.path.join(self.pfam_data_dir, file_name))
            expected_hash = checksums[file_name]

            if not expected_hash == hash_on_disk:
                # TO DO: Fix messages :(
                raise ConfigError("Please run with --reset, one file hash doesn't match. %s" % file_name)


    def decompress_files(self):
        """Decompresses and runs hmmpress on Pfam HMM profiles."""
        for file_name in self.files:
            full_path = os.path.join(self.pfam_data_dir, file_name)

            if full_path.endswith('.gz'):
                if not os.path.exists(full_path) and os.path.exists(full_path[:-3]):
                    self.run.warning("It seems the file at %s is already decompressed. You are probably seeing "
                                     "this message because Pfams was set up previously on this computer. Hakuna Matata. Anvi'o will "
                                     "simply skip decompressing this file at this time. But if you think there is an issue, you can "
                                     "re-do the Pfam setup by running `anvi-setup-pfams` again and using the --reset flag."
                                     % (full_path[:-3]))
                    continue
                elif not os.path.exists(full_path):
                    raise ConfigError("Oh no. The file at %s does not exist. Something is terribly wrong. :( Anvi'o suggests re-running "
                                      "`anvi-setup-pfams` using the --reset flag." % (full_path))
                utils.gzip_decompress_file(full_path)
                os.remove(full_path)

        for file_path in glob.glob(os.path.join(self.pfam_data_dir, '*.hmm')):
            cmd_line = ['hmmpress', file_path]
            log_file_path = os.path.join(self.pfam_data_dir, '00_hmmpress_log.txt')
            ret_val = utils.run_command(cmd_line, log_file_path)

            if ret_val:
                raise ConfigError("Hmm. There was an error while running `hmmpress` on the Pfam HMM profiles. "
                                  "Check out the log file ('%s') to see what went wrong." % (log_file_path))
            else:
                # getting rid of the log file because hmmpress was successful
                os.remove(log_file_path)

class Pfam(object):
    def __init__(self, args, run=run, progress=progress):
        self.args = args
        self.run = run
        self.progress = progress
        self.contigs_db_path = args.contigs_db
        self.num_threads = args.num_threads
        self.hmm_program = args.hmmer_program or 'hmmsearch'
        self.pfam_data_dir = args.pfam_data_dir

        # load_catalog will populate this
        self.function_catalog = {}

        filesnpaths.is_program_exists(self.hmm_program)
        utils.is_contigs_db(self.contigs_db_path)

        if not self.pfam_data_dir:
            self.pfam_data_dir = os.path.join(os.path.dirname(anvio.__file__), 'data/misc/Pfam')

        # here, in the process of checking whether Pfam has been downloaded into the pfam_data_dir,
        # we also decompress and hmmpress the profile if it is currently gzipped
        self.is_database_exists()

        self.run.info('Pfam database directory', self.pfam_data_dir)

        self.get_version()
        self.load_catalog()


    def is_database_exists(self):
        """
        This function verifies that pfam_data_dir contains the Pfam hmm profiles and checks whether they are compressed or not.

        If they are compressed, we decompress them and run hmmpress.
        """

        if not (os.path.exists(os.path.join(self.pfam_data_dir, 'Pfam-A.hmm.gz')) or os.path.exists(os.path.join(self.pfam_data_dir, 'Pfam-A.hmm'))):
            raise ConfigError("It seems you do not have Pfam database installed, please run 'anvi-setup-pfams' to download it.")
        # here we check if the HMM profile is compressed so we can decompress it for next time
        if os.path.exists(os.path.join(self.pfam_data_dir, 'Pfam-A.hmm.gz')):
            self.run.warning("Anvi'o has detected that your Pfam database is currently compressed. It will now be unpacked before "
                             "running HMMs.")
            utils.gzip_decompress_file(os.path.join(self.pfam_data_dir, 'Pfam-A.hmm.gz'), keep_original=False)

            cmd_line = ['hmmpress', os.path.join(self.pfam_data_dir, 'Pfam-A.hmm')]
            log_file_path = os.path.join(self.pfam_data_dir, '00_hmmpress_log.txt')
            ret_val = utils.run_command(cmd_line, log_file_path)

            if ret_val:
                raise ConfigError("Hmm. There was an error while running `hmmpress` on the Pfam HMM profiles. "
                                  "Check out the log file ('%s') to see what went wrong." % (log_file_path))
            else:
                # getting rid of the log file because hmmpress was successful
                os.remove(log_file_path)


    def get_version(self):
        with open(os.path.join(self.pfam_data_dir, 'Pfam.version')) as f:
            content = f.read()

        # below we are parsing this, not so elegant.
        # Pfam release       : 31.0
        # Pfam-A families    : 16712
        # Date               : 2017-02
        # Based on UniProtKB : 2016_10

        version = content.strip().split('\n')[0].split(':')[1].strip()
        release_date = content.strip().split('\n')[2].split(':')[1].strip()

        self.run.info("Pfam database version", "%s (%s)" % (version, release_date))


    def load_catalog(self):
        catalog_path = os.path.join(self.pfam_data_dir, 'Pfam-A.clans.tsv')
        self.function_catalog = utils.get_TAB_delimited_file_as_dictionary(catalog_path,
            column_names=['accession', 'clan', 'unknown_column1', 'unknown_column2', 'function'])


    def get_function_from_catalog(self, accession, ok_if_missing_from_catalog=False):
        if '.' in accession:
            accession = accession.split('.')[0]

        if not accession in self.function_catalog:
            if ok_if_missing_from_catalog:
                return "Unkown function with PFAM accession %s" % accession
            else:
                raise ConfigError("It seems hmmscan found an accession id that does not exists "
                                  "in Pfam catalog: %s" % accession)

        return self.function_catalog[accession]['function'] # maybe merge other columns too?


    def process(self):
        hmm_file = os.path.join(self.pfam_data_dir, 'Pfam-A.hmm')

        # initialize contigs database
        class Args: pass
        args = Args()
        args.contigs_db = self.contigs_db_path
        contigs_db = dbops.ContigsSuperclass(args)
        tmp_directory_path = filesnpaths.get_temp_directory_path()

        # get an instance of gene functions table
        gene_function_calls_table = TableForGeneFunctions(self.contigs_db_path, self.run, self.progress)

        # export AA sequences for genes
        target_files_dict = {'AA:GENE': os.path.join(tmp_directory_path, 'AA_gene_sequences.fa')}
        contigs_db.gen_FASTA_file_of_sequences_for_gene_caller_ids(output_file_path=target_files_dict['AA:GENE'],
                                                                   simple_headers=True,
                                                                   rna_alphabet=False,
                                                                   report_aa_sequences=True)

        # run hmmscan
        hmmer = HMMer(target_files_dict, num_threads_to_use=self.num_threads, program_to_use=self.hmm_program)
        hmm_hits_file = hmmer.run_hmmscan('Pfam', 'AA', 'GENE', None, None, len(self.function_catalog), hmm_file, None, '--cut_ga')

        if not hmm_hits_file:
            run.info_single("The HMM search returned no hits :/ So there is nothing to add to the contigs database. But "
                            "now anvi'o will add PFAMs as a functional source with no hits, clean the temporary directories "
                            "and gracefully quit.", nl_before=1, nl_after=1)
            shutil.rmtree(tmp_directory_path)
            hmmer.clean_tmp_dirs()
            gene_function_calls_table.add_empty_sources_to_functional_sources({'Pfam'})
            return

        # parse hmmscan output
        parser = parser_modules['search']['hmmscan'](hmm_hits_file, alphabet='AA', context='GENE', program=self.hmm_program)
        search_results_dict = parser.get_search_results()

        # add functions to database
        functions_dict = {}
        counter = 0
        for hmm_hit in search_results_dict.values():
            functions_dict[counter] = {
                'gene_callers_id': hmm_hit['gene_callers_id'],
                'source': 'Pfam',
                'accession': hmm_hit['gene_hmm_id'],
                'function': self.get_function_from_catalog(hmm_hit['gene_hmm_id'], ok_if_missing_from_catalog=True),
                'e_value': hmm_hit['e_value'],
            }

            counter += 1

        if functions_dict:
            gene_function_calls_table.create(functions_dict)
        else:
            self.run.warning("Pfam class has no hits to process. Returning empty handed, but still adding Pfam as "
                             "a functional source.")
            gene_function_calls_table.add_empty_sources_to_functional_sources({'Pfam'})

        if anvio.DEBUG:
            run.warning("The temp directories, '%s' and '%s' are kept. Please don't forget to clean those up "
                        "later" % (tmp_directory_path, ', '.join(hmmer.tmp_dirs)), header="Debug")
        else:
            run.info_single('Cleaning up the temp directory (you can use `--debug` if you would '
                            'like to keep it for testing purposes)', nl_before=1, nl_after=1)
            shutil.rmtree(tmp_directory_path)
            hmmer.clean_tmp_dirs()
