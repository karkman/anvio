#!/usr/bin/env python
"""A script to export a GenBank file from an anvi'o contigs database."""

import sys
from anvio.argparse import ArgumentParser

import anvio
import anvio.utils as utils
import anvio.terminal as terminal
import anvio.filesnpaths as filesnpaths

from anvio.errors import ConfigError, FilesNPathsError
from anvio.contigops import ExportGenbank


__copyright__ = "Copyleft 2015-2026, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__authors__ = ['karkman']
__requires__ = ['contigs-db']
__provides__ = ['genbank-file']
__description__ = "Export contigs and their features from an anvi'o contigs database as a GenBank file"


def main():
    args = get_args()
    run = terminal.Run()

    try:
        utils.is_contigs_db(args.contigs_db)
        filesnpaths.is_output_file_writable(args.output_file)

        exporter = ExportGenbank(args)
        exporter.export()

        run.info('Output GenBank', args.output_file)
    except ConfigError as e:
        print(e)
        sys.exit(-1)
    except FilesNPathsError as e:
        print(e)
        sys.exit(-1)
    except Exception as e:
        if anvio.DEBUG:
            raise
        else:
            print(f"An unexpected error occurred: {e}")
            sys.exit(-1)


def get_args():
    parser = ArgumentParser(description=__description__)

    parser.add_argument(*anvio.A('contigs-db'), **anvio.K('contigs-db', {'required': True}))
    parser.add_argument(*anvio.A('output-file'), **anvio.K('output-file', {'required': True, 'help': 'Output GenBank file path.'}))
    parser.add_argument(*anvio.A('contigs-of-interest'), **anvio.K('contigs-of-interest'))
    parser.add_argument('--annotation-sources', help='Comma-separated list of functional annotation sources to include in the GenBank file. If not provided, all sources will be included.')

    return parser.get_args(parser)


if __name__ == '__main__':
    main()
