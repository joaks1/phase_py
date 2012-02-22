#!/usr/bin/env python
import dendropy
import sys
import itertools
import re
import logging
logging.basicConfig(level=logging.WARNING)
_LOG = logging.getLogger('phase')

def parse_data_files(list_of_paths, schema='nexus'):
    ds = dendropy.DataSet()
    for path_name in list_of_paths:
        try:
            stream = open(path_name, 'rU')
        except IOError as e:
            sys.exit("Cannot open file '%s':\n%s\n" % (path_name, e))
        _LOG.info("Reading data from '%s'...\n" % stream.name)
        char_matrix = dendropy.DnaCharacterMatrix.get_from_stream(stream, schema = schema, preserve_underscores=True)
        char_matrix.label = path_name
        ds.add_char_matrix(char_matrix)
    return ds

def phase_sequences(data_set):
    # '-?ACBDGHKMNSRTWVYXacbdghkmnsrtwvyx'
    ambig_dict = {'M' : ['A', 'C'],
                  'R' : ['A', 'G'],
                  'W' : ['A', 'T'],
                  'S' : ['C', 'G'],
                  'Y' : ['C', 'T'],
                  'K' : ['G', 'T']}
    no_data_pattern = re.compile(r'^[?\-NX]+$')
    clean_pattern = re.compile(r'^[?\-NX]*[ACGT?\-NX]+[?\-NX]*$')
    two_state_pattern = re.compile(r'([MRWSYK])')
    two_plus_state_pattern = re.compile(r'([VHDBN])')
    ambig_pattern = re.compile(r'([MRWSYKVHDBN])')
    ds = dendropy.DataSet()
    for matrix in data_set.char_matrices:
        sys.stdout.write("\n### Parsing character matrix '%s' ###\n" % matrix.label)
        m = dendropy.DnaCharacterMatrix()
        for seq in matrix.iteritems():
            s = seq[1].symbols_as_string().upper()
            if no_data_pattern.match(s):
                sys.stdout.write("Removing '%s' from matrix '%s'\n" % (seq[0].label, matrix.label))
                continue
            elif clean_pattern.match(s):
                fasta1 = ">%s\n%s" % (seq[0].label + "a", s)
                fasta2 = ">%s\n%s" % (seq[0].label + "b", s)
                temp_matrix = dendropy.DnaCharacterMatrix()
                temp_matrix.read_from_string("%s\n%s\n" % (fasta1, fasta2), schema='fasta', data_type='dna')
                m.extend(temp_matrix)
            elif ambig_pattern.search(s):
                ambigs2 = two_state_pattern.findall(s)
                ambigs3 = two_plus_state_pattern.findall(s)
                if ambigs3:
                    sys.stdout.write("WARNING!: Greater than 2-state ambiguity found in '%s'\n" % seq[0].label)
                all_ambigs = ambig_pattern.findall(s)
                assert len(all_ambigs) == len(ambigs2) + len(ambigs3)
                total = len(all_ambigs)
                seq_list = list(s)
                indices = get_indices(seq_list, set(all_ambigs))
                message = ""
                if len(ambigs2) == total == 1:
                    states = ambig_dict[ambigs2[0]]
                    sequence1 = s.replace(ambigs2[0], states[0])
                    sequence2 = s.replace(ambigs2[0], states[1])
                    fasta1 = ">%s\n%s" % (seq[0].label + "a", sequence1)
                    fasta2 = ">%s\n%s" % (seq[0].label + "b", sequence2)
                    temp_matrix = dendropy.DnaCharacterMatrix()
                    temp_matrix.read_from_string("%s\n%s\n" % (fasta1, fasta2), schema='fasta', data_type='dna')
                    m.extend(temp_matrix)
                    message += "Found 1 ambiguity (%s) in taxon '%s'... sequences phased!\n" % (str(indices), seq[0].label)
                else:
                    fasta1 = ">%s_X\n%s" % (seq[0].label + "a", s)
                    fasta2 = ">%s_X\n%s" % (seq[0].label + "b", s)
                    temp_matrix = dendropy.DnaCharacterMatrix()
                    temp_matrix.read_from_string("%s\n%s\n" % (fasta1, fasta2), schema='fasta', data_type='dna')
                    m.extend(temp_matrix)
                    message += "Found %d ambiguities (%s) in taxon '%s'... sequences FLAGGED\n" % (total, str(indices), seq[0].label)
                sys.stdout.write(message)
            else:
                _LOG.error("Unexpected site pattern found in taxon '%s' in matrix '%s':\n%s\n" % (seq[0].label, matrix.label, s))
                sys.exit()
        m.label = matrix.label
        ds.add_char_matrix(m)
    return ds
            
def get_indices(target_list, set_of_elements):
    indices = {}
    for el in set_of_elements:
        indices[el] = []
    for index, element in enumerate(target_list):
        if element in set_of_elements:
            indices[element].append(index)
    return indices
    

def write_files(data_set, schema='nexus', extension=''):
    for matrix in data_set.char_matrices:
        path = "%s%s" % (matrix.label, extension)
        taxon_labels = matrix.taxon_set.labels()
        taxon_labels.sort()
        ts_sort = dendropy.TaxonSet(taxon_labels)
        matrix.reindex_taxa(taxon_set=ts_sort)
        matrix.write_to_path(path, schema=schema, simple=True)

if __name__ == '__main__':
    valid_schemas = ['nexus', 'phylip', 'fasta']
    from optparse import OptionParser
    usage = "Usage: %prog [ -v --schema=< nexus | phylip | fasta > ] data_file1 [ data_file2 data_file3 ...]"
    parser = OptionParser(usage = usage)
    parser.add_option("-v", "--verbose", dest="verbose", default=False, 
        action="store_true",
        help="Verbose output")
    parser.add_option("-d", "--debugging", dest="debugging", default=False, 
        action="store_true",
        help="Run in debugging mode.")
    parser.add_option("--schema", dest="schema", default="nexus", 
        action="store",
        type="string",
        help="File format. Default is nexus. Valid options: %s" % ", ".join(valid_schemas))
    (options, args) = parser.parse_args()

    if options.debugging or (options.debugging and options.verbose):
        _LOG.setLevel(logging.DEBUG)
    elif options.verbose:
        _LOG.setLevel(logging.INFO)
    else:
        _LOG.setLevel(logging.WARNING)

    if not args:
        sys.exit("Expecting a data file name")
    
    schema = options.schema
    if schema.lower() not in valid_schemas:
        sys.exit("Schema provided is not valid.\nValid options are: %s\nYou provided: %s\n" % (", ".join(valid_schemas), schema))
    
    data_set = parse_data_files(args, schema)
    
    new_data_set = phase_sequences(data_set)
    
    write_files(new_data_set, schema, extension='.phased')
    #write_files(data_set)
    
