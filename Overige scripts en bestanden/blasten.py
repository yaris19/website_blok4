# http://biopython.org/DIST/docs/tutorial/Tutorial.html#htoc97
# http://biopython.org/DIST/docs/api/Bio.Blast.Record.HSP-class.html
# https://biopython.org/DIST/docs/api/Bio.Blast.Record-pysrc.html
# https://www.biostars.org/p/95101/
# http://biopython.org/DIST/docs/tutorial/Tutorial.html
# https://www.biostars.org/p/180510/

from Bio.Blast import NCBIWWW
from Bio import SearchIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import warnings
from Bio import BiopythonWarning
from datetime import datetime


def main():
    seq_forward, seq_reverse = read_file()
    result_list1 = read_xml_file('read1_first50.xml', True, True)
    result_list2 = read_xml_file('read1_last50.xml', False, True)
    result_list3 = read_xml_file('read2_first50.xml', True, False)
    result_list4 = read_xml_file('read2_last50.xml', False, False)


def read_file():
    # get amount of sequences
    excel_file = open('data_groep7.csv', 'r')
    length_file = 0
    for length in excel_file:
        length_file += 1
    excel_file.close()

    excel_file = open('data_groep7.csv', 'r')
    seq_forward = []
    seq_reverse = []
    counter_forward = 1
    counter_reverse = length_file + 1
    for lines in excel_file:
        line = lines.split('\t')
        seq_forward.append([line[0], line[1], line[2], counter_forward])
        seq_reverse.append([line[3], line[4], line[5], counter_reverse])
        counter_forward += 1
        counter_reverse += 1

    file_forward_first = open('read1_sequences_first50.fasta', 'w')
    file_forward_last = open('read1_sequences_last50.fasta', 'w')
    seq_count_forward = 0
    for seq in seq_forward:
        if seq_count_forward < 50:
            record = SeqRecord(Seq(seq[1]), id=seq[0], description='')
            file_forward_first.write(record.format('fasta'))
        else:
            record = SeqRecord(Seq(seq[1]), id=seq[0], description='')
            file_forward_last.write(record.format('fasta'))
        seq_count_forward += 1
    file_forward_first.close()
    file_forward_last.close()

    file_reverse_first = open('read2_sequences_first50.fasta', 'w')
    file_reverse_last = open('read2_sequences_last50.fasta', 'w')
    seq_count_reverse = 0
    for seq in seq_reverse:
        if seq_count_reverse < 50:
            record = SeqRecord(Seq(seq[1]), id=seq[0], description='')
            file_reverse_first.write(record.format('fasta'))
        else:
            record = SeqRecord(Seq(seq[1]), id=seq[0], description='')
            file_reverse_last.write(record.format('fasta'))
        seq_count_reverse += 1

    file_reverse_first.close()
    file_reverse_last.close()
    return seq_forward, seq_reverse


def blastx(input_file_name, output_file_name):
    print('started at', datetime.now())
    file = open(input_file_name, 'r').read()
    result_handle = NCBIWWW.qblast('blastx', 'nr', file,
                                   word_size=6, gapcosts='11 1',
                                   expect=0.0001, format_type='XML')
    out_handle = open(output_file_name, 'w')
    out_handle.write(result_handle.read())
    out_handle.close()
    print('finished at ', datetime.now())


def read_xml_file(file_name, first50, read1):
    warnings.simplefilter('ignore', BiopythonWarning)
    blast_results = SearchIO.parse(file_name, 'blast-xml')
    result_list = []
    counter = 0

    sequence_id = 0
    if first50 and read1:
        sequence_id = 1
    elif not first50 and read1:
        sequence_id = 51
    elif first50 and not read1:
        sequence_id = 101
    elif not first50 and not read1:
        sequence_id = 151

    for result in blast_results:
        for i in range(len(result)):
            # print(result[1]._dict_)
            if counter < 5:
                accession = result[i].accession
                description = result[i].description
                id_ = result[i][0].query_id
                bitscore = result[i][0].bitscore
                evalue = result[i][0].evalue
                ident_num = result[i][0].ident_num
                pos_num = result[i][0].pos_num
                gap_num = result[i][0].gap_num
                result_list.append(
                    [sequence_id, accession, id_, description, bitscore,
                     evalue, ident_num, pos_num, gap_num])
                counter += 1
        sequence_id += 1
        counter = 0
    return result_list
