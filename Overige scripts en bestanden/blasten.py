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
    score_forward = 0
    score_reverse = 0
    counter_forward = 1
    counter_reverse = length_file + 1
    for lines in excel_file:
        line = lines.split('\t')
        for score_for in line[2]:
            score_forward += ord(score_for) - 33
        for score_rev in line[5]:
            score_reverse += ord(score_rev) - 33
        seq_forward.append(
            [line[0], line[1], round(score_forward / len(line[2])
                                     , 2), counter_forward])
        seq_reverse.append(
            [line[3], line[4], round(score_reverse / len(line[5])
                                     , 2), counter_reverse])
        counter_forward += 1
        counter_reverse += 1
        score_reverse = 0
        score_forward = 0

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
        for res in range(len(result)):
            # print(result[0][0].__dict__)
            if counter < 30:
                accession = result[res].accession
                description = result[res].description
                id_ = result[res][0].query_id
                bitscore = result[res][0].bitscore
                evalue = result[res][0].evalue
                ident_num = result[res][0].ident_num
                pos_num = result[res][0].pos_num
                gap_num = result[res][0].gap_num
                hit_range = [result[res][0].hit_range]
                align_len = result[res][0].hit_span
                ident_perc = round(ident_num / align_len * 100, 2)
                query_cov = round(
                    ([hit[1] - hit[0] for hit in hit_range][0]) / 301 * 100, 2)
                result_list.append(
                    [sequence_id, accession, id_, description, bitscore,
                     evalue, ident_num, pos_num, gap_num, ident_perc,
                     query_cov])

                counter += 1
        sequence_id += 1
        counter = 0
    return result_list


main()
