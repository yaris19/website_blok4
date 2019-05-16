# http://biopython.org/DIST/docs/tutorial/Tutorial.html#htoc97
# http://biopython.org/DIST/docs/api/Bio.Blast.Record.HSP-class.html
# https://biopython.org/DIST/docs/api/Bio.Blast.Record-pysrc.html
# https://www.biostars.org/p/95101/
# http://biopython.org/DIST/docs/tutorial/Tutorial.html
# https://www.biostars.org/p/180510/
# https://biopython.org/DIST/docs/api/Bio.SearchIO._model.hsp.HSP-class.html

from Bio.Blast import NCBIWWW
from Bio import SearchIO
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import Entrez
import warnings
from Bio import BiopythonWarning
from datetime import datetime
import re
import mysql.connector


def main():
    print(datetime.now())
    seq_forward, seq_reverse = read_file()
    insert_database_sequence(seq_forward, seq_reverse)
    replace('read1_first50.xml', 'read1_first50_adjusted.xml')
    replace('read1_last50.xml', 'read1_last50_adjusted.xml')
    replace('read2_first50.xml', 'read2_first50_adjusted.xml')
    replace('read2_last50.xml', 'read2_last50_adjusted.xml')
    result_list1 = read_xml_file('read1_first50_adjusted.xml', True, True)
    result_list2 = read_xml_file('read1_last50_adjusted.xml', False, True)
    result_list3 = read_xml_file('read2_first50_adjusted.xml', True, False)
    result_list4 = read_xml_file('read2_last50_adjusted.xml', False, False)
    insert_database_protein(result_list1)
    insert_database_protein(result_list2)
    insert_database_protein(result_list3)
    insert_database_protein(result_list4)
    print(datetime.now())


def read_file():
    """
    Reading the CSV file and doing multiple things:
    1. Formatting FASTQ scores to decimal scores
    2. Creating fasta files with the sequences from read 1
    and read 2. Each fasta file contains 50 sequences
    :return: List with headers, sequences and decimal scores
    """
    excel_file = open('data_groep7.csv', 'r')
    # get amount of sequences
    length_file = len(excel_file.readlines())
    excel_file.close()

    excel_file = open('data_groep7.csv', 'r')
    seq_forward = []
    seq_reverse = []
    score_forward = 0
    score_reverse = 0
    counter_forward = 1
    counter_reverse = length_file + 1

    # formatting FASTQ scores to decimal scores
    # using ASCII codes, and substractig 33
    for lines in excel_file:
        line = lines.split('\t')
        for score_for in line[2]:
            score_forward += ord(score_for) - 33
        for score_rev in line[5]:
            score_reverse += ord(score_rev) - 33
        seq_forward.append(
            [line[0], line[1], round(score_forward / len(line[2]),
                                     2), counter_forward])
        seq_reverse.append(
            [line[3], line[4], round(score_reverse / len(line[5]),
                                     2), counter_reverse])

        # resetting values, so it can be calculated again
        # in the next iteration
        counter_forward += 1
        counter_reverse += 1
        score_reverse = 0
        score_forward = 0

    # writing new fasta files with the read 1 sequences
    # each fasta file contains 50 sequences
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

    # writing new fasta files with the read 2 sequences
    # each fasta file contains 50 sequences
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
    """
    Blasting the fasta files
    :param input_file_name: name of the fasta file
    :param output_file_name:  name of the XML file that will be created
    :return: an XML file with the BLAST output
    """
    print('started at', datetime.now())
    file = open(input_file_name, 'r').read()
    result_handle = NCBIWWW.qblast('blastx', 'nr', file,
                                   word_size=6, gapcosts='11 1',
                                   expect=0.0001, format_type='XML')
    out_handle = open(output_file_name, 'w')
    out_handle.write(result_handle.read())
    out_handle.close()
    print('finished at ', datetime.now())


def replace(input_file, output_file):
    file = open(input_file)
    for line in file:
        match_ap = re.search('&apos;', line)
        match_gt = re.search('&gt;', line)
        match = re.search('&apos;-', line)
        if match:
            line = line.replace(match.group(), ' ')
        if match_ap:
            line = line.replace(match_ap.group(), ' ')
        if match_gt:
            line = line.replace(match_gt.group(), '>')
        fil = open(output_file, 'a')
        fil.write(line)


def read_xml_file(file_name, first50, read1):
    """
    Parsing the XML file, and getting values. And calculating values
    based on the values in the XML file
    :param file_name: Name of fasta file that should be used
    :param first50: boolean, True if the fasta file contains the first 50 sequences,
    False if the file contains the last 50 sequences
    :param read1: boolean, True if the fasta file contains the read 1 sequnces,
    False if the file contains the read 2 sequences
    :return: a list with all the results from the XML file
    """
    warnings.simplefilter('ignore', BiopythonWarning)
    blast_results = SearchIO.parse(file_name, 'blast-xml')
    result_list = []
    counter = 0

    # amout of hits that should be saved into the database
    saving_hits = 30

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
            if counter < saving_hits:
                accession = result[res].accession
                defenition = result[res].description
                id_ = result[res][0].query_id
                bitscore = result[res][0].bitscore
                evalue = result[res][0].evalue
                ident_num = result[res][0].ident_num
                pos_num = result[res][0].pos_num
                gap_num = result[res][0].gap_num
                hit_range = [result[res][0].hit_range]
                align_len = result[res][0].hit_span
                ident_perc = round(ident_num / align_len * 100, 2)
                query_cov = round(([hit[1] - hit[0] for hit in hit_range][0])
                                  / 301 * 100, 2)
                result_list.append(
                    [sequence_id, accession, id_, defenition, bitscore,
                     evalue, ident_num, pos_num, gap_num, ident_perc,
                     query_cov])

                counter += 1
        sequence_id += 1
        counter = 0
    return result_list


def insert_database_sequence(seq_forward, seq_reverse):
    """
    Inserting the header, sequence and score into the database
    :param seq_forward: list with the headers, sequences and scores of read 1
    :param seq_reverse: list with the headers, sequences and scores of read 2
    """
    connection = mysql.connector.connect(
        host='hannl-hlo-bioinformatica-mysqlsrv.mysql.database.azure.com',
        user='fifkv@hannl-hlo-bioinformatica-mysqlsrv',
        database='fifkv',
        password='613633')
    cursor = connection.cursor()
    for seq_for in seq_forward:
        cursor.execute(
            "insert into sequence(seq_id, sequence, header, score)"
            " values ('{}', '{}', '{}', '{}')".format(seq_for[3],
                                                      seq_for[1],
                                                      seq_for[0],
                                                      seq_for[2]))
    connection.commit()

    for seq_rev in seq_reverse:
        cursor.execute(
            "insert into sequence(seq_id, sequence, header, score)"
            " values ('{}', '{}', '{}', '{}')".format(seq_rev[3],
                                                      seq_rev[1],
                                                      seq_rev[0],
                                                      seq_rev[2]))
    connection.commit()
    connection.close()


def insert_database_protein(result_list):
    """
    Inserting all the results into the database
    :param result_list: List with all the results
    """
    # name of protein is between brackets [ ]
    prot = r'\[(.*?)\]'

    # defenition of protein is untill opening bracket [
    defenition = r'.*\['

    connection = mysql.connector.connect(
        host='hannl-hlo-bioinformatica-mysqlsrv.mysql.database.azure.com',
        user='fifkv@hannl-hlo-bioinformatica-mysqlsrv',
        database='fifkv',
        password='613633')
    cursor = connection.cursor()

    # getting current length of database
    cursor.execute("select count(*) from protein")
    result = cursor.fetchone()
    amount_res = [amount for amount in result][0]
    counter = amount_res + 1

    for result in result_list:
        match_def = re.search(defenition, result[3])
        match_prot = re.search(prot, result[3])
        if match_def:
            protein = match_def.group().replace('[', '')
            cursor.execute(
                "insert into protein(name_id, defenition, "
                "accession) values ('{}', '{}', '{}')".format(counter,
                                                              protein,
                                                              result[1]))

            if match_prot:
                organism = match_prot.group().replace('[', '').replace(']', '')
                cursor.execute(
                    "insert into organism(organism_id, organism_species, "
                    "organism_genus, organism_family)"
                    "values ('{}', '{}', null, null)".format(counter,
                                                             organism))
            cursor.execute(
                "insert into protein_attribute(protein_id, seq_id, "
                "organism_id, name_id, ident_num, pos_num, gap_num, e_value, "
                "bit_score, ident_perc, query_cov) values ('{}', '{}', '{}',"
                " '{}', '{}', '{}', '{}', '{}', '{}', '{}', '{}')".format(
                    counter, result[0], counter, counter, result[6], result[7],
                    result[8], result[5], result[4], result[9], result[10]))

        connection.commit()
        counter += 1
    connection.close()


def update_database_genus_family():
    print(datetime.now())
    Entrez.email = 'yarisvanthiel@gmail.com'
    connection = mysql.connector.connect(
        host='hannl-hlo-bioinformatica-mysqlsrv.mysql.database.azure.com',
        user='fifkv@hannl-hlo-bioinformatica-mysqlsrv',
        database='fifkv',
        password='613633')
    cursor = connection.cursor()
    cursor.execute("select name_id, accession from protein")
    accession_codes = cursor.fetchall()
    for code in accession_codes:
        handle = Entrez.efetch(db='nucleotide', id=code[1], rettype='gb',
                               retmode='text')
        output_handle = SeqIO.read(handle, 'genbank')
        try:
            cursor.execute("update organism set organism_family ='{}', "
                           "organism_genus = '{}' where organism_id = '{}'"
                           "".format(output_handle.annotations['taxonomy'][-2],
                                     output_handle.annotations['taxonomy'][-1],
                                     code[0]))
        except IndexError:
            cursor.execute(("update organism set organism_family = null, "
                            "organism_genus = null where organism_id = "
                            "{}".format(code[0])))
        connection.commit()
    connection.close()
    print(datetime.now())


update_database_genus_family()
