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
# from Bio.Alphabet import IUPAC
# from Bio.Blast import NCBIXML


def read_file():
    excel_file = open('data_groep7.csv', 'r')
    seq_forward = []
    seq_reverse = []
    for lines in excel_file:
        line = lines.split('\t')
        seq_forward.append([line[0], line[1], line[2]])
        seq_reverse.append([line[3], line[4], line[5]])

    file_forward = open('read_1_sequences.fasta', 'w')
    for seq in seq_forward:
        record = SeqRecord(Seq(seq[1]), id=seq[0], description='')
        file_forward.write(record.format('fasta'))
    file_forward.close()

    file_reverse = open('read_2_sequences.fasta', 'w')
    for seq in seq_reverse:
        record = SeqRecord(Seq(seq[1]), id=seq[0], description='')
        file_reverse.write(record.format('fasta'))
    file_reverse.close()

    return file_forward, file_reverse


def blastx():
    file_name = open('result.fasta', 'r').read()
    result_handle = NCBIWWW.qblast('blastx', 'nr', file_name,
                                   word_size=6, gapcosts='11 1',
                                   expect=0.0001)
    out_handle = open('my_blast6.xml', 'w')
    out_handle.write(result_handle.read())
    result_handle.close()
    print('finished')
    # seq_forward, seq_reverse = read_file()
    # for sequence in seq_forward:
    #     result_handle = NCBIWWW.qblast('blastx', 'nr', sequence[1],
    #                                    word_size=6, gapcosts='11 1',
    #                                    expect=0.00001)
    #     out_handle = open('my_blast4.xml', 'w')
    #     out_handle.write(result_handle.read())
    #     result_handle.close()
    #     print('finished')


def read_xml_file():
    warnings.simplefilter('ignore', BiopythonWarning)
    blast_results = SearchIO.parse('my_blast5.xml', 'blast-xml')
    result_list = []
    counter = 0
    for result in blast_results:
        # print(result.__dict__)
        for i in range(len(result)):
            # print(result[i][0].__dict__)
            if counter < 5:
                print(result[0][0].__dict__)
                accession = result[i].accession
                description = result[i].description
                id_ = result[i][0].query_id
                bitscore = result[i][0].bitscore
                evalue = result[i][0].evalue
                ident_num = result[i][0].ident_num
                pos_num = result[i][0].pos_num
                gap_num = result[i][0].gap_num
                result_list.append([accession, id_, description, bitscore,
                                    evalue,ident_num, pos_num, gap_num])
                counter += 1
        counter = 0
        # print('\n')
    return result_list


read_xml_file()

# test programs:
# def read():
#     result_handle = open('my_blast4.xml')
#     blast_record = NCBIXML.read(result_handle)
#     counter = 1
#     e_value_thresh = 0.00001
#     values = []
#     for alignment in blast_record.alignments:
#         for hsp in alignment.hsps:
#             # if float(hsp.expect) < e_value_thresh:
#             if counter >= 1000:
#                 break
#             else:
#                 identity_perc = (hsp.identities / hsp.align_length) * 100
#                 query_cover = hsp.align_length / alignment.length * 100
#                 values.append(
#                     [alignment.title, alignment.length, hsp.expect,
#                      alignment.accession, hsp.score, hsp.query, hsp.sbjct])
#                 print('***Alignment {}***'.format(counter))
#                 # print('Identity perc:', identity_perc)
#                 # print('Query cover:', query_cover)
#                 print('Sequence:', alignment.title)
#                 print('Length:', alignment.length)
#                 print('E value:', hsp.expect)
#                 print('Acession code:', alignment.accession)
#                 print(hsp.query)
#                 print(hsp.match)
#                 print(hsp.sbjct)
#                 print('\n')
#                 counter += 1
#     return values

# def read_3():
#     blast_results = SearchIO.index('my_blast5.xml', 'blast-xml')
#     for i in blast_results:
#         print(i)
#     print(blast_results)
