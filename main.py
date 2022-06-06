import os.path
import sys

from sanger import sequencing, sanger
from ngs import ngs


DATABASES_SEQ = sequencing.databases
DICT_READS = {}


def sanger_analysis(argv):
    '''Folder with sequencing reads'''
    # ab1_folder = '/home/flavia/Downloads/G_MART0122/magic_sanger_.ab1_files'
    # num_threads = 5
    ab1_folder = os.path.join(argv[1], 'Sanger_seq_test_data')
    num_threads = int(argv[2])
    quality_threshold = 20
    # main.mainloop()

    '''gRNA remove duplicity'''
    # sanger.remove_duplicity_gRNA_db(ab1_folder)

    '''Parse files'''
    # sanger.convert_ab12fasta(ab1_folder)
    sanger.convert_ab12fastq(ab1_folder, quality_threshold)
    sanger.concatenate_fwd_rev_fasta(ab1_folder)
    sanger.consensus_pairwise_alignment(ab1_folder)
    sanger.concatenate_consensus_db_identifier(ab1_folder)

    '''Identify database labels'''
    sanger.identify_database(ab1_folder)
    sanger.print_results_database(ab1_folder)

    '''Using consensus files add each one of the database identifier'''
    error, DICT_READS = sanger.load_read_from_db_result(ab1_folder)
    if DICT_READS is not None:
        temp_path = sanger.mount_virtual_path(ab1_folder)
        sanger.identify_gene_parallel(num_threads)
        sanger.print_results_gene(ab1_folder)

    else:
        print(error)


def ngs_analysis(argv):
    ngs_folder = '/Users/flavia/Documents/NGS_project/2022_0519_NGS_results/2-Merge-QC/'
    # ngs_folder = '/Users/flavia/Documents/NGS_project/2022_0519_NGS_results/'
    num_threads = 5
    # ngs_folder = os.path.join(argv[1], 'Sanger_seq_test_data')
    # num_threads = int(argv[2])

    '''Identify database labels'''
    for file in ngs.natural_sort(os.listdir(ngs_folder)):
        if file.__contains__('.assembled.'):
            output_path = os.path.join(ngs_folder, file.split('---')[0])
            if not os.path.exists(output_path):
                os.makedirs(output_path)
            ngs.identify_database_parallel(num_threads, ngs_folder, file, output_path)
            ngs.print_results_database(output_path)

    '''Using consensus files add each one of the database identifier'''
    for file in ngs.natural_sort(os.listdir(ngs_folder)):
        if file.__contains__('.assembled.'):
            database_path = os.path.join(ngs_folder, file.split('---')[0])
            error, DICT_READS = ngs.load_read_from_db_result(database_path)
            if DICT_READS is not None:
                sanger.identify_gene_parallel(num_threads)
                sanger.print_results_gene(database_path)


def main(argv):
    # sanger_analysis(argv)
    ngs_analysis(argv)


if __name__ == '__main__':
    main(sys.argv)
