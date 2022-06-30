import os, sys
from sanger import sanger
from data_parser import sequencing, parser
from ngs import ngs


DATABASES_SEQ = sequencing.databases


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
    parser.convert_ab12fastq(ab1_folder, quality_threshold)
    parser.concatenate_fwd_rev_fasta(ab1_folder)
    sanger.consensus_pairwise_alignment(ab1_folder)
    parser.concatenate_consensus_db_identifier(ab1_folder)

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
    num_threads = 8
    # ngs_folder = os.path.join(argv[1])
    # num_threads = int(argv[2])

    '''Identify database labels'''
    # for file in ngs.natural_sort(os.listdir(ngs_folder)):
    #     DICT_READS = {}
    #     if file.__contains__('.assembled.'):
    #         output_path = os.path.join(ngs_folder, file.split('---')[0])
    #         if not os.path.exists(output_path):
    #             os.makedirs(output_path)
    #         DICT_READS = ngs.identify_database_parallel(num_threads, ngs_folder, file, DICT_READS, output_path)
    #         ngs.print_results_database(output_path, DICT_READS)

    '''Identify gRNA by database in reads'''
    for file in ngs.natural_sort(os.listdir(ngs_folder)):
        # if file.__contains__('.assembled.'):
        if file.__contains__('MI.M06648_0266.001.IDT_i7_3---IDT_i5_1.M1_3_R1.fastq.gz.assembled'):
            database_path = os.path.join(ngs_folder, file.split('---')[0])
            '''Load result from database analysis, done in the previous step. 
            Every read is associated to one of three database label: Action, Deletion or Interference.
            With that information each read is aligned against the whole database associated to it, 
            to identify the gRNA sequences with highest match.
            '''
            error, DICT_READS = ngs.load_read_from_db_result(database_path)
            print(f'Load read from db result: {str(file.split("---")[0])}')

            if DICT_READS is not None:
                DICT_READS = ngs.identify_gene_parallel(num_threads, database_path, DICT_READS)
                ngs.print_results_gene(database_path, DICT_READS)


def main(argv):
    # sanger_analysis(argv)
    ngs_analysis(argv)
    # parser.split_result_files()
    # parser.merge_results_files()
    # parser.checking_missing_data_gRNA()
    # parser.remove_duplicity_gene_results()



if __name__ == '__main__':
    main(sys.argv)
