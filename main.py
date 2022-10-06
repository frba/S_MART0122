import os, sys
from sanger import sanger
from data_parser import sequencing, parser
from ngs import ngs
from nanopore import nanopore
from Bio import SeqIO, AlignIO


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
    ngs_folder = '/home/flavia/Documents/Concordia/compute_canada/NGS_DATA/'
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
    for file in parser.natural_sort(os.listdir(ngs_folder)):
        # if file.__contains__('.assembled.'):
        if file.__contains__('MI.M06648_0266.001.IDT_i7_9---IDT_i5_1.M2control_1_R1.fastq.gz.assembled'):
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


def nanopore_analysis():
    templates_folder = '/home/flavia/Documents/Concordia/Project/Nanopore/alignments/templates'
    barcodes_folder = '/home/flavia/Documents/Concordia/Project/Nanopore/alignments'
    assemblies_folder = '/home/flavia/Documents/Concordia/Project/Nanopore/assemblies'

    template_files = [f for f in parser.natural_sort(os.listdir(templates_folder)) if os.path.isfile(os.path.join(templates_folder, f))]

    template_sorted = ['W03','W02','W01','W96','W95','W94','W93','W92','W91','W90','W89','W88','W87','W86','W85','W84',
                       'W83','W82','W81','W80','W79','W78','W77','W76','W75','W74','W73','W72','W71','W70','W69','W68',
                       'W67','W66','W65','W64','W63','W62','W61','W60','W59','W58','W57','W56','W55','W54','W53','W52',
                       'W51','W50','W49','W48','W47','W46','W45','W44','W43','W42','W41','W40','W39','W38','W37','W36',
                       'W35','W34','W33','W32','W31','W30','W29','W28','W27','W26','W25','W24','W23','W22','W21','W20',
                       'W19','W18','W17','W16','W15','W14','W13','W12','W11','W10','W09','W08','W07','W06','W05','W04']

    barcode5_fasta = os.path.join(barcodes_folder, '5_barcode_sequence.fa')
    barcode3_fasta = os.path.join(barcodes_folder, '3_barcode_sequence.fa')
    assemblies_by_plate = parser.natural_sort([ f.path for f in os.scandir(assemblies_folder) if f.is_dir() ])

    # nanopore.barcode_alignment(barcode5_fasta, barcode3_fasta, assemblies_folder, assemblies_by_plate)

    # nanopore.template_alignment(templates_folder, template_files, assemblies_folder, assemblies_by_plate)
    nanopore.template_sorted_alignment(templates_folder, template_sorted, assemblies_folder, assemblies_by_plate)

def main(argv):
    # sanger_analysis(argv)
    # ngs_analysis(argv)
    # parser.split_result_files()
    # parser.merge_results_files()
    # parser.checking_missing_data_gRNA()
    # parser.remove_duplicity_gene_results()
    # parser.count_gRNA_byline()
    nanopore_analysis()



if __name__ == '__main__':
    main(sys.argv)
