import os, pandas, re, numpy, csv
from Bio import SeqIO, AlignIO
from . import sequencing
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import Counter

DATABASES_SEQ = sequencing.databases
DICT_READS = {}
ARRAY_DF_GENES = []


''' Import gRNA database sequences'''
def load_gene_db(db_label):
    '''Local path Dell'''
    activation_file_path = '/home/flavia/Documents/Concordia/compute_canada/NGS-gRNA-DB/Activation_new_db.csv'
    deletion_file_path = '/home/flavia/Documents/Concordia/compute_canada/NGS-gRNA-DB/Deletion_new_db.csv'
    interference_file_path = '/home/flavia/Documents/Concordia/compute_canada/NGS-gRNA-DB/Interference_new_db.csv'

    '''Local path Mac'''
    # activation_file_path = '/Users/flavia/Documents/NGS_project/2022_0519_NGS_results/NGS-gRNA-DB/Activation_new_db_2.csv'
    # deletion_file_path = '/Users/flavia/Documents/NGS_project/2022_0519_NGS_results/NGS-gRNA-DB/Deletion_new_db_2.csv'
    # interference_file_path = '/Users/flavia/Documents/NGS_project/2022_0519_NGS_results/NGS-gRNA-DB/Interference_new_db_2.csv'

    '''Compute canada path'''
    # activation_file_path = '/home/frba/scratch/gRNA_barcode information-20220120T180349Z-001/gRNA_barcode information/Activation_new_db.csv'
    # deletion_file_path = '/home/frba/scratch/gRNA_barcode information-20220120T180349Z-001/gRNA_barcode information/Deletion_new_db.csv'
    # interference_file_path = '/home/frba/scratch/gRNA_barcode information-20220120T180349Z-001/gRNA_barcode information/Interference_new_db.csv'

    if db_label == 'Activation':
        df = pandas.read_csv(activation_file_path)
        return df
    elif db_label == 'Deletion':
        df = pandas.read_csv(deletion_file_path)
        return df
    else:
        df = pandas.read_csv(interference_file_path)
        return df


def fix_my_stuff(x):
    x = x.tolist()
    x = ', '.join([str(y) for y in x])
    return(x)


'''gRNA remove duplicity'''
def remove_duplicity_gRNA_db(ab1_folder):
    output_folder = os.path.join(ab1_folder, 'result')
    for db in DATABASES_SEQ:
        df_genes = load_gene_db(db.name)
        new_df = df_genes.groupby('Sequence').agg(lambda x: fix_my_stuff(x)).reset_index()
        new_df = new_df[['Number', 'Name', 'Sequence', 'Identifier']].sort_values(['Number'])
        new_df.to_csv(os.path.join(output_folder, str(db.name) + '_new_db.csv'), index=False)


def checking_missing_data_gRNA():
    database_path = '/home/flavia/Documents/Concordia/compute_canada/NGS_DATA/MI.M06648_0266.001.IDT_i7_9/'
    input_file = os.path.join(database_path, 'result', 'gene_result-1-incomplete.csv')
    df = pandas.read_csv(input_file).sort_values('Database')

    # test = df['Number'].isna()
    test = df.dropna(subset=['Number'])

    print(test)


def remove_duplicity_gene_results():
    database_path = '/home/flavia/Documents/Concordia/compute_canada/NGS_DATA/MI.M06648_0266.001.IDT_i7_3/'

    input_file = os.path.join(database_path, 'result', 'gene_result-3.csv')
    df = pandas.read_csv(input_file).sort_values('Database')
    test = df.drop_duplicates(subset=["Read"], inplace=True)

    print(test)


def split_result_files():
    error = ''
    database_path = '/home/flavia/Documents/Concordia/compute_canada/NGS_DATA/MI.M06648_0266.001.IDT_i7_15/'
    input_file = os.path.join(database_path, 'result', 'db_result.csv')
    if not os.path.exists(input_file):
        error = 'File with Database results not found!'
        return error, None
    else:
        df = pandas.read_csv(input_file).sort_values('Database')

    num_split = 6
    df_split = numpy.array_split(df, num_split)

    for i in range(0, num_split):
        df_split[i].to_csv(os.path.join(database_path, 'result', 'db_result-'+ str(i+1) +'.csv'), index=False)


def merge_results_files():
    # database_path = '/home/flavia/Documents/Concordia/compute_canada/NGS_DATA/MI.M06648_0266.001.IDT_i7_3/'
    database_path = '/home/flavia/Documents/Concordia/compute_canada/NGS_DATA/MI.M06648_0266.001.IDT_i7_15/'
    input_file_1 = os.path.join(database_path, 'result', 'gene_result_2022-07-01.csv')
    input_file_2 = os.path.join(database_path, 'result', 'gene_result_2022-06-30.csv')
    # input_file_3 = os.path.join(database_path, 'result', 'gene_result-3.csv')
    # input_file_4 = os.path.join(database_path, 'result', 'gene_result-1_2022-06-25_22.30.29.csv')

    df1 = pandas.read_csv(input_file_1).sort_values('Read')
    df2 = pandas.read_csv(input_file_2).sort_values('Read')
    # df3 = pandas.read_csv(input_file_3).sort_values('Read')
    # df4 = pandas.read_csv(input_file_4).sort_values('Database')

    # df = pandas.concat([df1, df2, df3]).sort_values('Read')
    df = pandas.concat([df1, df2]).sort_values('Read')
    print(len(df))
    df.to_csv(os.path.join(database_path, 'result', 'gene_result.csv'), index=False)


def count_gRNA_byline():
    for i in range(1, 2):
        database_path = f'/home/flavia/Documents/Concordia/compute_canada/NGS_DATA_2/MI.M06648_0269.001.IDT_i7_{i}/'
        col_name = ['Line', 'Number', 'Gene', 'Count', 'Database', 'Sequence']
        input_file = os.path.join(database_path, 'result', 'gene_result.csv')
        df = pandas.read_csv(input_file)

        all_df = []
        for db in DATABASES_SEQ:
            dict_db = {}
            all_gRNAs_db = []
            df_db = df[df['Database']==db.name]
            if len(df_db) > 0:
                df_genes = load_gene_db(db.name)

                # for idx, row in df_db.iterrows():
                #     ids = str(row['Number']).split(', ')
                #     for l in range(0, len(ids)):
                #         all_gRNAs_db.append(ids[l])
                # dict_count_db = Counter(all_gRNAs_db)

                for idx, row in df_db.iterrows():
                    sequence = str(row['Gene_Sequence'])
                    all_gRNAs_db.append(sequence)
                dict_count_db = Counter(all_gRNAs_db)
                print(dict_count_db)

                for idx, row in df_db.iterrows():
                    ids = str(row['Number']).split(', ')
                    genes = str(row['Gene']).split(', ')
                    sequence = str(row['Gene_Sequence'])
                    filter_df_genes = df_genes[df_genes["Sequence"] == sequence]
                    for j in range(0, len(ids)):
                        if ids[j] not in dict_db:
                            print([filter_df_genes.index.tolist()[0], ids[j], genes[j], dict_count_db[sequence], row['Database'], sequence])
                            dict_db[ids[j]] = [filter_df_genes.index.tolist()[0], ids[j], genes[j], dict_count_db[sequence], row['Database'], sequence]

                if len(dict_db) > 0:
                    new_df_db = pandas.DataFrame.from_dict(dict_db)
                    new_df_db = new_df_db.T
                    new_df_db.columns = col_name
                    if len(new_df_db) > 0:
                        all_df.append(new_df_db)

        df = pandas.concat(all_df).sort_values('Line')
        df.sort_values(['Database', 'Count'], ascending = [True, False]).to_csv(os.path.join(database_path, 'result', "count-byline_gRNA_IDT_i7_" +str(i)+".csv"),
                                                       index=False)
        print(f'MI.M06648_0266.001.IDT_i7_{i} done')


def count_gRNA():
    for i in range(14, 15):
        database_path = f'/home/flavia/Documents/Concordia/compute_canada/NGS_DATA_2/MI.M06648_0269.001.IDT_i7_{i}/'
        col_name = ['Number', 'Gene', 'Count', 'Database']
        input_file = os.path.join(database_path, 'result', 'gene_result.csv')
        df = pandas.read_csv(input_file)

        all_df = []
        for db in DATABASES_SEQ:
            dict_db = {}
            all_gRNAs_db = []
            df_db = df[df['Database']==db.name]

            for idx, row in df_db.iterrows():
                ids = str(row['Number']).split(', ')
                for l in range(0, len(ids)):
                    all_gRNAs_db.append(ids[l])
            dict_count_db = Counter(all_gRNAs_db)


            for idx, row in df_db.iterrows():
                ids = str(row['Number']).split(', ')
                genes = str(row['Gene']).split(', ')
                for j in range(0, len(ids)):
                    if ids[j] not in dict_db:
                        dict_db[ids[j]] = [ids[j], genes[j], dict_count_db[ids[j]], row['Database']]

            new_df_db = pandas.DataFrame.from_dict(dict_db)
            new_df_db = new_df_db.T
            new_df_db.columns = col_name
            if len(new_df_db) > 0:
                all_df.append(new_df_db)

        df = pandas.concat(all_df).sort_values('Number')
        df.sort_values(['Database', 'Count'], ascending = [True, False]).to_csv(os.path.join(database_path, 'result', "count_gRNA_IDT_i7_" +str(i)+".csv"),
                                                       index=False)
        print(f'MI.M06648_0266.001.IDT_i7_{i} done')


'''Parse files'''
def natural_sort(l):
    convert = lambda text: int(text) if text.isdigit() else text.lower()
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ]
    return sorted(l, key = alphanum_key)


def convert_ab12fasta(ab1_folder):
    '''Function converts ab1 files to fasta'''
    '''Check if output directory exists'''
    output_folder = os.path.join(ab1_folder, 'fasta')
    if not os.path.exists(output_folder):
        os.mkdir(output_folder)

    '''Convert ab1 to fasta using biopython'''
    for ab1_file in natural_sort(os.listdir(ab1_folder)):
        if ab1_file.endswith('.ab1'):
            try:
                records = SeqIO.parse(os.path.join(ab1_folder, ab1_file), "abi-trim")
                count = SeqIO.write(records, os.path.join(output_folder, ab1_file[:-4] + '.fa'), "fasta")
                # print("%s Converted %i records" % (ab1_file, count))
            except:
                print(ab1_file + ' could not be converted.')

def get_trim_idx_left(read_quals, quality_threshold):
    '''
    :param read_quals:
    :param quality_threshold:
    :return:
    https://cutadapt.readthedocs.io/en/stable/algorithms.html#quality-trimming-algorithm
    '''
    read_quals_reduced = [i - quality_threshold for i in read_quals]
    sum = 0
    for x in range(0, len(read_quals_reduced)):
        sum += read_quals_reduced[x]
        if sum > 0:
            return x
    return 0


def get_trim_idx_right(read_quals, quality_threshold):
    read_quals_reduced = []
    sum = 0
    for i in range(len(read_quals)-1, 0, -1):
        quality = read_quals[i] - quality_threshold
        read_quals_reduced.append(quality)
    for x in range(0, len(read_quals_reduced)):
        sum += read_quals_reduced[x]
        if sum > 0:
            idx = (len(read_quals_reduced)-1) - x
            return idx
    return len(read_quals_reduced)-1


def convert_ab12fastq(ab1_folder, quality_threshold):
    '''Files are created in output folder named fastq'''

    '''Check if output directory exists'''
    output_folder = os.path.join(ab1_folder, 'fastq')
    if not os.path.exists(output_folder):
        os.mkdir(output_folder)

    '''Check if output directory exists'''
    output_trimmed_folder = os.path.join(ab1_folder, 'fasta_trimmed')
    if not os.path.exists(output_trimmed_folder):
        os.mkdir(output_trimmed_folder)

    '''Check if quality output directory exists'''
    output_quality_trimmed_folder = os.path.join(ab1_folder, 'fasta_trimmed_quality')
    if not os.path.exists(output_quality_trimmed_folder):
        os.mkdir(output_quality_trimmed_folder)

    '''Convert ab1 to fasta using biopython'''
    for ab1_file in natural_sort(os.listdir(ab1_folder)):
        if ab1_file.endswith('.ab1'):
            try:
                records = SeqIO.parse(os.path.join(ab1_folder, ab1_file), "abi-trim")
                count = SeqIO.write(records, os.path.join(output_folder, ab1_file[:-4] + '.fastq'), "fastq")

                for rec in SeqIO.parse(os.path.join(output_folder, ab1_file[:-4] + '.fastq'), "fastq"):
                    new_rec = ''
                    read_quals = rec.letter_annotations['phred_quality']

                    l_idx_trimm = get_trim_idx_left(read_quals, quality_threshold)
                    r_idx_trimm = get_trim_idx_right(read_quals, quality_threshold)

                    for idx in range(l_idx_trimm, r_idx_trimm+1):
                        new_rec += rec.seq[idx]

                    trimmed_rec = SeqRecord(Seq(new_rec),
                                                id=rec.id,
                                                name=rec.name,
                                                description=ab1_file[:-4])

                    if len(trimmed_rec.seq) == len(read_quals[l_idx_trimm:r_idx_trimm+1]):
                        print(ab1_file, len(rec.seq), len(trimmed_rec.seq), len(read_quals[l_idx_trimm:r_idx_trimm+1]))
                    else:
                        print(ab1_file, len(rec.seq), len(trimmed_rec.seq), len(read_quals[l_idx_trimm:r_idx_trimm+1]) + ' ERROR')

                    count = SeqIO.write(trimmed_rec, os.path.join(output_trimmed_folder, ab1_file[:-4] + '.fa'), "fasta")
                    with open(os.path.join(output_quality_trimmed_folder, ab1_file[:-4] + '.qual'), "w") as quality_handle:
                        quality_handle.write(str(read_quals[l_idx_trimm:r_idx_trimm+1])[+1:-1])
                # print("%s Converted %i records" % (ab1_file, count))
            except:
                print(ab1_file + ' could not be converted.')
                # pass


def concatenate_fwd_rev_fasta(ab1_folder):
    '''Using fasta files to copy forward and reverse sequencing read in an unique fasta file
    the output files are created in the folder named pair
    Forward reads MUST BE added first in the file then reverse!!!
    MAFFT assumes the first added sequence as the correct direction.
    '''
    filename_lenght_to_compare = 14
    '''Check if output directory exists'''
    output_folder = os.path.join(ab1_folder, 'pair')
    if not os.path.exists(output_folder):
        os.mkdir(output_folder)

    fasta_folder = os.path.join(ab1_folder, 'fasta_trimmed')
    '''Mafft is affected by the order the sequences are added in file'''
    for fasta_file in natural_sort(os.listdir(fasta_folder)):
        read_pair = []
        if len(read_pair) == 0 and fasta_file.__contains__('_518_'):
            read_pair.append(fasta_file)
            for fasta_file in natural_sort(os.listdir(fasta_folder)):
                if fasta_file.__contains__('_520_') \
                        and read_pair[0][-filename_lenght_to_compare:] == fasta_file[-filename_lenght_to_compare:]:
                    read_pair.append(fasta_file)

                    fasta_read1 = SeqIO.read(open(os.path.join(ab1_folder, 'fasta_trimmed', read_pair[0])), "fasta")
                    fasta_read2 = SeqIO.read(open(os.path.join(ab1_folder, 'fasta_trimmed', read_pair[1])), "fasta")

                    with open(os.path.join(output_folder, read_pair[0][:-3] + '_' + read_pair[1][:-3] + '.fa'), 'w') \
                        as handle:
                        SeqIO.write([fasta_read1, fasta_read2], handle, 'fasta')
                    print('Concatenate files: %s - %s' % (read_pair[0], read_pair[1]))


# sanger.consensus_pairwise_alignment(ab1_folder)
def concatenate_consensus_db_identifier(ab1_folder):
    fasta_consensus_folder = os.path.join(ab1_folder, 'fasta_consensus')
    for consensus_fasta_file in natural_sort(os.listdir(fasta_consensus_folder)):
        if os.path.isfile(os.path.join(fasta_consensus_folder, consensus_fasta_file)):
            consensus_fasta_read = SeqIO.read(open(os.path.join(fasta_consensus_folder, consensus_fasta_file)), "fasta")
            for db in DATABASES_SEQ:
                output_folder = os.path.join(fasta_consensus_folder, db.name)
                if not os.path.exists(output_folder):
                    os.mkdir(output_folder)

                with open(os.path.join(output_folder, consensus_fasta_read.name + '_' + db.name + '.fa'), 'w') \
                        as handle:
                    SeqIO.write([consensus_fasta_read, db], handle, 'fasta')