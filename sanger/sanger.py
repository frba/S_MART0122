import concurrent.futures
import os, pandas, numpy, time
from sanger import sequencing
from io import StringIO
from Bio import SeqIO, AlignIO
from Bio.Align.Applications import MafftCommandline
from Bio.Align import AlignInfo, PairwiseAligner
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from tkinter import *
from tkinter import filedialog

DATABASES_SEQ = sequencing.databases
DICT_READS = {}


def convert_ab12fasta(ab1_folder):
    '''Function converts ab1 files to fasta'''
    '''Files are created in output folder named fasta'''

    #'''Check if output directory exists'''
    output_folder = os.path.join(ab1_folder, 'fasta')
    if not os.path.exists(output_folder):
        os.mkdir(output_folder)

    #'''Convert ab1 to fasta using biopython'''
    for ab1_file in os.listdir(ab1_folder):
        if ab1_file.endswith('.ab1'):
            try:
                records = SeqIO.parse(os.path.join(ab1_folder, ab1_file), "abi-trim")
                count = SeqIO.write(records, os.path.join(output_folder, ab1_file[:-4] + '.fa'), "fasta")
                print("%s Converted %i records" % (ab1_file, count))
            except:
                print(ab1_file + ' could not be converted.')


def convert_ab12fastq(ab1_folder):
    '''Function converts ab1 files to fasta'''
    '''Files are created in output folder named fastq'''

    #'''Check if output directory exists'''
    output_folder = os.path.join(ab1_folder, 'fastq')
    if not os.path.exists(output_folder):
        os.mkdir(output_folder)

    #'''Convert ab1 to fasta using biopython'''
    for ab1_file in os.listdir(ab1_folder):
        if ab1_file.endswith('.ab1'):
            try:
                records = SeqIO.parse(os.path.join(ab1_folder, ab1_file), "abi-trim")
                count = SeqIO.write(records, os.path.join(output_folder, ab1_file[:-4] + '.fastq'), "fastq")
                print("%s Converted %i records" % (ab1_file, count))
            except:
                print(ab1_file + ' could not be converted.')


def concatenate_fwd_rev_fasta(ab1_folder):
    '''Using fasta files to copy forward and reverse sequencing read in an unique fasta file
    the output files are created in the folder named pair'''
    #'''check if output directory exists'''
    output_folder = os.path.join(ab1_folder, 'pair')
    if not os.path.exists(output_folder):
        os.mkdir(output_folder)

    fasta_folder = os.path.join(ab1_folder, 'fasta')
    total_pair_reads = int(len(os.listdir(fasta_folder)) / 2) + 1

    for i in range(1, total_pair_reads):
        file_label = 'M' + str(i) + '_'
        read_pair = []
        for fasta_file in os.listdir(fasta_folder):
            if fasta_file.__contains__(file_label):
                read_pair.append(fasta_file)

        fasta_read1 = SeqIO.read(open(os.path.join(ab1_folder, 'fasta', read_pair[0])), "fasta")
        fasta_read2 = SeqIO.read(open(os.path.join(ab1_folder, 'fasta', read_pair[1])), "fasta")

        with open(os.path.join(output_folder, fasta_read1.name + '_' + fasta_read2.name + '.fa'), 'w') \
                as handle:
            SeqIO.write([fasta_read1, fasta_read2], handle, 'fasta')
            print('Concatenate files: %s - %s' % (read_pair[0], read_pair[1]))


def score_alignment(read_seq, identifier_seq):
    aligner = PairwiseAligner()
    aligner.mode = 'local'
    aligner.match_score = 1
    aligner.mismatch_score = -1
    aligner.open_gap_score = -1
    aligner.extend_gap_score = -.1
    aligner.query_end_gap_score = 0
    alignments = aligner.align(read_seq.upper(), identifier_seq.upper())
    align = alignments[0]

    return align


def get_consensus(sequence_1, sequence_2):
    consensus_sequence = ''
    for index, nucleotideo in enumerate(sequence_1):
        if nucleotideo == sequence_2[index]:
            consensus_sequence += nucleotideo
        elif nucleotideo == '-' and sequence_2[index] != '-':
            consensus_sequence += sequence_2[index]
        elif nucleotideo == 'N' and sequence_2[index] != '-':
            consensus_sequence += sequence_2[index]
        elif sequence_2[index] == '-' and nucleotideo != '-':
            consensus_sequence += nucleotideo
        elif sequence_2[index] == 'N' and nucleotideo != '-':
            consensus_sequence += nucleotideo
        elif nucleotideo != sequence_2[index]:
            consensus_sequence += sequence_2[index]

    return consensus_sequence


def consensus_pairwise_alignment(ab1_folder):
    ''''A pairwise alignment is applied in a fasta file with two sequencing reads'''
    pair_fasta_folder = os.path.join(ab1_folder, 'pair')
    '''Check if output directory exists'''
    output_folder = os.path.join(ab1_folder, 'alignment_reads')
    if not os.path.exists(output_folder):
        os.mkdir(output_folder)
    output_folder_fasta = os.path.join(ab1_folder, 'fasta_consensus')
    if not os.path.exists(output_folder_fasta):
        os.mkdir(output_folder_fasta)

    for pair_fasta_file in sorted(os.listdir(pair_fasta_folder)):
        file_path = os.path.join(pair_fasta_folder, pair_fasta_file)
        mafft_cline = MafftCommandline(adjustdirection=True, input=file_path, lep=-0.5)
        stdout, stderr = mafft_cline()
        align = AlignIO.read(StringIO(stdout.upper()), "fasta")

        consensus = get_consensus(align[0].seq, align[1].seq)
        consensus_rec = SeqRecord(Seq(consensus), id=pair_fasta_file[:-3] + '_consensus')

        '''Save alignment_consensus in file'''
        with open(os.path.join(output_folder, pair_fasta_file[:-3]+'.align'), 'w') as handle_align:
            AlignIO.write(align, handle_align, 'fasta')

        with open(os.path.join(output_folder_fasta, pair_fasta_file[:-3] + '_consensus.fa'), 'w') as handle_consensus:
            SeqIO.write(consensus_rec, handle_consensus, 'fasta')


def concatenate_consensus_db_identifier(ab1_folder):
    fasta_consensus_folder = os.path.join(ab1_folder, 'fasta_consensus')
    for consensus_fasta_file in sorted(os.listdir(fasta_consensus_folder)):
        if os.path.isfile(os.path.join(fasta_consensus_folder, consensus_fasta_file)):
            consensus_fasta_read = SeqIO.read(open(os.path.join(fasta_consensus_folder, consensus_fasta_file)), "fasta")
            for db in DATABASES_SEQ:
                output_folder = os.path.join(fasta_consensus_folder, db.name)
                if not os.path.exists(output_folder):
                    os.mkdir(output_folder)

                with open(os.path.join(output_folder, consensus_fasta_read.name + '_' + db.name + '.fa'), 'w') \
                        as handle:
                    SeqIO.write([consensus_fasta_read, db], handle, 'fasta')


def identify_database(ab1_folder):
    start_time = round(time.time()*1000)
    count_num_alignments = 0
    fasta_consensus_folder = os.path.join(ab1_folder, 'fasta_consensus')

    for db in DATABASES_SEQ:
        fasta_consensus_db_folder = os.path.join(ab1_folder, 'fasta_consensus', db.name)

        for consensus_fasta_db_file in sorted(os.listdir(fasta_consensus_db_folder)):
            read_label = str(consensus_fasta_db_file).replace(str('_' + db.name), '')
            consensus_fasta_read = SeqIO.read(open(os.path.join(fasta_consensus_folder, read_label)), "fasta")

            if read_label in DICT_READS:
                read_temp = DICT_READS[read_label]
            else:
                read_temp = sequencing.Read(read_label, consensus_fasta_read.seq)

            if read_temp.db_score < 100:
                count_num_alignments +=1
                align1 = score_alignment(consensus_fasta_read.seq,db.seq)
                score1 = (align1.score/len(str(db.seq)))*100
                # print(str(score1) +'% '+ str(consensus_fasta_read.name) + ' + ' + str(db.name) + ': ' + str(align1))
                if score1 < 100:
                    align2 = score_alignment(consensus_fasta_read.reverse_complement().seq, db.seq)
                    score2 = (align2.score / len(str(db.seq))) * 100
                    # print(str(score2) +'% '+ str(consensus_fasta_read.name) + ' + ' + str(db.name) + ': ' + str(align2))
                score = max(score1, score2)

                if read_temp.db_score < score:
                    read_temp.db_score = round(score,0)
                    read_temp.db_name = db.name
                    read_temp.db_sequence = str(db.seq).upper()
                    DICT_READS[read_label] = read_temp

                delta = (round(time.time()*1000) - start_time)
                print('\r' + 'Total alignments ' + str(count_num_alignments) + ' performed in ' + str(delta)
                      + 'ms. Working on database: ' + str(db.name) + ' speed: '
                      + str(round(count_num_alignments * 1000 / delta, 0)) + ' alignments/sec', end='')

    return DICT_READS


def print_results_database(ab1_folder):
    '''Functions receive a dictionary with Read object'''
    '''Every key in dictionary is associated to consensus read sequence'''
    data = []
    output_folder = os.path.join(ab1_folder, 'result')
    if not os.path.exists(output_folder):
        os.mkdir(output_folder)

    for read in DICT_READS:
        seq = DICT_READS[read]
        data.append([seq.name, seq.sequence, seq.db_name, seq.db_sequence, seq.db_score])

    df = pandas.DataFrame(data, columns=['Read', 'R_Sequence', 'Database', 'DB_Sequence', 'High_score'])
    df = df.sort_values(by=['Database'])
    df.to_excel(os.path.join(output_folder,'db_result.xlsx'), index=False)

    print('\nDB result file placed at: %s' % output_folder)


def load_read_from_db_result(ab1_folder):
    error = ''
    input_file = os.path.join(ab1_folder, 'result', 'db_result.xlsx')
    if not os.path.exists(input_file):
        error = 'File with Database results not found!'
        return error, None
    else:
        df = pandas.read_excel(input_file)

        for idx, row in df.iterrows():
            read_temp = sequencing.Read(row[0], row[1])
            read_temp.db_name = row[2]
            read_temp.db_sequence = row[3]
            read_temp.db_score = row[4]
            DICT_READS[row[0]] = read_temp

    return error, DICT_READS


def load_gene_db(db_label):
    activation_file_path = '/home/flavia/Downloads/G_MART0122/gRNA_barcode information-20220120T180349Z-001/gRNA_barcode information/Activation gRNA database.csv'
    deletion_file_path = '/home/flavia/Downloads/G_MART0122/gRNA_barcode information-20220120T180349Z-001/gRNA_barcode information/Deletion gRNA datatbase.csv'
    interference_file_path = '/home/flavia/Downloads/G_MART0122/gRNA_barcode information-20220120T180349Z-001/gRNA_barcode information/Interference gRNA database.csv'

    if db_label == 'Activation':
        df = pandas.read_csv(activation_file_path)
        return df
    elif db_label == 'Deletion':
        df = pandas.read_csv(deletion_file_path)
        return df
    else:
        df = pandas.read_csv(interference_file_path)
        return df

def mount_virtual_path(ab1_folder):
    # path = '/mnt/ram'
    # cmd = 'mount -t ramfs -o size=20m ramfs /mnt/ram'
    # os.system(cmd)

    path = '/dev/shm'
    if not os.path.exists(path):
        os.mkdir(path)

    return path


def print_results_gene(ab1_folder):
    '''Functions receive a dictionary with Read object'''
    '''Every key in dictionary is associated to consensus read sequence'''
    data = []
    output_folder = os.path.join(ab1_folder, 'result')
    if not os.path.exists(output_folder):
        os.mkdir(output_folder)

    for read in DICT_READS:
        seq = DICT_READS[read]
        data.append([seq.name, seq.sequence, seq.db_name, seq.db_sequence, seq.db_score, seq.gene_name, seq.gene_number,
                     seq.gene_score, seq.gene_sequence])

    df = pandas.DataFrame(data, columns=['Read', 'R_Sequence', 'Database', 'DB_Sequence', 'DB_High_score', 'Gene',
                                         'Number', 'Gene_High_Score', 'Gene_Sequence'])
    df = df.sort_values(by=['Database'])
    df.to_excel(os.path.join(output_folder,'gene_result.xlsx'), index=False)


def identify_gene(ab1_folder, temp_path, start_time):
    output_filepath = os.path.join(temp_path, 'temp.fa')
    open(output_filepath, 'w').close()
    count_num_alignments = 0

    with open(output_filepath, 'r+') as handle:
        for db in DATABASES_SEQ:
            df_genes = load_gene_db(db.name)
            for read in DICT_READS:
                seq_count = 0
                '''Create a seq object for a seq'''
                seq = DICT_READS[read]
                consensus_fasta_read = SeqRecord(Seq(str(seq.sequence).upper()),
                                     name=seq.name,
                                     description="")

                if (str(seq.db_name) == str(db.name)) and seq.gene_score < 100:
                    for idx, row in df_genes.iterrows():
                        seq_count += 1
                        '''Create a seq object for a gene'''
                        gene_rec = SeqRecord(Seq(str(row['Sequence']).upper()),
                                                id=str(row['Number']),
                                                name=str(row['Name']),
                                                description="")

                        '''Write output file'''
                        handle.seek(0)
                        handle.truncate(0)
                        SeqIO.write([consensus_fasta_read, gene_rec], handle, 'fasta')
                        handle.flush()

                        mafft_cline = MafftCommandline(adjustdirection=True, input=output_filepath, lep=-0.5)
                        stdout, stderr = mafft_cline()
                        align_file = AlignIO.read(StringIO(stdout.upper()), "fasta")
                        align = score_alignment(align_file[0].seq, align_file[1].seq)
                        score = round((align.score/len(str(align_file[1].seq).replace('-','')))*100,0)
                        count_num_alignments +=1

                        if seq.gene_score < score:
                            seq.gene_score = score
                            seq.gene_name = str(gene_rec.name)
                            seq.gene_number = str(gene_rec.id)
                            seq.gene_sequence = str(gene_rec.seq).upper()
                            # print('\n', seq.name, seq.db_name, seq.db_score, seq.gene_number, seq.gene_name, seq.gene_score, seq.gene_sequence)
                        delta = int(time.time() - start_time)
                        print('\r' + 'Total alignments ' + str(count_num_alignments) + ' performed. '
                              + 'Working on database: ' + str(db.name) + ' #seq: ' + str(count_num_alignments)
                              + str(count_num_alignments*60/delta) + ' alignments/min', end='')
                print_results_gene(ab1_folder)

    return DICT_READS


def process_job_biopython(job_description):
    job_id, read, db, output_filepath, count_num_alignments = job_description
    df_genes = load_gene_db(db.name)
    '''Create a seq object for a consensus read'''
    consensus_fasta_read = SeqRecord(Seq(str(read.sequence).upper()),
                                     name=read.name,
                                     description="")

    for idx, row in df_genes.iterrows():
        '''Create a seq object for a gene'''
        gene_rec = SeqRecord(Seq(str(row['Sequence']).upper()),
                             id=str(row['Number']),
                             name=str(row['Name']),
                             description="")

        if read.gene_score < 100:
            count_num_alignments += 1
            '''Biopython alignment'''
            align1 = score_alignment(consensus_fasta_read.seq, gene_rec.seq)
            score1 = round((align1.score / len(str(gene_rec.seq))) * 100, 0)

            if score1 < 100:
                align2 = score_alignment(consensus_fasta_read.reverse_complement().seq, gene_rec.seq)
                score2 = round((align2.score / len(str(gene_rec.seq))) * 100, 0)

            score = max(score1, score2)

            if read.gene_score < score:
                read.gene_score = score
                read.gene_name = str(gene_rec.name)
                read.gene_number = str(gene_rec.id)
                read.gene_sequence = str(gene_rec.seq).upper()
                DICT_READS[read.name] = read

    job_description[4] = count_num_alignments
    return job_description


def process_job_mafft(job_description):
    count_num_alignments = 0
    job_id, read, db, output_filepath = job_description
    df_genes = load_gene_db(db.name)
    '''Create a seq object for a consensus read'''
    consensus_fasta_read = SeqRecord(Seq(str(read.sequence).upper()),
                         name=read.name,
                         description="")
    with open(output_filepath, 'w+') as handle:
        for idx, row in df_genes.iterrows():
            '''Create a seq object for a gene'''
            gene_rec = SeqRecord(Seq(str(row['Sequence']).upper()),
                         id=str(row['Number']),
                         name=str(row['Name']),
                         description="")

            '''Write output file'''
            handle.seek(0)
            handle.truncate(0)
            SeqIO.write([consensus_fasta_read, gene_rec], handle, 'fasta')
            handle.flush()

            if read.gene_score < score:
                read.gene_score = score
                read.gene_name = str(gene_rec.name)
                read.gene_number = str(gene_rec.id)
                read.gene_sequence = str(gene_rec.seq).upper()
                DICT_READS[read.name] = read

            mafft_cline = MafftCommandline(adjustdirection=True, input=output_filepath, lep=-0.5)
            stdout, stderr = mafft_cline()
            align_file = AlignIO.read(StringIO(stdout.upper()), "fasta")
            align = score_alignment(align_file[0].seq, align_file[1].seq)
            score = round((align.score / len(str(align_file[1].seq).replace('-', ''))) * 100, 0)

    os.remove(output_filepath)

    return count_num_alignments


def job_descriptions_generator():
    job_id = 0
    for db in DATABASES_SEQ:
        for read_name in DICT_READS:
            read = DICT_READS[read_name]
            output_filepath = '/dev/shm/gene/temp' + str(job_id)
            if (str(read.db_name) == str(db.name)) and read.gene_score < 100:
                job_description = [job_id, read, db, output_filepath, 0]
                job_id+=1
                yield job_description


def identify_gene_parallel():
    start_time = round(time.time()*1000)
    count_total_alignments = 0

    job_descriptions = job_descriptions_generator()
    with concurrent.futures.ProcessPoolExecutor() as executor:
        jobs = [executor.submit(process_job_biopython, job_description) for job_description in job_descriptions]

        for f in concurrent.futures.as_completed(jobs):
            job_description = f.result()
            job_id, read, db, output_filepath, alignments_done = job_description
            DICT_READS[read.name] = read
            count_total_alignments+=alignments_done

            delta = (round(time.time() * 1000) - start_time)
            print('\r' + 'Total alignments ' + str(count_total_alignments) + ' performed in ' + str(round(delta/60000,0)) + 'min'
                  + ' Speed: ' + str(round(count_total_alignments * 1000 / delta, 0)) + ' alignments/sec '
                  + str(read.name) + ' ' + str(read.db_name) + ' ' + str(read.db_score) + ' ' + str(read.gene_name)
                  + ' ' + str(read.gene_number) + ' ' + str(read.gene_score), end='' )


# if __name__ == '__main__':
#     '''Folder with sequencing reads'''
#     ab1_folder = '/home/flavia/Downloads/G_MART0122/Sanger_seq_test data'
#     main.mainloop()
    #
    # '''Parse files'''
    # convert_ab12fasta(ab1_folder)
    # convert_ab12fastq(ab1_folder)
    # concatenate_fwd_rev_fasta(ab1_folder)
    # consensus_pairwise_alignment(ab1_folder)
    # concatenate_consensus_db_identifier(ab1_folder)

    # '''Identify database labels'''
    # identify_database(ab1_folder)
    # print_results_database(ab1_folder)

    # '''Using consensus files add each one of the database identifier'''
    # error, DICT_READS = load_read_from_db_result(ab1_folder)
    # if DICT_READS is not None:
    #     temp_path = mount_virtual_path(ab1_folder)
    #     identify_gene_parallel()
    #     print_results_gene(ab1_folder)
    #
    # else:
    #     print(error)