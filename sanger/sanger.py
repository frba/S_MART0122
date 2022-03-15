import concurrent.futures
import os, pandas, numpy, time, re
from . import sequencing
from io import StringIO
from Bio import SeqIO, AlignIO
from Bio.Align.Applications import MafftCommandline
from Bio.Align import PairwiseAligner
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from tkinter import *
from tkinter import filedialog

DATABASES_SEQ = sequencing.databases
DICT_READS = {}


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

    '''Check if output directory exists'''
    output_folder = os.path.join(ab1_folder, 'pair')
    if not os.path.exists(output_folder):
        os.mkdir(output_folder)

    fasta_folder = os.path.join(ab1_folder, 'fasta_trimmed')
    total_pair_reads = int(len(os.listdir(fasta_folder)) / 2) + 1

    for i in range(1, total_pair_reads):
        file_label = 'M' + str(i) + '_'
        read_pair = []
        '''Mafft is affected by the order the sequences are added in file'''
        for fasta_file in natural_sort(os.listdir(fasta_folder)):
            if fasta_file.__contains__(file_label):
                read_pair.append(fasta_file)

        fasta_read1 = SeqIO.read(open(os.path.join(ab1_folder, 'fasta_trimmed', read_pair[0])), "fasta")
        fasta_read2 = SeqIO.read(open(os.path.join(ab1_folder, 'fasta_trimmed', read_pair[1])), "fasta")

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


def get_consensus(sequence_1, sequence_2, sequence_1_qual, sequence_2_qual):
    idx_qual_seq_1 = 0
    idx_qual_seq_2 = 0
    consensus_sequence = ''
    for index, nucleotideo in enumerate(sequence_1):
        if nucleotideo == sequence_2[index]:
            consensus_sequence += nucleotideo
            idx_qual_seq_1 += 1
            idx_qual_seq_2 += 1
        elif nucleotideo == '-' and sequence_2[index] != '-':
            consensus_sequence += sequence_2[index]
            idx_qual_seq_2 += 1

        elif sequence_2[index] == '-' and nucleotideo != '-':
            consensus_sequence += nucleotideo
            idx_qual_seq_1 += 1

        elif nucleotideo != sequence_2[index]:
            if sequence_1_qual[idx_qual_seq_1] > sequence_2_qual[idx_qual_seq_2]:
                consensus_sequence += sequence_1[index]
            else:
                consensus_sequence += sequence_2[index]
            idx_qual_seq_1 += 1
            idx_qual_seq_2 += 1

    return consensus_sequence


def consensus_pairwise_alignment(ab1_folder):
    ''''A pairwise alignment is applied in a fasta file with two sequencing reads'''
    pair_fasta_folder = os.path.join(ab1_folder, 'pair')
    trimmed_qual_folder = os.path.join(ab1_folder, 'fasta_trimmed_quality')
    '''Check if output directory exists'''
    output_folder = os.path.join(ab1_folder, 'alignment_reads')
    if not os.path.exists(output_folder):
        os.mkdir(output_folder)
    output_folder_fasta = os.path.join(ab1_folder, 'fasta_consensus')
    if not os.path.exists(output_folder_fasta):
        os.mkdir(output_folder_fasta)

    for pair_fasta_file in natural_sort(os.listdir(pair_fasta_folder)):
        file_path = os.path.join(pair_fasta_folder, pair_fasta_file)
        #lep is used to penalize internal gaps more strongly than terminal gaps
        mafft_cline = MafftCommandline(adjustdirection=True, localpair=True, lep=-0.5, input=file_path)
        stdout, stderr = mafft_cline()

        align = AlignIO.read(StringIO(stdout), "fasta")
        read1_qual_name = align[0].description[+len(align[0].name)+1:]
        read2_qual_name = align[1].description[+len(align[1].name)+1:]
        read1_qual = numpy.loadtxt(os.path.join(ab1_folder, 'fasta_trimmed_quality', read1_qual_name+'.qual'), comments="#", delimiter=",", unpack=False)
        read2_qual = numpy.loadtxt(os.path.join(ab1_folder, 'fasta_trimmed_quality', read2_qual_name+'.qual'), comments="#", delimiter=",", unpack=False)

        '''Array with Phred quality needs to be revert for reverse reads'''
        read1_qual = read1_qual[::-1] if str(align[0].id).startswith('_R') else read1_qual
        read2_qual = read2_qual[::-1] if str(align[1].id).startswith('_R') else read2_qual

        consensus = get_consensus(align[0].seq, align[1].seq, read1_qual, read2_qual)
        consensus_rec = SeqRecord(Seq(consensus), id=pair_fasta_file[:-3] + '_consensus')

        '''Save alignment_consensus in file'''
        with open(os.path.join(output_folder, pair_fasta_file[:-3]+'.align'), 'w') as handle_align:
            AlignIO.write(align, handle_align, 'fasta')

        with open(os.path.join(output_folder_fasta, pair_fasta_file[:-3] + '_consensus.fa'), 'w') as handle_consensus:
            SeqIO.write(consensus_rec, handle_consensus, 'fasta')


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


def identify_database(ab1_folder):
    start_time = round(time.time()*1000)
    count_num_alignments = 0
    fasta_consensus_folder = os.path.join(ab1_folder, 'fasta_consensus')

    for db in DATABASES_SEQ:
        fasta_consensus_db_folder = os.path.join(ab1_folder, 'fasta_consensus', db.name)

        for consensus_fasta_db_file in natural_sort(os.listdir(fasta_consensus_db_folder)):
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
                score2 = -1

                if score1 < 100:
                    count_num_alignments += 1
                    align2 = score_alignment(consensus_fasta_read.reverse_complement().seq, db.seq)
                    score2 = (align2.score / len(str(db.seq))) * 100

                if score1 > score2:
                    score = score1
                    l_end = min(align1.aligned[0])[0]
                    r_end = max(align1.aligned[0])[1]
                else:
                    score = score2
                    l_end = min(align2.aligned[0])[0]
                    r_end = max(align2.aligned[0])[1]

                if read_temp.db_score < score:
                    read_temp.db_score = round(score,0)
                    read_temp.db_name = db.name
                    read_temp.db_sequence = str(db.seq).upper()
                    read_temp.db_right_end = r_end
                    read_temp.db_left_end = l_end
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

        right = seq.db_right_end + 25 if seq.db_right_end + 25 < len(seq.sequence) else len(seq.sequence)-1
        left = seq.db_left_end - 25 if seq.db_left_end - 25 > 0 else 0
        if seq.db_name != 'Deletion':
            seq.sequence = seq.sequence[left:right]
        data.append([seq.name, seq.sequence, seq.db_name, seq.db_sequence, seq.db_score])

    df = pandas.DataFrame(data, columns=['Read', 'R_Sequence', 'Database', 'DB_Sequence', 'High_score'])
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
    # activation_file_path = '/home/flavia/Downloads/G_MART0122/gRNA_barcode information-20220120T180349Z-001/gRNA_barcode information/Activation gRNA database.csv'
    # deletion_file_path = '/home/flavia/Downloads/G_MART0122/gRNA_barcode information-20220120T180349Z-001/gRNA_barcode information/Deletion gRNA datatbase.csv'
    # interference_file_path = '/home/flavia/Downloads/G_MART0122/gRNA_barcode information-20220120T180349Z-001/gRNA_barcode information/Interference gRNA database.csv'
    activation_file_path = '/home/frba/scratch/gRNA_barcode information-20220120T180349Z-001/gRNA_barcode information/Activation gRNA database.csv'
    deletion_file_path = '/home/frba/scratch/gRNA_barcode information-20220120T180349Z-001/gRNA_barcode information/Deletion gRNA datatbase.csv'
    interference_file_path = '/home/frba/scratch/gRNA_barcode information-20220120T180349Z-001/gRNA_barcode information/Interference gRNA database.csv'

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


def identify_gene_parallel(num_threads):
    start_time = round(time.time()*1000)
    count_total_alignments = 0

    job_descriptions = job_descriptions_generator()
    with concurrent.futures.ProcessPoolExecutor(max_workers=num_threads) as executor:
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