import os, pandas, numpy, time, re, concurrent.futures, gc, queue
import multiprocessing as mp
import linecache, tracemalloc
from datetime import datetime
from data_parser import sequencing, parser
from io import StringIO
from Bio import SeqIO, AlignIO
from Bio.Align.Applications import MafftCommandline
from Bio.Align import PairwiseAligner
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

DATABASES_SEQ = sequencing.databases
q = queue.Queue()


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

    for pair_fasta_file in parser.natural_sort(os.listdir(pair_fasta_folder)):
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


def identify_database(ab1_folder):
    start_time = round(time.time()*1000)
    count_num_alignments = 0
    fasta_consensus_folder = os.path.join(ab1_folder, 'fasta_consensus')

    for db in DATABASES_SEQ:
        fasta_consensus_db_folder = os.path.join(ab1_folder, 'fasta_consensus', db.name)

        for consensus_fasta_db_file in parser.natural_sort(os.listdir(fasta_consensus_db_folder)):
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
                    align2 = score_alignment(consensus_fasta_read.seq, db.reverse_complement().seq)
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
                print('Total alignments ' + str(count_num_alignments) + ' performed in ' + str(delta)
                      + 'ms. Working on database: ' + str(db.name) + ' speed: '
                      + str(round(count_num_alignments * 1000 / delta, 0)) + ' alignments/sec')

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
            plate_well = str(row[0]).upper().split('PLATE')[1]
            read_temp.plate = plate_well.split('_')[1]
            read_temp.well = plate_well.split('_')[2]
            DICT_READS[row[0]] = read_temp

    return error, DICT_READS


def mount_virtual_path(ab1_folder):
    # path = '/mnt/ram'
    # cmd = 'mount -t ramfs -o size=20m ramfs /mnt/ram'
    # os.system(cmd)
    path = '/dev/shm'
    if not os.path.exists(path):
        os.mkdir(path)

    return path


def identify_gene(ab1_folder, temp_path, start_time):
    output_filepath = os.path.join(temp_path, 'temp.fa')
    open(output_filepath, 'w').close()
    count_num_alignments = 0

    with open(output_filepath, 'r+') as handle:
        for db in DATABASES_SEQ:
            df_genes = parser.load_gene_db(db.name)
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
                            if str(seq.name).__contains__('MagicSanger286'):
                                print(score, seq.gene_score, seq.gene_name, gene_rec.name)
                            seq.gene_score = score
                            seq.gene_name = str(gene_rec.name)
                            seq.gene_number = str(gene_rec.id)
                            seq.gene_sequence = str(gene_rec.seq).upper()

                        elif seq.gene_score == score:
                            if str(seq.name).__contains__('MagicSanger286'):
                                print(score, seq.gene_score, seq.gene_name, gene_rec.name)
                            seq.gene_name = str(seq.gene_name) + ', ' + str(gene_rec.name)
                            seq.gene_number = str(seq.gene_number) + ', ' + str(gene_rec.id)
                            seq.gene_sequence = str(seq.gene_sequence) + ', ' + str(gene_rec.seq).upper()

                            if str(seq.name).__contains__('MagicSanger286'):
                                print(seq.name, seq.db_name, seq.db_score, seq.gene_number, seq.gene_name, seq.gene_score)
                        delta = int(time.time() - start_time)
                        print('Total alignments ' + str(count_num_alignments) + ' performed. '
                              + 'Working on database: ' + str(db.name) + ' #seq: ' + str(count_num_alignments)
                              + str(count_num_alignments*60/delta) + ' alignments/min')
                print_results_gene(ab1_folder)

    return DICT_READS


def process_job_mafft(job_description):
    count_num_alignments = 0
    job_id, read, db, output_filepath = job_description
    df_genes = parser.load_gene_db(db.name)
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


def print_results_from_queue(ab1_folder, read):
    data = []
    date_time = datetime.now().strftime("%Y-%m-%d")
    output_folder = os.path.join(ab1_folder, 'result')
    filepath = os.path.join(output_folder, 'gene_result_3_' + date_time + '.csv')
    if not os.path.exists(output_folder):
        os.mkdir(output_folder)

    data.append([read.name, read.sequence, read.db_name, read.db_sequence, read.db_score, read.gene_name, read.gene_number, read.gene_score, read.gene_sequence])

    df = pandas.DataFrame(data, columns=['Read', 'R_Sequence', 'Database', 'DB_Sequence', 'DB_High_score', 'Gene', 'Number', 'Gene_High_Score', 'Gene_Sequence'])

    if os.path.isfile(filepath):
        with open(filepath, 'a') as f:
            df.to_csv(f, header=False, index=False)
    else:
        df.to_csv(filepath, index=False)


def print_results_gene(ab1_folder, DICT_READS):
    '''Functions receive a dictionary with Read object'''
    '''Every key in dictionary is associated to consensus read sequence'''
    data = []
    date_time = datetime.now().strftime("%Y-%m-%d_%H.%M.%S")
    output_folder = os.path.join(ab1_folder, 'result')
    if not os.path.exists(output_folder):
        os.mkdir(output_folder)

    for read in DICT_READS:
        seq = DICT_READS[read]
        # data.append([seq.name, seq.plate, seq.well, seq.sequence, seq.db_name, seq.db_sequence, seq.db_score, seq.gene_name, seq.gene_number, seq.gene_score, seq.gene_sequence])
        data.append([seq.name, seq.sequence, seq.db_name, seq.db_sequence, seq.db_score, seq.gene_name, seq.gene_number, seq.gene_score, seq.gene_sequence])

    df = pandas.DataFrame(data, columns=['Read', 'R_Sequence', 'Database', 'DB_Sequence', 'DB_High_score', 'Gene', 'Number', 'Gene_High_Score', 'Gene_Sequence'])
    # df = pandas.DataFrame(data, columns=['Read', 'Plate', 'Well', 'R_Sequence', 'Database', 'DB_Sequence', 'DB_High_score', 'Gene', 'Number', 'Gene_High_Score', 'Gene_Sequence'])

    df.to_csv(os.path.join(output_folder, 'gene_result_1_'+date_time+'.csv'), index=False)
    DICT_READS.clear()
    print('\nDB result file placed at: %s' % output_folder)


def display_top(snapshot, key_type='lineno', limit=3):
    snapshot = snapshot.filter_traces((
        tracemalloc.Filter(False, "<frozen importlib._bootstrap>"),
        tracemalloc.Filter(False, "<unknown>"),
    ))
    top_stats = snapshot.statistics(key_type)

    print("Top %s lines" % limit)
    for index, stat in enumerate(top_stats[:limit], 1):
        frame = stat.traceback[0]
        # replace "/path/to/module/file.py" with "module/file.py"
        filename = os.sep.join(frame.filename.split(os.sep)[-2:])
        print("#%s: %s:%s: %.1f KiB"
              % (index, filename, frame.lineno, stat.size / 1024))
        line = linecache.getline(frame.filename, frame.lineno).strip()
        if line:
            print('    %s' % line)

    other = top_stats[limit:]
    if other:
        size = sum(stat.size for stat in other)
        print("%s other: %.1f KiB" % (len(other), size / 1024))
    total = sum(stat.size for stat in top_stats)
    print("Total allocated size: %.1f KiB" % (total / 1024))