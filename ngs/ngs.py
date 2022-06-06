import os, pandas, gzip, numpy, time, re, concurrent.futures
from concurrent.futures import ProcessPoolExecutor
from Bio import SeqIO, AlignIO
from Bio.Align.Applications import MafftCommandline
from Bio.Align import PairwiseAligner
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from sanger import sequencing

MAX_WORKERS = 10


DATABASES_SEQ = sequencing.databases
DICT_READS = {}


def natural_sort(l):
    convert = lambda text: int(text) if text.isdigit() else text.lower()
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ]
    return sorted(l, key = alphanum_key)


def process_task(task_description):
    task_id, record, output_path = task_description
    with open(os.path.join(output_path, str(record.id) + '.fasta'), "w") as output_handle:
        SeqIO.write(record, output_handle, "fasta")
        return f'Finished file {task_id}'


def create_task(path):
    task_id = 1
    for file in natural_sort(os.listdir(path)):
        if file.__contains__('.assembled.'):
            output_path = os.path.join(path, file.split('---')[0])
            if not os.path.exists(output_path):
                os.mkdir(output_path)

            with open(os.path.join(path, file)) as handle:
                for record in SeqIO.parse(handle, "fastq"):
                    task_description = [task_id, record, output_path]
                    task_id+=1
                    yield task_description


# def fastq2fasta():
#     path = '/Users/flavia/Documents/NGS_project/2022_0519_NGS_results/2-Merge-QC'
#     task_description = create_task(path)
#     with ProcessPoolExecutor(max_workers=MAX_WORKERS) as executor:
#         tasks = [executor.submit(process_task, task) for task in task_description]
#
#         for f in concurrent.futures.as_completed(tasks):
#             result = f.result()
#             print(result)


def fastq2fasta(path):
    # path = '/Users/flavia/Documents/NGS_project/2022_0519_NGS_results/2-Merge-QC'

    for file in natural_sort(os.listdir(path)):
        if file.__contains__('.assembled.'):
            output_path = os.path.join(path, file.split('---')[0])
            if not os.path.exists(output_path):
                os.makedirs(output_path)

            with open(os.path.join(path, file)) as handle:
                count = 1
                for record in SeqIO.parse(handle, "fastq"):
                    with open(os.path.join(output_path, str(record.id)+'.fasta'), "w") as output_handle:
                        SeqIO.write(record, output_handle, "fasta")
                    count+=1
                    print('\r'+str(count))


'''Functions copied and changed form sanger'''


def concatenate_consensus_db_identifier(ngs_folder):
    for ngs_exp_folder in natural_sort(os.listdir(ngs_folder)):
        if os.path.isdir(os.path.join(ngs_folder, ngs_exp_folder)):
            fasta_consensus_folder = os.path.join(ngs_folder, ngs_exp_folder)
            print(f'Working on folder: {ngs_exp_folder}')
            count = 0
            for consensus_fasta_file in natural_sort(os.listdir(fasta_consensus_folder)):
                if os.path.isfile(os.path.join(fasta_consensus_folder, consensus_fasta_file)):
                    count = count + 1
                    print(count)
                    consensus_fasta_read = SeqIO.read(open(os.path.join(fasta_consensus_folder, consensus_fasta_file)), "fasta")
                    for db in DATABASES_SEQ:
                        output_folder = os.path.join(fasta_consensus_folder, db.name)
                        if not os.path.exists(output_folder):
                            os.mkdir(output_folder)

                        with open(os.path.join(output_folder, consensus_fasta_read.name + '_' + db.name + '.fa'), 'w') \
                                as handle:
                            SeqIO.write([consensus_fasta_read, db], handle, 'fasta')


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


def process_job_identify_database(job_description):
    job_id, read_temp, record, ngs_folder, output_path, count_num_alignments = job_description

    for db in DATABASES_SEQ:
        # if read_temp.name in DICT_READS:
        #     read_temp = DICT_READS[read_temp.name]
        # else:
        #     read_temp = sequencing.Read(read_temp.name, read_temp.sequence)

        if read_temp.db_score < 100:
            count_num_alignments += 1
            align1 = score_alignment(record.seq, db.seq)
            score1 = (align1.score / len(str(db.seq))) * 100
            score2 = -1

            if score1 < 100:
                count_num_alignments += 1
                align2 = score_alignment(record.reverse_complement().seq, db.seq)
                score2 = (align2.score / len(str(db.seq))) * 100

            if score1 > score2:
                score = score1
                l_end = min(align1.aligned[0])[0]
                r_end = max(align1.aligned[0])[1]
            else:
                score = score2
                l_end = min(align2.aligned[0])[0]
                r_end = max(align2.aligned[0])[1]
            # print(read_temp.name, db.name, score1, score2, score)
            if read_temp.db_score < score:
                read_temp.db_score = round(score, 0)
                read_temp.db_name = db.name
                read_temp.db_sequence = str(db.seq).upper()
                read_temp.db_right_end = r_end
                read_temp.db_left_end = l_end
                DICT_READS[read_temp.name] = read_temp

    job_description = job_id, read_temp, record, ngs_folder, output_path, count_num_alignments
    return job_description


def job_descriptions_database_generator(ngs_folder, file, output_path):
    job_id = 0
    # with gzip.open(os.path.join(ngs_folder, file), "rt") as handle:
    with open(os.path.join(ngs_folder, file)) as handle:
        for record in SeqIO.parse(handle, "fastq"):
            read_temp = sequencing.Read(record.name, record.seq)
            job_description = [job_id, read_temp, record, ngs_folder, output_path, 0]
            job_id += 1
            yield job_description


def identify_database_parallel(num_threads, ngs_folder, file, output_path):
    start_time = round(time.time()*1000)
    count_total_alignments = 0

    job_descriptions = job_descriptions_database_generator(ngs_folder, file, output_path)
    with concurrent.futures.ProcessPoolExecutor(max_workers=num_threads) as executor:
        jobs = [executor.submit(process_job_identify_database, job_description) for job_description in job_descriptions]

        for f in concurrent.futures.as_completed(jobs):
            job_description = f.result()
            job_id, read_temp, record, ngs_folder, output_path, alignments_done = job_description
            DICT_READS[read_temp.name] = read_temp
            count_total_alignments += alignments_done

            delta = (round(time.time() * 1000) - start_time)
            print('Total alignments ' + str(count_total_alignments) + ' performed in ' + str(round(delta/60000,0)) + 'min'
                  + ' Speed: ' + str(round(count_total_alignments * 1000 / delta, 0)) + ' alignments/sec '
                  + str(read_temp.name) + ' ' + str(read_temp.db_name) + ' ' + str(read_temp.db_score),)


def print_results_database(output_path):
    '''Functions receive a dictionary with Read object'''
    '''Every key in dictionary is associated to consensus read sequence'''
    data = []
    output_folder = os.path.join(output_path, 'result')
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
    df.to_csv(os.path.join(output_folder,'db_result.csv'), index=False)

    print('\nDB result file placed at: %s' % output_folder)


def load_read_from_db_result(database_path):
    error = ''
    input_file = os.path.join(database_path, 'result', 'db_result.csv')
    if not os.path.exists(input_file):
        error = 'File with Database results not found!'
        return error, None
    else:
        df = pandas.read_csv(input_file)

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