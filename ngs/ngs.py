import os, re, concurrent.futures
from Bio import SeqIO
from concurrent.futures import ProcessPoolExecutor
MAX_WORKERS = 10


def natural_sort(l):
    convert = lambda text: int(text) if text.isdigit() else text.lower()
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ]
    return sorted(l, key = alphanum_key)


def fastq2fasta():
    path = '/Users/flavia/Documents/NGS_project/2022_0519_NGS_results/2-Merge-QC'

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