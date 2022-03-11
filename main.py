import sys

from sanger import sequencing, sanger
from tkinter import *
from tkinter import filedialog

DATABASES_SEQ = sequencing.databases
DICT_READS = {}

# def select_path():
#     folder_selected = filedialog.askdirectory(initialdir = "/",title = "Select file",)
#
# #Creating a object window
# main = Tk()
# main.geometry('500x500')
#
# frame1 = Frame(main)
# frame1.pack(expand=True, fill=BOTH)
#
# button1 = Button(frame1, text='Select file', command=select_path)
# button1.grid(row=2, column=0)
#
# button2 = Button(frame1, text='Select file')
# button2.grid(row=4, column=0)

def main(argsv):
    '''Folder with sequencing reads'''
    ab1_folder = '/home/flavia/Downloads/G_MART0122/Sanger_seq_test data'
    quality_threshold = 20
    # main.mainloop()

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
        sanger.identify_gene_parallel()
        sanger.print_results_gene(ab1_folder)

    else:
        print(error)

if __name__ == '__main__':
    main(sys.argv)
