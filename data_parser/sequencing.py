'''Class to represent the information that the read holds'''
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

db_activation_rec = SeqRecord(Seq("AATTTCTACTAAGTGTAGAT"),
                                    id="Activation_identifier",
                                    name="Activation",
                                    description="")

db_deletion_rec = SeqRecord(Seq("GTTTTAGTACTCTGTAATTT"),
                                  id="Deletion_identifier",
                                  name="Deletion",
                                  description="")

db_interference_rec = SeqRecord(Seq("GTTTTAGAGCTAGAAATAGC"),
                                    id="Interference_identifier",
                                    name="Interference",
                                    description="")

databases = [db_activation_rec, db_deletion_rec, db_interference_rec,]


class Read:

    def __init__(self, name, sequence):
        self.name = name
        self.sequence = sequence
        self.plate = ''
        self.well = ''
        self.db_name = ''
        self.db_score = 0
        self.db_sequence = ''
        self.db_right_end = 0
        self.db_left_end = 0
        self.gene_name = ''
        self.gene_score = 0
        self.gene_number = ''
        self.gene_sequence = ''

    def add_db_info(self, db_name, db_score, db_sequence):
        self.db_name = db_name
        self.db_score = db_score
        self.db_sequence = db_sequence

    def add_gene_info(self, gene_name, gene_number, gene_score, gene_sequence):
        self.gene_name = gene_name
        self.gene_number = gene_number
        self.gene_score = gene_score
        self.gene_sequence = gene_sequence


class DB_gene:
    def __init__(self, db_name, gene_name, gene_number):
        self.db_name = db_name
        self.gene_name = gene_name
        self.gene_number = gene_number
