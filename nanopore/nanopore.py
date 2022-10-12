import os, sys, fnmatch
from sanger import sanger
from data_parser import sequencing, parser
from ngs import ngs
from Bio import SeqIO, AlignIO


def template_sorted_alignment(templates_folder, template_sorted, assemblies_folder, assemblies_by_plate):

    folder_count = 0
    for folder in assemblies_by_plate:
        count_index = 0
        folder_count = folder_count + 1
        for assembly in parser.natural_sort(os.listdir(folder)):
            if os.path.isfile(os.path.join(folder, assembly)):
                consensus_fasta_read = SeqIO.read(open(os.path.join(folder, assembly)), "fasta")

                template_name = 'P0' + str(folder_count) + template_sorted[count_index]
                if template_name == 'P01W92':
                    template_name = 'P01W91'
                    count_index = count_index + 1

                if template_name == 'P03W93':
                    template_name = 'P03W92'
                    count_index = count_index + 1

                if template_name == 'P03W27':
                    template_name = 'P03W26'
                    count_index = count_index + 1

                if template_name == 'P04W33':
                    template_name = 'P04W32'
                    count_index = count_index + 1

                for template_file in os.listdir(templates_folder):
                    if fnmatch.fnmatch(template_file, template_name+'_*'):
                        template_fasta_read = SeqIO.read(open(os.path.join(templates_folder, template_file)),
                                                         "fasta")

                        align = ngs.score_alignment(template_fasta_read.seq, consensus_fasta_read.seq)
                        align_rev = ngs.score_alignment(template_fasta_read.seq, consensus_fasta_read.reverse_complement().seq)

                        if align_rev.score > align.score:
                            align = align_rev
                        score = (align.score / (len(template_fasta_read.seq)-32))* 100
                        print(template_file[:-3] +'\t'+ os.path.basename(folder) +'-'+ assembly[:-6] +'\t' + str(round(align.score)) +'\t' + str(round(score))+'%')
                        if score < 100:
                            print(align)

                count_index = count_index + 1


def template_alignment(templates_folder, template_files, assemblies_folder, assemblies_by_plate):
    count_index = 0
    for folder in assemblies_by_plate:
        for assembly in parser.natural_sort(os.listdir(folder)):
            if os.path.isfile(os.path.join(folder, assembly)):
                # print(template, assembly)
                consensus_fasta_read = SeqIO.read(open(os.path.join(folder, assembly)), "fasta")

                # print(os.path.join(templates_folder, template))
                template_fasta_read = SeqIO.read(open(os.path.join(templates_folder, template_files[count_index])),
                                                 "fasta")
                align = ngs.local_alignment(template_fasta_read.seq, consensus_fasta_read.seq)

                align1, align2, score, begin, end = align[0]
                print(template_files[count_index], os.path.basename(folder), assembly, round(score))
                header = str(template_files[count_index]) + ' ' + str(os.path.basename(folder)) + ' ' + str(
                    assembly) + ' ' + str(round(score))
                a, b = align[0][:2]

                n = 200
                # Chunk each sequence according to size
                chunked = [[a[i: i + n], b[i: i + n]] for i in range(0, len(a), n)]
                # Asses matches for each chunk
                matches = [''.join(['|' if chunk[0][i] == chunk[1][i] else ' ' for i in range(len(chunk[0]))]) for chunk
                           in chunked[:]]
                # Print
                with open(os.path.join('/home/flavia/Documents/Concordia/Project/Nanopore/alignments', 'all.align'),
                          'a') as handle_align:
                    handle_align.write(header)

                    for c, m, i in zip(chunked, matches, range(len(chunked))):
                        txt = '%s to %s' % (i * n + 1, i * n + n)
                        print(txt)
                        print(c[0])
                        print(m)
                        print(c[1])
                        print()

                        handle_align.write(txt)
                        handle_align.write('\n')
                        handle_align.write(c[0])
                        handle_align.write('\n')
                        handle_align.write(m)
                        handle_align.write('\n')
                        handle_align.write(c[1])
                        handle_align.write('\n')

                handle_align.close()
                count_index = count_index + 1

                # align1, align2, score, begin, end = align
                # score1 = round((score / len(str(consensus_fasta_read.seq))) * 100, 0)
                # print(template_files[count_index], os.path.basename(folder), assembly, len(str(consensus_fasta_read.seq)), round(score), score1)
                # print(align)
                # count_index = count_index + 1
                # with open(os.path.join('/home/flavia/Documents/Concordia/Project/Nanopore/alignments', 'all.align'), 'w+') as handle_align:
                #     AlignIO.write(align, handle_align, "fasta")
                # handle_align.close()

        # for template in template_files:
    #     max_score = 0
    #     list_assemblies = []
    #     for folder in assemblies_by_plate:
    #         for assembly in parser.natural_sort(os.listdir(folder)):
    #             if os.path.isfile(os.path.join(folder, assembly)):
    #                 # print(template, assembly)
    #                 consensus_fasta_read = SeqIO.read(open(os.path.join(folder, assembly)), "fasta")
    #                 # print(os.path.join(templates_folder, template))
    #                 template_fasta_read = SeqIO.read(open(os.path.join(templates_folder, template)), "fasta")
    #                 align1 = ngs.score_alignment(template_fasta_read.seq, consensus_fasta_read.seq)
    #                 score1 = round((align1.score / len(str(template_fasta_read.seq))) * 100, 0)
    #                 if align1.score > max_score:
    #                     max_score = align1.score
    #                     list_assemblies = []
    #                     list_assemblies.append([os.path.basename(folder), assembly])
    #                 elif align1.score == max_score:
    #                     list_assemblies.append([os.path.basename(folder), assembly])
    #                 print(template, os.path.basename(folder), assembly, align1.score)
    #                 print(align1.score)
    #     print(template, list_assemblies, max_score)

def my_format_alignment(result):
    barcode_id, assembly_folder, assembly, align, align1, align2, score, perc_score, begin, end, direction = result
    s = [] #All alignment and result
    r = [] # resume result
    s.append("%s\t" % barcode_id)
    r.append("%s\t" % barcode_id)
    s.append("%s\t" % assembly_folder)
    r.append("%s\t" % assembly_folder)
    s.append("%s\n" % assembly)
    r.append("%s\t" % assembly)
    gap=0
    mismatch=0
    match=0
    for a, b in zip(align1[begin:end], align2[begin:end]):
        if a == b:
            # s.append("|")  # match
            match+=1
        elif a == "-" or b == "-":
            # s.append(" ")  # gap
            gap+=1
        else:
            # s.append(".")  # mismatch
            mismatch+=1

    # s.append("\n")

    a, b = align[0][:2]
    n = 200
    # Chunk each sequence according to size
    chunked = [[a[i: i + n], b[i: i + n]] for i in range(0, len(a), n)]
    # Asses matches for each chunk
    matches = [''.join(['|' if chunk[0][i] == chunk[1][i] else ' ' for i in range(len(chunk[0]))]) for
               chunk
               in chunked[:]]
    # Print
    for c, m, i in zip(chunked, matches, range(len(chunked))):
        txt = '%s to %s' % (i * n + 1, i * n + n)
        s.append("%s\n" % txt)
        s.append("%s\n" % c[0])
        s.append("%s\n" % m)
        s.append("%s\n" % c[1])
        s.append("\n")

    s.append("  Score=%g\n" % score)
    r.append("  %g\t" % score)
    s.append("  Perc_Score=%s\n" % perc_score)
    r.append("  %s\t" % perc_score)
    s.append("  Match=%g\n" % match)
    r.append("  %g\t" % match)
    s.append("  Gap=%g\n" % gap)
    r.append("  %g\t" % gap)
    s.append("  Mismatch=%g\n" % mismatch)
    r.append("  %g\t" % mismatch)
    s.append("  Direction=%s\n" % direction)
    r.append("  %s\n" % direction)
    s.append("\n")

    write_result = ''.join(s)
    resume_result = ''.join(r)
    return write_result, resume_result


def barcode_alignment(barcode5_fasta, barcode3_fasta, assemblies_folder, assemblies_by_plate):
    count_index = 0

    dict = {}
    rev = False

    for folder in assemblies_by_plate:
        for assembly in parser.natural_sort(os.listdir(folder)):
            max_result = 0
            result = []
            barcode_fasta_read = SeqIO.parse(barcode3_fasta, "fasta")
            if os.path.isfile(os.path.join(folder, assembly)):
                for barcode in barcode_fasta_read:
                    consensus_fasta_read = SeqIO.read(open(os.path.join(folder, assembly)), "fasta")
                    align = ngs.local_alignment(consensus_fasta_read.seq, barcode.seq)
                    align1, align2, score, begin, end = align[0]
                    perc_score = score/len(barcode.seq)*100

                    if perc_score > max_result:
                        max_result = perc_score
                        result = barcode.id, os.path.basename(folder), assembly, align, align1, align2, score, "{:.2f}".format(perc_score), begin, end, 'F'

                    # if perc_score < 100:
                    align = ngs.local_alignment(consensus_fasta_read.reverse_complement().seq, barcode.seq)
                    align1, align2, score, begin, end = align[0]
                    rev_perc_score = score/len(barcode.seq)*100

                    if rev_perc_score > max_result:
                        max_result = rev_perc_score
                        result = barcode.id, os.path.basename(folder), assembly, align, align1, align2, score, "{:.2f}".format(rev_perc_score), begin, end, 'R'

                write_result, resume_result = my_format_alignment(result)
                print(resume_result)
                with open(os.path.join('/home/flavia/Documents/Concordia/Project/Nanopore/alignments', 'barcode3.align'), 'a') as handle_align:
                    handle_align.write(write_result)
                    handle_align.close()
                with open(os.path.join('/home/flavia/Documents/Concordia/Project/Nanopore/alignments', 'barcode3-resume.align'), 'a') as handle_align:
                    handle_align.write(resume_result)
                    handle_align.close()