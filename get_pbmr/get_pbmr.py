import pysam
import os
import pandas as pd
import matplotlib.pyplot as plt
import argparse

def add_to_dict(MD_string, cigar, counts):
    key_conv = {
        'I': 'Insertion', 'S': 'Soft Clipping', '=': 'Match',
        'X': 'Mismatch'}
    pos = 0
    num = 0
    MD_pos = []
    for i in list(MD_string):
        try:
            tmp = int(i)
            num = num*10 + tmp
        except:
            pos += num
            pos += 1
            if i != '^':
                MD_pos.append(pos)
            num = 0
    pos = 1
    pos_md = 1
    num = 0
    # Read CIGAR ranges and increment the catogories' values
    c_string = cigar
    for i in list(c_string):
        try:
            tmp = int(i)
            num = num*10 + tmp
        except:
            if i == 'M':
                for j in range(pos, pos+num):
                    if j not in counts:
                        counts[j] = {
                            'Match': 0, 'Mismatch': 0,
                            'Insertion': 0, 'Soft Clipping': 0}
                    pos_md += 1
                    if pos_md in MD_pos:
                        counts[j]['Mismatch'] += 1
                    else:
                        counts[j]['Match'] += 1
                pos = pos+num
            elif i in ['X', '=']:
                for j in range(pos, pos+num):
                    if j not in counts:
                        counts[j] = {
                            'Match': 0, 'Mismatch': 0,
                            'Insertion': 0, 'Soft Clipping': 0}
                    pos_md += 1
                    counts[j][key_conv[i]] += 1
                pos = pos+num
            elif i in ['S', 'I']:
                for j in range(pos, pos+num):
                    if j not in counts:
                        counts[j] = {
                            'Match': 0, 'Mismatch': 0,
                            'Insertion': 0, 'Soft Clipping': 0}
                    counts[j][key_conv[i]] += 1
                pos = pos+num
            num = 0

    return counts


def get_bam_mutation_rate_df(filename, counts):
    bamfile = pysam.AlignmentFile(filename, "rb")
    iter = bamfile.fetch(until_eof=True)
    for line_number, x in enumerate(iter):
        print(line_number, end="\r")
        for tag in x.get_tags():
            if tag[0] == 'MD':
                MD_string = tag[1]
        if MD_string and x.cigarstring:
            add_to_dict(MD_string, x.cigarstring, counts)


def get_sam_mutation_rate_df(filename, counts):
    try:
        with open(filename, "r") as input_ob:
            for line_num, file_line in enumerate(read_file(input_ob)):
                print(line_num, end="\r")
                file_line = file_line.strip()
                if file_line[0] == '@':
                    continue
                file_line_arr = file_line.split()
                for c in file_line_arr[11:]:
                    if c.startswith('MD'):
                        MD_string = c.split(':')[2]
                c_string = file_line_arr[5]
                if MD_string and c_string:
                    add_to_dict(MD_string, c_string, counts)

    except (IOError, OSError):
        print("Error opening / processing file")


def gen_mutation_rate_graph(filename=None, prefix=None):

    if not filename and not prefix:
        parser = argparse.ArgumentParser(
            description=('Get Per Base Mutation Rate from sam files. Takes in a'
                         ' SAM file or a BAM file and outputs a PNG file.'),
            add_help=True)

        parser.add_argument(
            '-i', '--input', metavar="<SAM/BAM file>", type=str, action='store',
            dest='input', required=True,
            help=('Input file in sam format.'))

        parser.add_argument(
            '-o', '--output', metavar="PREFIX", type=str, action='store',
            dest='prefix', required=True,
            help=('Prefix for PNG output file.'))

        args = parser.parse_args()
        filename = args.input
        prefix = args.prefix

    file_prefix, file_extension = os.path.splitext(filename)
    counts = {}
    if file_extension == '.bam':
        get_bam_mutation_rate_df(filename, counts)
    elif file_extension == '.sam':
        get_sam_mutation_rate_df(filename, counts)
    else:
        raise Exception("Expected a file ending with .sam or .bam")
    df = pd.DataFrame(counts).T
    df['total'] = df.sum(axis=1)
    for i in ['Insertion', 'Match', 'Mismatch', 'Soft Clipping']:
        df['{0} (%)'.format(i)] = df.apply(lambda x: x[i]/x['total'], axis=1)
    fig, (ax1, ax2, ax3, ax4) = plt.subplots(4,1)

    ax1.plot(df.index, df['Match (%)']*100, 'k.')
    ax1.set_title('Match')
    ax1.set_xlabel('Position (base)')
    ax1.set_ylabel('Match (%)')

    ax2.plot(df.index, df['Mismatch (%)']*100, 'k.')
    ax2.set_title('Mismatch')
    ax2.set_xlabel('Position (base)')
    ax2.set_ylabel('Mismatch (%)')

    ax3.plot(df.index, df['Insertion (%)']*100, 'k.')
    ax3.set_title('Insertion')
    ax3.set_xlabel('Position (base)')
    ax3.set_ylabel('Insertion (%)')

    ax4.plot(df.index, df['Soft Clipping (%)']*100, 'k.')
    ax4.set_title('Soft Clipping')
    ax4.set_xlabel('Position (base)')
    ax4.set_ylabel('Soft Clipping (%)')

    fig.set_figheight(30)
    fig.set_figwidth(len(df)/10)

    plt.savefig("%s.png" % prefix, bbox_inches='tight')

if __name__== "__main__":
    gen_mutation_rate_graph()
