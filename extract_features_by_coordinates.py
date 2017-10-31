from Bio import Seq, SeqIO, SeqRecord
from Bio.SeqRecord import SeqRecord
import csv, sys

def extract_by_coords(ref_file, ref_file_type, feature_table, outfile):

    record = SeqIO.read(ref_file, ref_file_type)
    sRNA_records = []
    with open(feature_table, 'r') as csv_in:
        reader = csv.DictReader(csv_in)
        for row in reader:
            start = int(row['start']) - 1
            end = int(row['end'])
            if int(row['strand']) == 1:
                sRNA_seq = record.seq[start:end]
            elif int(row['strand']) == -1:
                sRNA_seq = record.seq[start:end].reverse_complement()
            else: 
                print('Feature not stranded: %s' % row['name'])
                continue
            descrip = 'sRNA len: ' + str(len(sRNA_seq))
            sRNA_record = SeqRecord(sRNA_seq, id=row['name'], description=descrip)
            sRNA_records.append(sRNA_record)

    with open(outfile, 'w') as fasta_out:
        SeqIO.write(sRNA_records, fasta_out, 'fasta')
        print('features written to fasta: ', len(sRNA_records))

if __name__ == '__main__':
    if len(sys.argv) == 5:
         extract_by_coords(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
    else:
         print("Usage: extract_features_by_coordinates.py reference_file_in ref_file_type('fasta' or 'genbank') feature_table.csv outpath")
         sys.exit(0)
