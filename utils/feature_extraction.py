import tempfile
import os
from Bio import SeqIO

mathFeature_dir = '/home/alisson/work/MathFeature'

def os_system(cmd):
  print(cmd)
  os.system(cmd)

# Remove invalid sequences and duplicate sequences
# Adapted from
# https://biopython.org/wiki/Sequence_Cleaner
def sequence_cleaner(fasta_file, fasta_output, min_length=0):
  print('fe.sequence_cleaner', fasta_file)
  
  sequences = {}

  # Using the Biopython fasta parse we can read our fasta input
  for seq_record in SeqIO.parse(fasta_file, "fasta"):
    sequence = str(seq_record.seq).upper()

    if (
      len(sequence) >= min_length
      and set(sequence) <= set("GATC")
    ):
      sequences[seq_record.id] = sequence
    #else:
    #  print('Removed:', seq_record.id)

  # Write the clean sequences
  # Create a file in the same directory where you ran this script
  with open(fasta_output, "w+") as output_file:
    # Just read the hash table and write on the file as a fasta format
    for id, sequence in sequences.items():
      output_file.write(">" + id + "\n" + sequence + "\n")

def convert_to_dna(fasta_rna, fasta_dna):
  print('fe.convert_to_dna', fasta_rna)
  
  cmd = "sed '/^[^>]/s/u/t/g' {fasta_rna} | sed '/^[^>]/s/U/T/g' > {fasta_dna}" \
    .format(fasta_rna=fasta_rna, fasta_dna = fasta_dna)
  os_system(cmd)

def calc_length(label, fasta_file, output):
  with tempfile.TemporaryDirectory() as tmpdirname:
    print('fe.calc_length', output)

    tmp1 = '{tmpdirname}/tmp1.csv'.format(tmpdirname=tmpdirname)
    # os_system('rm ' + tmp1)

    cmd = 'printf "nameseq\\tlength\\n" > {tmp1}' \
      .format(tmp1 = tmp1)
    os_system(cmd)

    cmd = 'seqkit fx2tab -lin {fasta_file} >> {tmp1}' \
      .format(fasta_file=fasta_file, tmp1 = tmp1)
    os_system(cmd)
    
    cmd = 'awk \'BEGIN{{FS="\\t";OFS=","}} {{print $1,$2}}\' {tmp1} > {output}' \
      .format(tmp1 = tmp1, output=output)
    os_system(cmd)

def calc_kmers(label, fasta_file, output, size =3):
  print('fe.calc_kmers', output)

  os_system('[ -e {output} ] && rm {output}'.format(output = output))

  cnco = "--no-capture-output"
  cmd = '(echo {size} | (conda run -n mathfeature-terminal {cnco} python {mathFeature_dir}/methods/ExtractionTechniques.py -o {output} -l {label} -t kmer -seq 1 -i {fasta_file} > /dev/null))' \
    .format(size=size, cnco=cnco, mathFeature_dir=mathFeature_dir,
            output=output, label=label, fasta_file=fasta_file)
  os_system(cmd)

def calc_fickettScore(label, fasta_file, output):
  print('fe.calc_fickettScore', output)
  
  os_system('[ -e {output} ] && rm {output}'.format(output = output))

  cmd = 'conda run -n mathfeature-terminal python {mathFeature_dir}/methods/FickettScore.py -o {output} -l {label} -seq 1 -i {fasta_file} > /dev/null' \
    .format(mathFeature_dir = mathFeature_dir, output=output, label=label, fasta_file=fasta_file)
  os_system(cmd)

def calc_codingClass(label, fasta_file, output):
  print('fe.calc_codingClass', output)

  os_system('[ -e {output} ] && rm {output}'.format(output = output))

  cmd = 'conda run -n mathfeature-terminal python {mathFeature_dir}/methods/CodingClass.py -o {output} -l {label} -i {fasta_file} > /dev/null' \
    .format(mathFeature_dir = mathFeature_dir, output=output, label=label, fasta_file=fasta_file)
  os_system(cmd)

def merge(csv_inputs, output):
  print('fe.merge')
  
  os_system('[ -e {output} ] && rm {output}'.format(output = output))

  cnco = "--no-capture-output"
  cmd1 = ';'.join(list(map(lambda input : 'echo ' + input, csv_inputs)))
  cmd2 = 'conda run -n mathfeature-terminal {cnco} python {mathFeature_dir}/preprocessing/concatenate.py -n {n} -o {output}' \
    .format(cnco=cnco, mathFeature_dir=mathFeature_dir, output=output, n=len(csv_inputs))

  os_system('({cmd1}) | ({cmd2}) > /dev/null'.format(cmd1=cmd1, cmd2=cmd2))

def rmfiles(inputs):
  print('fe.rmfiles')

  cmd = 'rm -v ' + (' '.join(inputs))
  os_system(cmd)