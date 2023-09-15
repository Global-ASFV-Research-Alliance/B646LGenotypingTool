import pandas as pd
from Bio import SeqIO

def p72finder(queryfile, subjectfile, GenomeColumn = 'Genome', SequenceColumn = 'Sequence', IsolateColumn = 'Isolate', GenotypeColumn = 'GenotypeNumber', outputname_fasta = 'P72.fasta', outputlocation = 'blastxout.txt'):
    
    subject_excel = pd.read_excel(subjectfile, usecols=[GenomeColumn,SequenceColumn,IsolateColumn,GenotypeColumn])
    with open(outputname_fasta, "w") as f:
        pass

    for index, row in subject_excel.iterrows():      
        with open(outputname_fasta, "a") as f:
            print(">" + row[GenomeColumn] + "\n" + str(row[SequenceColumn]).replace("[", "").replace("]", "").replace("'", ""), file = f)

    #Identify if Nucleotide or Ammino Acid Sequence
    
    for record in SeqIO.parse(queryfile, "fasta"):
        ok = r'ACTGUWSMKRYBDHVN*'
        Type = all(c in ok for c in str(record.seq))
    
    import os
    if Type == True:
        from Bio.Blast.Applications import NcbiblastxCommandline
        blastcmdline = NcbiblastxCommandline(cmd='blastx', query = queryfile, subject= outputname_fasta, max_hsps = 1, max_target_seqs = 43, out = outputlocation, outfmt = "6 qseqid sseqid evalue pident bitscore length qlen slen qstart qend sstart send qseq sseq", evalue=0.001)
        os.system(str(blastcmdline))
    if Type == False:
        from Bio.Blast.Applications import NcbiblastpCommandline
        blastcmdline = NcbiblastpCommandline(cmd='blastp', query = queryfile, subject= outputname_fasta, max_hsps = 1, max_target_seqs = 43, out = outputlocation, outfmt = "6 qseqid sseqid evalue pident bitscore length qlen slen qstart qend sstart send qseq sseq", evalue=0.001)
        os.system(str(blastcmdline))

    #ReadBlastResults
    prediction = pd.read_csv(outputlocation, sep = '\t', names=('qseqid', 'sseqid', 'evalue', 'pident', 'bitscore', 'length', 'qlen', 'slen', 'qstart', 'qend', 'sstart', 'send', 'qseq', 'sseq')).drop_duplicates(["qseq"])
    result = prediction.merge(subject_excel, left_on='sseqid', right_on=GenomeColumn)
    print('Predicted Genotype: ' + str(result.iloc[0][GenotypeColumn]) + "\n" + "The B646L(P72) encoded by your genome is " + str(result.iloc[0]['pident']) + "%" + " identical to the " + str(result.iloc[0]['length']) + ' ammino acids in ' + str(result.iloc[0][IsolateColumn]) + ' (' + str(result.iloc[0][GenomeColumn]) + ')')
    
    if result.iloc[0]['length'] <= 600:
        print("WARNING!!! No full length B646L found in the submitted sequence. A manual double-check of raw blast output is recommended, and genotyping results are highly likely to be illegitamite.")
    elif result.iloc[0]['length'] <= 635:
        print("Warning - No full length B646L found in submitted sequence. A manual double-check of the raw blast output is recommended, there are potentially multiple genotypes overlapping.")
