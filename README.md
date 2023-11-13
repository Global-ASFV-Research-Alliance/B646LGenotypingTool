# B646LGenotypingTool

Note that this tool requires that BLAST is installed on your computer, and the custom file includes the P72 genotypes. The .\\23-11-13_Representative_p72.xlsx file must always be referenced as excel for usage. Query should be the fasta file you are examining.

### Help Printout

options:
  
  -i, --input, Input Fasta File Location
  
  -e, --excel, Input Excel File Location 

### Example Usage
This is an example of using the files located within this folder. The ASFVG file is an example one - the Representative_p72 file should always be used.

    python .\BlastAlgorithm.py -i .\ASFVG_protein.fa -e .\23-11-13_Representative_p72.xlsx
  
