#This example is using GRCh38.
#The GRCh37.fa is not provided in this directory due to the size of it.
#If you install INTEGRATE-Neo and have the setup.ini set correctly, you should be able to run this:
#(but make sure to change the "/PATH/TO"s to proper paths before running).

python /PATH/TO/integrate-neo.py -t HLA_alleles.tsv -f fusions.bedpe -r /PATH/TO/GRCh37.fa -g /PATH/TO/GRCh37.genePhred -k

#and get the output result.bedpe in the fusion_antigen_out directory. 
