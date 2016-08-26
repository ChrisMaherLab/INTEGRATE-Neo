#This example is using GRCh37.
#The GRCh37.fa is not provided in this directory due to the size of it.  
#If you installed INTEGRATE-Neo and had the setup.ini setted correcty, you should be able to run this:
#(but make sure to change the "/PATH/TO"s to proper paths before running).

python /PATH/TO/integrate-neo.py -t HLA_alleles.tsv -f fusions.bedpe -r /PATH/TO/GRCh37.fa -g transcripts.tsv -k

#and get the output result.bedpe in the fusion_antigen_out directory. 
