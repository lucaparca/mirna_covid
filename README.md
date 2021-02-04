# mirna_covid
Code snippets associated to the paper by Natarelli, Parca et al. (2021)
The complete procedure, with command line instructions, parameters and final script for the merging of the results of the three different prediction methods can be found at the following link: https://github.com/lucaparca/mirna_covid

Results for three different miRNA-RNA interaction binding site prediction methods (IntaRNA, RNAplex and RNAup) are available. These programs have been run to predict interaction binding site between miRNAs and portion of the SARS-CoV-2 virus genome (5'UTR, with a focus on the leader sequence, 3'UTR and spike-coding region):
results_intarna_covid_3utr.out
results_intarna_covid_5utr.out
results_rnaup_covid_3utr.out  
results_rnaup_covid_5utr.out  
results_rnaplex_covid_3utr.out
results_rnaplex_covid_5utr.out

Moreover the raw output of these methods for the prediction of binding sites between the Spike-coding region and lncRNA are available here. These files have not been automatically merged given the diversity and low overlap of the predictions.
output_intarna_lncrna_covid_spike.out
output_rnaplex_lncrna_covid_spike.out
output_rnaup_lncrna_covid_spike.out  


These results are predicted using the command-line version of the methods with the following commands and options.
<mirna_lenghth>, which determines the lenght of the binding sites between the miRNA and the target RNA sequence, has been set to the miRNA length in the case of miRNAs and to 100 in the case of lncRNA.

IntaRNA:
IntaRNA -t <input_file_query> -q <input_file_target> > <output_file>

RNAplex:
RNAplfold -W <mirna_length> -u <mirna_length> -O --plex_output < <input_file>
RNAplex -l <mirna_length> -q <input_file_query> -t <input_file_target> -a ./

RNAup:
RNAup -w <mirna_length> -b -o -3 -5 --interaction_first < <input_file>


The "merge.py" script merges the miRNA results for the different methods. This script has been commented in order to follow the logic of the different parts.
It requires the python "numpy" and "itertools" libraries.
