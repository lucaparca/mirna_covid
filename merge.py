import os,sys,numpy as np
from itertools import combinations

#pick matches for either 3'UTR or 5'UTR

utr='5'
intarna=open('results_intarna_covid_%sutr.out'%(utr)).readlines()
rnaup=open('results_rnaup_covid_%sutr.out'%(utr)).readlines()
rnaplex=open('results_rnaplex_covid_%sutr.out'%(utr)).readlines()

#results coming from the different methods will be stored in these dictionaries
intarna_d={}
rnaup_d={}
rnaplex_d={}

#read and store the results for the three methods
for line in intarna:
	mirna,energy,start,end,start2,end2,region=line.split()
	if mirna not in intarna_d:intarna_d[mirna]={}
	intarna_d[mirna]['energy']=float(energy)
	intarna_d[mirna]['start']=start
	intarna_d[mirna]['end']=end
	intarna_d[mirna]['start2']=start2
	intarna_d[mirna]['end2']=end2
	intarna_d[mirna]['region']=region

for line in rnaup:
	mirna,energy,start,end,start2,end2,region=line.split()
	if mirna not in rnaup_d:rnaup_d[mirna]={}
	rnaup_d[mirna]['energy']=float(energy)
	rnaup_d[mirna]['start']=start
	rnaup_d[mirna]['end']=end
	rnaup_d[mirna]['start2']=start2
	rnaup_d[mirna]['end2']=end2
	rnaup_d[mirna]['region']=region

for line in rnaplex:
	mirna,energy1,energy2,start,end,start2,end2,region=line.split()
	if mirna not in rnaplex_d:rnaplex_d[mirna]={}
	rnaplex_d[mirna]['energy1']=float(energy1)
	rnaplex_d[mirna]['energy2']=float(energy2)
	rnaplex_d[mirna]['start']=start
	rnaplex_d[mirna]['end']=end
	rnaplex_d[mirna]['start2']=start2
	rnaplex_d[mirna]['end2']=end2
	rnaplex_d[mirna]['region']=region

#results will be written here
out=open('complete_results_%sUTR.tab'%(utr),'w')
out.write('\t'.join(['mirna','intarna_energy','rnaup_energy','rnaplex_energy1','rnaplex_energy2','mean_energy','intarna_covid_start','intarna_covid_end','rnaup_covid_start','rnaup_covid_end','rnaplex_covid_start','rnaplex_covid_end','intarna_mirna_start','intarna_mirna_end','rnaup_mirna_start','rnaup_mirna_end','rnaplex_mirna_start','rnaplex_mirna_end','covid_site_jaccard','mirna_site_jaccard','covid_site_jaccard_best_method_pair','mirna_site_jaccard_best_method_pair','intarna_region','rnaup_region','rnaplex_region'])+'\n')

for mirna in set(list(intarna_d.keys()))|set(list(rnaup_d.keys()))|set(list(rnaplex_d.keys())):
	mean_energy=[]
	covid_starts=[]
	covid_ends=[]
	mirna_starts=[]
	mirna_ends=[]
	#read results miRNA-wise from the different methods
	if mirna in intarna_d:
		energy1=intarna_d[mirna]['energy'];mean_energy.append(energy1)
		covid_start1=intarna_d[mirna]['start'];covid_starts.append(int(covid_start1))
		covid_end1=intarna_d[mirna]['end'];covid_ends.append(int(covid_end1))
		mirna_start1=intarna_d[mirna]['start2'];mirna_starts.append(int(mirna_start1))
		mirna_end1=intarna_d[mirna]['end2'];mirna_ends.append(int(mirna_end1))
		region1=intarna_d[mirna]['region']
	else:energy1='NA';covid_start1='NA';covid_end1='NA';mirna_start1='NA';mirna_end1='NA';region1='NA'
	if mirna in rnaup_d:
		energy2=rnaup_d[mirna]['energy'];mean_energy.append(energy2)
		covid_start2=rnaup_d[mirna]['start'];covid_starts.append(int(covid_start2))
		covid_end2=rnaup_d[mirna]['end'];covid_ends.append(int(covid_end2))
		mirna_start2=rnaup_d[mirna]['start2'];mirna_starts.append(int(mirna_start2))
		mirna_end2=rnaup_d[mirna]['end2'];mirna_ends.append(int(mirna_end2))
		region2=rnaup_d[mirna]['region']
	else:energy2='NA';covid_start2='NA';covid_end2='NA';mirna_start2='NA';mirna_end2='NA';region2='NA'
	if mirna in rnaplex_d:
		energy3=rnaplex_d[mirna]['energy1'];mean_energy.append(energy3)
		energy4=rnaplex_d[mirna]['energy2']
		covid_start3=rnaplex_d[mirna]['start'];covid_starts.append(int(covid_start3))
		covid_end3=rnaplex_d[mirna]['end'];covid_ends.append(int(covid_end3))
		mirna_start3=rnaplex_d[mirna]['start2'];mirna_starts.append(int(mirna_start3))
		mirna_end3=rnaplex_d[mirna]['end2'];mirna_ends.append(int(mirna_end3))
		region3=rnaplex_d[mirna]['region']
	else:energy3='NA';covid_start3='NA';covid_end3='NA';mirna_start3='NA';mirna_end3='NA';region3='NA'
	#calculate the mean energy for the binding site (if possible, otherwise NA)
	try:mean_energy=np.mean(mean_energy)
	except:mean_energy='NA'

	#calculate the Jaccard index between the interaction ranges on the viral genome sequence
	#calculate the Jaccard index when there are 3 complete matches in 3 different methods
	if len(covid_starts)==3:
		covid_range1=list(range(min([covid_starts[0],covid_ends[0]]),max([covid_starts[0],covid_ends[0]])))
		covid_range2=list(range(min([covid_starts[1],covid_ends[1]]),max([covid_starts[1],covid_ends[1]])))
		covid_range3=list(range(min([covid_starts[2],covid_ends[2]]),max([covid_starts[2],covid_ends[2]])))
		i=len(set(covid_range1)&set(covid_range2)&set(covid_range3))
		u=len(set(covid_range1)|set(covid_range2)|set(covid_range3))
		j1=float(i)/u

		j1_secondary=-1
		#calculate the best Jaccard considering pairs of methods
		for pair in combinations([covid_range1,covid_range2,covid_range3],2):
			i=len(set(pair[0])&set(pair[1]))
			u=len(set(pair[0])|set(pair[1]))
			j_tmp=float(i)/u
			if j_tmp>j1_secondary:j1_secondary=j_tmp

	#same as above but only 2 matches are available out of 3 methods
	elif len(covid_starts)==2:
		covid_range1=list(range(min([covid_starts[0],covid_ends[0]]),max([covid_starts[0],covid_ends[0]])))
		covid_range2=list(range(min([covid_starts[1],covid_ends[1]]),max([covid_starts[1],covid_ends[1]])))
		i=len(set(covid_range1)&set(covid_range2))
		u=len(set(covid_range1)|set(covid_range2))
		j1=float(i)/u
		j1_secondary=j1
	else:j1='NA'

	#same as above but for the ranges on the miRNA sequence
	if len(mirna_starts)==3:
		mirna_range1=list(range(min([mirna_starts[0],mirna_ends[0]]),max([mirna_starts[0],mirna_ends[0]])))
		mirna_range2=list(range(min([mirna_starts[1],mirna_ends[1]]),max([mirna_starts[1],mirna_ends[1]])))
		mirna_range3=list(range(min([mirna_starts[2],mirna_ends[2]]),max([mirna_starts[2],mirna_ends[2]])))
		i=len(set(mirna_range1)&set(mirna_range2)&set(mirna_range3))
		u=len(set(mirna_range1)|set(mirna_range2)|set(mirna_range3))
		j2=float(i)/u
		j2_secondary=-1
		for pair in combinations([mirna_range1,mirna_range2,mirna_range3],2):
			i=len(set(pair[0])&set(pair[1]))
			u=len(set(pair[0])|set(pair[1]))
			j_tmp=float(i)/u
			if j_tmp>j2_secondary:j2_secondary=j_tmp

	elif len(mirna_starts)==2:
		mirna_range1=list(range(min([mirna_starts[0],mirna_ends[0]]),max([mirna_starts[0],mirna_ends[0]])))
		mirna_range2=list(range(min([mirna_starts[1],mirna_ends[1]]),max([mirna_starts[1],mirna_ends[1]])))
		i=len(set(mirna_range1)&set(mirna_range2))
		u=len(set(mirna_range1)|set(mirna_range2))
		j2=float(i)/u
		j2_secondary=j2
	else:j2='NA'

	out.write('\t'.join([mirna,str(energy1),str(energy2),str(energy3),str(energy4),str(mean_energy),covid_start1,covid_end1,covid_start2,covid_end2,covid_start3,covid_end3,mirna_start1,mirna_end1,mirna_start2,mirna_end2,mirna_start3,mirna_end3,str(j1),str(j2),str(j1_secondary),str(j2_secondary),region1,region2,region3])+'\n')
out.close()
