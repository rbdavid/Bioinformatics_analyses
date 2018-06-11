#!/home/rbdavid/bin/python
##!/Library/Frameworks/Python.framework/Versions/2.7/bin/python
#------------------
# USAGE:
#------------------

#------------------
# PREAMBLE:
#------------------
import sys
import numpy as np
from Bio import SeqIO
from Bio import AlignIO
from Bio.Align import AlignInfo
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.ticker import MultipleLocator

alignment_file = sys.argv[1]
output_file_name_root = sys.argv[2]
psuedocounts = 0.05	# 1/n where n is 20 for standard amino acid alphabet

#------------------
# MAIN:
#------------------
alignment_results = AlignIO.read(alignment_file,'fasta')
alignment_array = np.array([list(rec) for rec in alignment_results],np.character)
nSequences = len(alignment_results)
nSequences_range = range(nSequences)

alignment_summary_info = AlignInfo.SummaryInfo(alignment_results)
consensus = alignment_summary_info.dumb_consensus()
nResidues = len(consensus)
nResidues_range = range(nResidues)

position_frequency_matrix = alignment_summary_info.pos_specific_score_matrix(consensus,chars_to_ignore = ['X','B','Z'])

nResidue_types = len(list(position_frequency_matrix[0]))
nResidue_types_range = range(nResidue_types)
nResidue_types_keys = np.array(sorted(list(position_frequency_matrix[0])))

#print nResidue_types_keys

pfm = np.full((nResidues,nResidue_types),psuedocounts)
for i in nResidues_range:
	temp = 0
	for j in nResidue_types_keys:
		pfm[i][temp] += position_frequency_matrix[i][j]
		temp += 1

ppm = np.copy(pfm)
ppm /= nSequences

pwm = np.copy(ppm)
b = psuedocounts		# simplified some algebra... b is a numerical value describing the background model; in this case, I am using the simplest background model possible, where there is a 1 in 20 chance (0.05) of observing an amino acid at a specific position; therefore, I divide my observed results by the model (1/0.05)
pwm /= b
pwm = np.log2(pwm)

#------------------
# OUTPUT CONCENSUS, PFM, PPM, PWM 
#------------------

np.savetxt(output_file_name_root+'.consensus.txt',consensus,fmt='%s')
np.savetxt(output_file_name_root+'.pfm.dat',pfm,fmt='%.3f')
np.savetxt(output_file_name_root+'.ppm.dat',ppm,fmt='%f')
np.savetxt(output_file_name_root+'.pwm.dat',pwm,fmt='%f')

#------------------
# Post-analysis of position frequency matrix that highlights which residues are observed at which locations. Calculate the probability of observing a residue that is not descibed by the consensus sequence; 'Variance' away from consensus; ignores residue positions that have 'X' as their consensus residue type
#------------------
dash_index = np.where(nResidue_types_keys=='-')
#print dash_index
with open(output_file_name_root+'.pfm.txt','w') as W, open(output_file_name_root+'.prob_nonconsensus.dat','w') as X, open(output_file_name_root+'.prob_nonconsensus_skipped.dat','w') as Y:
	for i in nResidues_range:
                if np.argmax(pfm[i]) != dash_index:
		        Y.write('Resid ' + '%3s'%(i) + ' Consensus residue: ' + consensus[i] + '    Variance from Consensus: %f \n'%(100*np.sum([pfm[i][j]-psuedocounts for j in nResidue_types_range if j != np.argmax(pfm[i])])/nSequences))
                        
		W.write('Resid ' + '%3s'%(i) + ' Consensus residue: ' + consensus[i] + ' Observed residues:  ')
		X.write('Resid ' + '%3s'%(i) + ' Consensus residue: ' + consensus[i] + '    Variance from Consensus: %f \n'%(np.sum([pfm[i][j]-psuedocounts for j in nResidue_types_range if j != np.argmax(pfm[i])])/nSequences))
		for j in nResidue_types_range:
			if pfm[i][j] > 1.0:
				W.write(nResidue_types_keys[j] + '%5d   ' %(pfm[i][j]))
		W.write('\n')

#------------------
# Calculate Mutual Information
#------------------
numerical_alignment_array = np.zeros((alignment_array.shape),dtype=int)
for seq in nSequences_range:
	for res in nResidues_range:
		if alignment_array[seq][res] in ['X','B','Z']:
			numerical_alignment_array[seq][res] = 9999
		else:
			numerical_alignment_array[seq][res] = np.argwhere(nResidue_types_keys == alignment_array[seq][res])[0][0]

M_ij = np.zeros((nResidues,nResidues))
for res1 in nResidues_range[:-1]:
	for res2 in nResidues_range[res1+1:]:
		p_ij = np.zeros((21,21))	# probability array of finding res type i at position res1 and res type j at position res2
		for seq in nSequences_range:
			if numerical_alignment_array[seq][res1] != 9999 and numerical_alignment_array[seq][res2] != 9999:
				p_ij[numerical_alignment_array[seq][res1]][numerical_alignment_array[seq][res2]] += 1
		p_ij /= nSequences	# turn counts into probabilities
		for i in nResidue_types_range:
			for j in nResidue_types_range:
				if p_ij[i][j] > 0.:
					M_ij[res1][res2] += p_ij[i][j]*np.log2(p_ij[i][j]/(ppm[res1][i]*ppm[res2][j]))
	        M_ij[res2][res1] = M_ij[res1][res2]
        print 'Finished calculating mutual information for residue position ', res1

np.savetxt(output_file_name_root+'.mutual_info.dat',M_ij,fmt='%f')

fig, ax = plt.subplots()
ax.tick_params(which='major',length=6,width=1)
ax.tick_params(which='minor',length=3,width=0.5)
ax.xaxis.set_minor_locator(MultipleLocator(10))
ax.xaxis.set_major_locator(MultipleLocator(100))
ax.yaxis.set_minor_locator(MultipleLocator(10))
ax.yaxis.set_major_locator(MultipleLocator(100))

temp = plt.pcolormesh(nResidues_range,nResidues_range,M_ij,cmap='Blues')
cb1 = plt.colorbar()
cb1.set_label(r'Mutual Information (bits)',size=14)
xlabels = [str(int(x)) for x in temp.axes.get_xticks()[:]]
ylabels = [str(int(y)) for y in temp.axes.get_yticks()[:]]
temp.axes.set_xticks(temp.axes.get_xticks(minor=True)[:]+0.5,minor=True)
temp.axes.set_xticks(temp.axes.get_xticks()[:]+0.5)
temp.axes.set_yticks(temp.axes.get_yticks(minor=True)[:]+0.5,minor=True)
temp.axes.set_yticks(temp.axes.get_yticks()[:]+0.5)
temp.axes.set_xticklabels(xlabels)
temp.axes.set_yticklabels(ylabels)

plt.xlim((-0.5,nResidues+0.5))
plt.ylim((-0.5,nResidues+0.5))
plt.xlabel('Residue Number',size=14)
plt.ylabel('Residue Number',size=14)
ax.set_aspect('equal')
plt.tight_layout()
plt.savefig(output_file_name_root+'.mutual_info.png',dpi=600,transparent=True)
plt.close()

histogram = np.zeros(nResidues)
M_ij_cutoff = np.mean(M_ij)+4*np.std(M_ij)
for i in nResidues_range[:-1]:
        for j in nResidues_range[i+1:]:
                if M_ij[i][j] > M_ij_cutoff:
                        histogram[i] += 1
                        histogram[j] += 1

plt.plot(histogram[:],'b-')
plt.grid(b=True,which='major',axis='both',color='#808080',linestyle='--')
plt.xlabel('Residue Number')
plt.ylabel('High Mutual Information Counts')
plt.savefig(output_file_name_root+'.mi_counts.png',dpi=600,transparent=True)
plt.close()

