"""
EBprotV2_PLOT.R coded in Python
"""
import glob
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from os import listdir
import pandas as pd
from math import ceil

def getFig(fm, nrow):
    if (len(fm)==1):
        fig = plt.figure(figsize=(6,6))
        subplotsize=1
    elif (nrow>1):
        fig = plt.figure(figsize=(12,12))
        subplotsize=4
    elif (nrow==1 and len(fm)>1):
        fig = plt.figure(figsize=(12,6))
        subplotsize=2
    return fig,subplotsize

def getFig2(nums):
    if (nums==1):
        fig = plt.figure(figsize=(6,6))
        subplotsize = 1
    elif(nums>2):
        fig = plt.figure(figsize=(12,12))
        subplotsize = 4
    elif(nums==2):
        fig = plt.figure(figsize=(12,6))
        subplotsize = 2
    return fig,subplotsize

if __name__ == "__main__":

	###density plots
	fm = glob.glob('density*.txt')
	fm.sort()
	nrow = ceil(len(fm)/2)
	subplotsize=0
	fig, subplotsize = getFig(fm,nrow)
	with PdfPages("DensityPlots_EBprotV2.pdf") as pdf:
		for i in range(len(fm)):
			label = fm[i].strip('.txt').strip('densityfile_')
			denfile = pd.read_csv(fm[i], sep = '\t')
			if (subplotsize == 1):
				ax = fig.add_subplot(1,1,1)
			elif (subplotsize == 2):
				ax = fig.add_subplot(1,2,i%subplotsize+1)
			elif (subplotsize == 4):
				ax = fig.add_subplot(2,2,i%subplotsize+1)
			ax.plot(denfile.Ratio, denfile.Density, 'k', label = "overall density", linewidth=2)
			ax.plot(denfile.Ratio, denfile.Null_density, "tab:gray", label = "null-density", linewidth=1.5)
			ax.plot(denfile.Ratio, denfile.Pos_density, 'r', label = "up-density", linewidth=1.5)
			ax.plot(denfile.Ratio, denfile.Neg_density, 'b', label = "down-density", linewidth=1.5)
			ax.set_title(label,fontweight='bold')
			ax.set_xlabel("Ratio")
			ax.set_ylabel("Density")
			ax.legend(loc = "upper left")
			if (i%subplotsize+1 == subplotsize or i==len(fm)-1):
				pdf.savefig()
				plt.clf()
		plt.close()

	###PP plots
	results = pd.read_csv("EBprot_results.txt", sep = '\t')
	###if %6 == 2 then it's rep case
	nums = int((len(results.columns)-8)/6) if len(results.columns)%6 == 2 else int((len(results.columns)-1)/6)
	subplotsize=0
	fig, subplotsize = getFig2(nums)
	with PdfPages("EBprotV2_PPscoreplot.pdf") as pdf:
		for i in range(1,nums+1):
			tmp_l = results.columns[(i-1)*6+1].split('_')[1]
			index_r = ((i-1)*6)+1
			index_p = ((i-1)*6)+4
			index_b = (i*6)
			upperL = np.nanmax(results.iloc[:,index_r]) + 0.05
			lowerL = np.nanmin(results.iloc[:,index_r]) - 0.05
			if (subplotsize == 1):
				ax = fig.add_subplot(1,1,1)
			elif (subplotsize == 2):
				ax = fig.add_subplot(1,2,(i-1)%subplotsize+1)
			elif (subplotsize == 4):
				ax = fig.add_subplot(2,2,(i-1)%subplotsize+1)
			ax.scatter(results.iloc[:,index_r], results.iloc[:,index_p], s = 5, c = 'k', label='_nolegend_')
			ax.set_title(tmp_l,fontweight='bold')
			ax.set_xlabel("MedianLog2ratio")
			ax.set_ylabel("PPscore")
			ax.set_xlim(lowerL, upperL)
			ax.set_ylim(-1.05,1.05)
			
			cidup = np.intersect1d(np.where(results.iloc[:,index_b] < 0.05)[0],np.where(results.iloc[:,index_r] > 0)[0])
			ciddown = np.intersect1d(np.where(results.iloc[:,index_b] < 0.05)[0],np.where(results.iloc[:,index_r] < 0)[0])
			###handling size 0 cases
			if (cidup.size!=0):
				cutoffup = np.nanmin(np.abs(results.iloc[cidup,index_p]))
				ax.axhline(cutoffup, linestyle='dashed', color = 'g')
			if (ciddown.size!=0):
				cutoffdown = np.nanmin(np.abs(results.iloc[ciddown,index_p]))
				ax.axhline(-1*cutoffdown, linestyle='dashed', color = 'g')
			if (cidup.size!=0 or ciddown.size!=0):
				conc = np.concatenate((cidup,ciddown))
				ax.scatter(results.iloc[conc,index_r], results.iloc[conc,index_p], s = 5, c = 'r', label ='sig at 5% FDR')
        			
			ax.legend(loc = "upper left", frameon=False)
			
			if ((i-1)%subplotsize+1 == subplotsize or i==nums):
				pdf.savefig()
				plt.clf()
		plt.close()

