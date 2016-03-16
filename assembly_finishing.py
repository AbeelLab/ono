# -*- coding: utf-8 -*-
"""
Created on Sun Feb 22 12:14:15 2015

@author: jasperlinthorst and c_georgescu
"""

import GSA
import matplotlib.pyplot as plt
import os
import argparse
import pickle
import assembly_finishing_objects

def output_to_fasta(output, compt, seq, prune):
	n = False
	for char in seq:
		if (prune == True and (char == 'n' or char == 'N')):
			if (n == True):
				continue
			n = True
		else:
			n = False
		if compt >= 60:
			output.write('\n')
			compt = 0
		output.write(char)
		compt += 1
	return compt

def fasta_reader(fn,truncN=False):
    seq=""
    with open(fn,'r') as ff:
        for line in ff:
            if line.startswith(">"):
                if seq!="":
                    yield name,seq
                name=line.rstrip().replace(">","")
                seq=""
            else:
                if truncN:
                    for base in line.rstrip():
                        if base.upper()=='N':
                            if len(seq)==0:
                                continue
                            elif seq[-1]=='N':
                                continue
                            else:
                                seq+='N'
                        else:
                            seq+=base
                else:
                    seq+=line.rstrip()
        if seq!="":
            yield name,seq

def fasta_writer(fn,name_seq,lw=100):
    seq=""
    with open(fn,'w') as ff:
        for name,seq in name_seq:
            if not name.startswith(">"):
                name=">"+name+"\n"
            ff.write(name)
            for i in range( (len(seq)/lw)+(len(seq) % lw > 0)):
                ff.write(seq[i*lw:(i+1)*lw]+"\n")

def main():
	usage = "usage: assembly_finishing.py [options] <reference.fasta> <scaffolds.fasta>"
	parser = argparse.ArgumentParser(usage)
	parser.add_argument('fastas', nargs='*', help='Two fasta files for which mumplot is to be made.')
	parser.add_argument("--minmum", type=int, dest="minmum", default=1000, help="Mums below this size wont be used for the analysis nor shown in the plot.")
	parser.add_argument("-i","--interactive", action="store_true", dest="interactive", default=False, help="Whether or not to show interactive plots.")
	parser.add_argument("-p", action="store_true", dest="ploting", default=False, help="Use to output an image of the mumplot in png format.")
	parser.add_argument("-b", action="store_true", dest="breakgaps", default=False, help="Break the contigs at gaps, before scaffolding.")
	parser.add_argument("-o", default = None, help = "The name of the output.")
	parser.add_argument("-n", type = int, default = 1000, help = "Number of 'N's to be inserted between concatenated contigs.")
	parser.add_argument("-step", type = int, default = 1000, help = "How far to look for neighbour MUMs.")
	parser.add_argument("-smallest", type = int, default = 1000, help = "Smallest allowed MUM size for contigs that only have that MUM.")
	parser.add_argument("-discard", action="store_true", dest = "discard", default = False, help = "If used, discards the contigs that couldn't be aligned (with the current parameters) instead of appending them to the end of the output file.")
	parser.add_argument("-prune", action="store_true", dest = "prune", default = True)
	args = parser.parse_args()
	
	if len(args.fastas)!=2:
		parser.error('Specify the fasta file to finish assembling and the reference fasta file.')
	
	f1=args.fastas[0]
	f2=args.fastas[1]
	
        if args.breakgaps:
            for f in [f1,f2]:
                contigs=[]
                for name,seq in fasta_reader(f,truncN=True):
	            seqs=seq.split('N')
                    for seq in seqs:
                        contigs.append((name+" ("+str(len(seq))+")",seq))
                print f,len(contigs)
                fasta_writer('_'+f,contigs)
	    f1='_'+f1
	    f2='_'+f2
        
	global index
	index=GSA.index(f1,f2,1)
	g=index.graph
	mums=index.get_mums(args.minmum)

##########################################################

	seq = assembly_finishing_objects.Sequence(index, mums, args.step, args.smallest)
		
	if (seq.contigs[0].futur == -1):
		print "Nothing to do. Exiting."
		return

	print "New sequence generated", len(seq.contigs)
	print "Writing to file....",

	if (args.o is None):
		new_seq = "treated_" + g.vertices[1].origin.copy().pop()
	else:
		new_seq = args.o

	output = open(new_seq, 'w')
	output.write(">" + g.vertices[1].contig_origin.copy().pop() + "\n")
	compt = 0
	rolling = seq.contigs[0].mum_sequences[0].mums[0][1]
	v = g.vertices.values()[seq.contigs[0].id]
	if (seq.contigs[0].futur == 0):
		compt = output_to_fasta(output, compt, index.T[rolling : v.saoffset + v.contig_end - v.contig_start], args.prune)
	else:
		compt = output_to_fasta(output, compt, index.T[2*v.rcsaoffset - rolling : v.rcsaoffset + v.contig_end - v.contig_start], args.prune)
	
        compt = output_to_fasta(output, compt, "X"*args.n, args.prune)
	for c in seq.contigs[1:]:
		# [1]+1 for id
		# [2] to know if to reverse or not (take rcsaoffset instead of saoffset)
		# v = g.vertices.values()[c[1]+1]
		v = g.vertices.values()[c.id]
		if c.futur == 0 or c.futur == -1:
			if (c.futur == -1 and args.discard == True):
				break
			compt = output_to_fasta(output, compt, index.T[v.saoffset:v.saoffset + v.contig_end - v.contig_start], args.prune)
		else:
			compt = output_to_fasta(output, compt, index.T[v.rcsaoffset:v.rcsaoffset + v.contig_end - v.contig_start], args.prune)
		compt = output_to_fasta(output, compt, "Y"*args.n, args.prune)
	
        v = g.vertices.values()[seq.contigs[0].id]
	if (seq.contigs[0].futur == 0):
		compt = output_to_fasta(output, compt, index.T[v.saoffset : rolling], args.prune)
	else:
		compt = output_to_fasta(output, compt, index.T[v.rcsaoffset : 2*v.rcsaoffset - rolling], args.prune)
	output.write("\n")
	if (args.discard == False):
		for c in seq.discarded_contigs:
			v = g.vertices.values()[c.id]
			output.write(">" + v.contig_origin.copy().pop() + "\n")
			output_to_fasta(output, 0, index.T[v.saoffset:v.saoffset + v.contig_end - v.contig_start], args.prune)
			output.write("\n")
	output.close()
	
	print "Done."

        if (args.ploting == True):
            plot(f1,args,f2)
            plot(f1,args,new_seq)

def plot(f1,args,new_seq):
        #############################################################
        ################ Restarting the analysis ####################
        #############################################################
                index=GSA.index(f1,new_seq,1)
                g=index.graph
                mums=index.get_mums(args.minmum)

	#    import collinearity
	#    mums,falsemums=collinearity.filter_collinear_mums(mums)
		
		ax=plt.subplot(111)
		x=0
		y=0
		plot_offsets={}
		contigs1=[]
		contigs2=[]
		
		#for every vertex/contig in the graph create a domain on the x or y axis
		for vertex in g.vertices.values():
			if vertex.input_origin==0: #y-axis
				l=vertex.contig_end-vertex.contig_start
				plot_offsets[vertex]=y
				y+=l
				contigs1.append(vertex)
			else: #x-axis for input_origin == 1
				l=vertex.contig_end-vertex.contig_start
				plot_offsets[vertex]=x
				x+=l
				contigs2.append(vertex)
		
		#Create a dictionary with all oberserved contigpairs as keys and a list of mums as values
		mums_per_contigpair={}
		for mum in mums:
			if mums_per_contigpair.has_key((mum[3],mum[4])):
				start1=mum[0]-mum[3].saoffset+plot_offsets[mum[3]]
				start2=mum[1]-mum[4].saoffset+plot_offsets[mum[4]]
				mums_per_contigpair[(mum[3],mum[4])].append((start1,start2,mum[2],mum[5],0))
			else:
				start1=mum[0]-mum[3].saoffset+plot_offsets[mum[3]]
				start2=mum[1]-mum[4].saoffset+plot_offsets[mum[4]]
				mums_per_contigpair[(mum[3],mum[4])]=[(start1,start2,mum[2],mum[5],0)]
		
		plot_false_mums=False
		if plot_false_mums:
			for mum in falsemums:
				if mums_per_contigpair.has_key((mum[3],mum[4])):
					start1=mum[0]-mum[3].saoffset+plot_offsets[mum[3]]
					start2=mum[1]-mum[4].saoffset+plot_offsets[mum[4]]
					mums_per_contigpair[(mum[3],mum[4])].append((start1,start2,mum[2],mum[5],1))
				else:
					start1=mum[0]-mum[3].saoffset+plot_offsets[mum[3]]
					start2=mum[1]-mum[4].saoffset+plot_offsets[mum[4]]
					mums_per_contigpair[(mum[3],mum[4])]=[(start1,start2,mum[2],mum[5],1)]
			
		
		for c1 in contigs1:
			for c2 in contigs2:
				if mums_per_contigpair.has_key((c1,c2)):
					for mum in mums_per_contigpair[(c1,c2)]:
						#draw a line!
						if mum[4]==0: #correct mums
							if mum[3]==1: #rc, so draw opposite direction
								x=plt.plot([mum[0], mum[0]+mum[2]], [mum[1]+mum[2], mum[1]], color='g', linestyle='-', linewidth=1, rasterized=True)
							else:
								x=plt.plot([mum[0],mum[0]+mum[2]], [mum[1],mum[1]+mum[2]], color='r', linestyle='-', linewidth=1, rasterized=True)
						else:
							if mum[3]==1: #rc, so draw opposite direction
								x=plt.plot([mum[0], mum[0]+mum[2]], [mum[1]+mum[2], mum[1]], color='k', linestyle='-', linewidth=1, rasterized=True)
							else:
								x=plt.plot([mum[0],mum[0]+mum[2]], [mum[1],mum[1]+mum[2]], color='k', linestyle='-', linewidth=1, rasterized=True)

		xticks=[plot_offsets[v]+(v.contig_end-v.contig_start) for v in contigs1]
		xticklabels=[v.contig_origin.copy().pop() for v in contigs1]
		plt.xticks(xticks, xticklabels, rotation='25', fontsize=8)
		yticks=[plot_offsets[v]+(v.contig_end-v.contig_start) for v in contigs2]
		yticklabels=[v.contig_origin.copy().pop() for v in contigs2]
		plt.yticks(yticks, yticklabels, rotation='horizontal', fontsize=8)
		plt.xlim(0, max(xticks))
		plt.ylim(0, max(yticks))
		plt.xlabel(os.path.basename(f1))
		plt.ylabel(os.path.basename(new_seq))
		#ax=plt.twinx()
		
		plt.grid(linestyle='-',linewidth=1.)
		plt.title(os.path.basename(f1)+' vs. '+os.path.basename(new_seq))
		# pickle.dump(ax, file(os.path.basename(f1)+'_'+os.path.basename(f2)+'.mumplot.pickle', 'w'))
		plt.savefig(file(os.path.basename(f1)+'_'+os.path.basename(new_seq)+'.png','w'))

		if args.interactive==True:
			plt.show()
                plt.clf()





if __name__ == "__main__":
	main()



# replacement for smallest
