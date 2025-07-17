import os, sys, io
import numpy as np
from subprocess import call
from Bio import Phylo, Nexus, AlignIO, SeqIO
import random

def writeDateTrait(dates,  output):

	traits = []	
	for i in dates.keys():
		traits.append('\n\t\t\t\t\t%s=%.10f' %(i, dates[i]))
	output.write(','.join(traits)+'\n')
#

def writeAlignment(newicktree, output):

	tree   = Nexus.Trees.Tree(newicktree)	
	taxa   = tree.get_taxa()

	for i in range(0,len(taxa)):
		output.write('\n\t\t\t\t\t<sequence spec="Sequence" taxon="%d" value="?" />' % (i+1))
#


# Read sampling times from tree
# 
# If forwards = True  gives times as time from tMRCA
# If forwards = False gives times as time from most recent sample
#
def getSamplingTimes(newicktree, forwards=False):

	# Sampling times from tree
	tree    = Nexus.Trees.Tree(newicktree)
	leaves  = tree.get_terminals()
	heights = np.zeros(len(leaves))	
	for i in range(0,len(leaves)):
		heights[i] = tree.sum_branchlength(node=leaves[i])

	# Forwards or backwards
	if (forwards == True):
		treetimes = heights
	else:
		treetimes = max(heights) - heights

	# Associate with label
	times = dict()
	for i in range(0,len(treetimes)):
		label = tree.get_taxa(node_id=leaves[i])[0]
		times[label] = treetimes[i]		
	#	

	return times	
#


def makeXMLFile(pars, template, outputfile="", outputpath=""):

	if (outputpath == ""):
		outputpath = pars['outputpath']

	if (not os.path.exists(outputpath)):
		os.makedirs(outputpath)

	if (outputfile == ""):
		outputfile = pars["name"]

	sys.stdout.write(outputfile+"...\n")

	formatpars = dict()
	for par in pars:
		formatpars['$'+par] = pars[par]
	output = template.format(**formatpars)

	outfile = open(outputpath+"/"+outputfile+".xml", 'w')
	outfile.write(output)
	outfile.close()

def subsample_fasta_columns(input_fasta, output_fasta, max_sites):
    records = list(SeqIO.parse(input_fasta, "fasta"))
    if not records:
        raise ValueError("No sequences found in FASTA file: " + input_fasta)

    alignment_length = len(records[0].seq)
    keep_sites = list(range(alignment_length))

    # subsamples if there are more sites than max_sites, otherwise keeps all sites
    if alignment_length > max_sites:
        keep_sites = sorted(random.sample(keep_sites, max_sites))

    # Subsample and write to output
    with open(output_fasta, "w") as out_f:
        for record in records:
            record.seq = record.seq.__class__("".join([record.seq[i] for i in keep_sites]))
            SeqIO.write(record, out_f, "fasta")

    return len(keep_sites)/alignment_length # Return fraction of sites kept



def pars_alignment(aln_path, pars):
    """
    Read a FASTA alignment and add it to pars["taxa"] in BEAST-compatible XML format.
    """
    alignment = AlignIO.read(aln_path, "fasta")
    output_align = io.StringIO()

    for record in alignment:
        output_align.write(f'\n\t\t<sequence><taxon idref="{record.id}"/>{str(record.seq)}</sequence>')

    pars["taxa"] = output_align.getvalue()

def pars_dates(tree, pars):
    """
    Add dates to the tree in the BEAST XML parameters.
    """
    output_dates = io.StringIO()

    tree_ = Phylo.read(io.StringIO(tree), "newick")
    leaves = tree_.get_terminals()

    # Get homochronous setting from pars (default to False if not present)
    is_homochronous = pars.get("homochronous", False)

    if is_homochronous:
        for leaf in leaves:
            output_dates.write(f'\n\t\t<taxon id = "{leaf.name}"> <date value="0.0" direction="backwards" units="years"/> </taxon>')
    else:
        heights = np.array([tree_.distance(leaf) for leaf in leaves])
        treetimes = np.max(heights) - heights
        for i, leaf in enumerate(leaves):
            output_dates.write(f'\n\t\t<taxon id = "{leaf.name}"> <date value="{treetimes[i]}" direction="backwards" units="years"/> </taxon>')

    pars["dates"] = output_dates.getvalue()

def get_beast_call(beast, filename, seed):	
	return ["java", "-jar", beast, "-seed", str(seed), "-overwrite", filename]
