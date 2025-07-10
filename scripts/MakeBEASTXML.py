import os, sys, io, yaml
import numpy as np
from fnmatch import fnmatch
from subprocess import call
from optparse import OptionParser

from beastutils import *

usage = "usage: %prog [option]"
parser = OptionParser(usage=usage)
seqGenPath = "../../Seq-Gen-1.3.5/source/seq-gen"  # Path to seq-gen executable, adjust if necessary

parser.add_option("-i","--inputpath",
                  dest = "inputpath",
                  default = "",
                  metavar = "path",
                  help = "Path to input files [required]")

parser.add_option("-c","--config",
                  dest = "config",
                  default = "*.cfg",
                  metavar = "",
                  help = "Pattern to match for config files [default = %default]")

parser.add_option("-b","--beast",
                  dest = "beast",
                  default = "",
                  metavar = "path",
                  help = "Path to BEAST2 jar file [required]")

parser.add_option("-x","--template",
                  dest = "template",
                  default = "",
                  metavar = "path",
                  help = "Path to template XML file [required]")

parser.add_option("-o","--outputpath",
                  dest = "outputpath",
                  default = "",
                  metavar = "path",
                  help = "Path to save output file in [required]")

parser.add_option("-n","--name",
                  dest = "name",
                  default = "",
                  metavar = "path",
                  help = "Name of the runs [required]")

parser.add_option("-s","--seed",
                  dest = "seed",
                  default = "127",
                  metavar = "integer",
                  help = "Seed or comma separated list of seeds to use [required]")

(options,args) = parser.parse_args()

if (options.inputpath != ""):
	config         = options.config
	inputpath      = os.path.abspath(options.inputpath)+"/"
else:
	config         = options.config[options.config.rfind("/")+1:]
	inputpath      = os.path.abspath(options.config[:options.config.rfind("/")])+"/"


################################################################################################################################  


	

################################################################################################################################  

for filename in sorted(os.listdir(inputpath)):
    if (fnmatch(filename,config)):

        sys.stdout.write(filename+"\t"+config+"\n")

        # Load config file
        configfile = open(inputpath+filename, 'r').read().replace("\t"," ")
        pars 	   = yaml.load(configfile, Loader=yaml.SafeLoader)

        # Set BEAST specific parameters	
        seeds      = list(map(int, options.seed.split(',')))				
        basename   = pars["name"] if options.name == '' else pars["name"]+"_"+options.name
        outputpath = os.path.abspath(pars["outputpath"] if options.outputpath == '' else options.outputpath)
        template   = open(os.path.abspath(pars["template"] if options.template == '' else options.template), 'r').read()
        treefile   = open(pars["trees"], 'r')

        # Output scripts
        if (not os.path.exists(outputpath)):
            os.makedirs(outputpath)
        scriptfile = open(outputpath+"/"+basename+".sh",'w')

        i = 0		
        for tree in treefile:
            
            # Save tree to temporary file
            tree_path = os.path.join(outputpath, f"tmp_tree_{i}.nwk")
            with open(tree_path, "w") as tf:
                tf.write(tree.strip() + "\n")

            # Run Seq-Gen to generate alignment
            aln_path = os.path.join(outputpath, f"aln_{i}.fasta")
            cmd = f"{seqGenPath} -mHKY -t0.5 -f0.25,0.25,0.25,0.25 -l{pars['seqLength']} -s{pars['clockRate']} -n1 < {tree_path} > {aln_path}"
            os.system(cmd)

            # Parse alignment and add to parameters
            pars_alignment(aln_path, pars)

            pars["tree"]  = tree
            pars_dates(tree, pars)
            

            # Create XML file			
            pars["name"] = basename+".T"+str(i)			
            makeXMLFile(pars, template, outputfile=pars["name"], outputpath=outputpath)

            # Write command to scripts
            for seed in seeds:
                cmd ="java -jar -Xms2G -Xmx4G $1 -overwrite -seed %d %s" % (seed, pars["name"]+".xml") 
                scriptfile.write("%s\n" % (cmd))

            i += 1

            # Append alignment to combined file with tree-specific headers
            combined_aln_path = os.path.join(outputpath, "combined_alignments.fasta")
            with open(aln_path, "r") as aln_in, open(combined_aln_path, "a") as aln_out:
                for line in aln_in:
                    if line.startswith(">"):
                        line = f">T{i}_{line[1:]}"  # Prefix with tree ID
                    aln_out.write(line)

            # Delete temporary files
            os.remove(aln_path)
            os.remove(tree_path)
            
            
        treefile.close()		
        scriptfile.close()
    #
    #
