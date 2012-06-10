synbioinformatica
=================

Command line interface for computable molecular / synthetic biology construction files


##Example of new construction file, simple ecobam PCR
##this should be all read in by a wrapper file that loads the appropriate modules
##ca998, ca1034F and pSB1A2-Bca9033 should be auto-loaded by the db
## we should also not be strict on upper/lowercase... not sure how to do that 

prod = pcr(ca998, ca1034F, pSB1A2-Bca9033)
gp = gel_purify(prod)
ins = digest( [ecori, spei], zymo)

vector = gel_digest([ecori, spei], pbca1020-bca9022, large)

pbca1020-bth2039 = ligate_transform( [vector, ins], strain(MC1061) )

new_bglb_part ( pbca1020-bth2039 )  #this returns the dna that's inbetween the bglii and bam




##determines optimal set (1-2) of enzymes to use to differentiate between a list of plasmids
enzyme_set = find_differentiating_enzymes( list_of_dnas )
## want to find a restriction digest for the dnas that results in a maximal difference
## needs to be at least one dna band that is differentiable... hrmm
## the differentiable dna band should be 1.5x up/down of any dna bands in the other lane
## the differentiable dna band should not be < 300bp
##


##These may eventually be used to write a generalized function that determines which procedure to do


#####################
##SOEing example
#####################
prod1 = pcr(TH123, TH125, pth7035k-Bth8055)
prod2 = pcr(TH343, TH255, pth7035K-Bth8099)
final = SOE(TH123, TH255, [prod1, prod2])

####################
##primer design, designs the primers for making new_dna from pth7035K-Bth8199
####################
primer_list = design_eipcr(pth7035K-Bth8199, new_dna)
##identifies the region that is different
##checks if there is <2 BsaI or BseRI sites present, if not, chooses a "good" restriction enzyme and throws warning
##        --NNN----->   -----> needs to have Tm of 55C, -- tail needs to be long enough for enzyme site


primer_list = design_SOE(old_dna, new_dna)
##identifies two good restriction sites (defaults to ecobam first) that flank the mutation, finds up and downstream 
##  primers, and designs primers for the mutation




