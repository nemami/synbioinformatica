#!/usr/bin/python -tt
# Copyright Nima Emami, 2012

# PCR function tester script
# Takes 3 invocation arguments: 
# template, forward primer, reverse primer (in that order)
# returns tuple: (PCR product, start index, end index)

import synbioinformatica, sys

def main(template, fwd_primer, rev_primer):
	print PCR(template, fwd_primer, rev_primer)

if __name__ == "__main__":
	PCR = getattr(synbioinformatica,'PCR')
	main(sys.argv[1],sys.argv[2],sys.argv[3])