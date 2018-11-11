#!/bin/python3
# NOTES
# Can't use * to unpack collections because the notation is broken on 3.4 (aka the virtual machine test env)

# CONSTANTS
VERSION = 'v0.2.1'
RSCRIPT = './VCF_whisper_0_2_1.r'

# IMPORTS
# Standard libraries
import argparse
import csv
import subprocess
import re
import time

# Custom template
import Rmd_template

# CLASSES
# Well actually the only class, and it's more of a giant function TBH

class VCF:
    """
    Class for VCF analysis. 
    """

    def __init__(self, vcf_file):
        """
        TODO: init
        """
        self.vcf_file = vcf_file
        # stuff
        self.fileformat = ''
        self.reference = ''
        self.format_lines = 0
        self.info_lines = 0
        self.filters = []
        self.samples = []
        self.contigs = {}
        self.variants = 0
        self.contigs_tot = 0
        self.bin_size = 0
        # results
        self.variant_density = []
        self.variant_summary = []
        self.variant_overlap = []
        self.base_change = []
        self.ti_tv = []
        self.allele_frequency = {}
        self.allele_frequency_density = []

    def controller(self):
        """
        Main loop for the VCF reading program.
        Basically just calls all the appropriate functions in order.
        'Handles' what errors are even checked for.
        """
        # Read Header (n stuff)
        start_time = time.perf_counter()
        return_code = self.header_parse()
        if return_code == 1:
            print("File does not appear to be a VCF. Check if the header is formatted correctly.")
        header_time = time.perf_counter()
        print("Header time: {}\n".format(header_time - start_time))
        self.summary()
        summary_time = time.perf_counter()
        print("Summary time: {}\n".format(summary_time - header_time))
        self.printer()
        self.rmd_writer()
        print_time = time.perf_counter()
        print("File write time: {}\n".format(print_time - summary_time))
        self.plotter()
        plot_time = time.perf_counter()
        print("Rmarkdown time: {}\n".format(plot_time - print_time))

    def header_parse(self):
        """
        Reads the header and also figures out the contig lengths
        TODO: potentially figure out samples here?
        """
        with open(self.vcf_file, 'r') as f:
            line = f.readline()
            if line[:12] != '##fileformat':
                return 1
            while line[:2] == '##':
                if line[:13] == '##fileformat=':
                    self.fileformat = line[13:]
                if line[:12] == '##reference=':
                    self.reference = line[12:]
                if line[:8] == '##FORMAT':
                    self.format_lines += 1
                if line[:6] == '##INFO':
                    self.info_lines += 1
                if line[:8] == '##FILTER':
                    elements = line[line.find('<')+1:line.find('>')].split(',')
                    for e in elements:
                        if e[:3] == 'ID=':
                            ID = e[3:]
                        elif e[:12] == 'Description=':
                            desc = e[12:]
                    self.filters.append((ID, desc))
                if line[:8] == '##contig':
                    elements = line[line.find('<')+1:line.find('>')].split(',')
                    for e in elements:
                        if e[:3] == 'ID=':
                            ID = e[3:]
                        elif e[:7] == 'length=':
                            length = int(e[7:])
                    self.contigs[ID] = [length, 0]
                line = f.readline()
            # Skip over the column headers
            line = f.readline()
            # This only triggers if no information on contigs found in the header
            if len(self.contigs) == 0:
                current_contig = line.split('\t')[0]
                current_length = 0
                while line:
                    chrom, pos = line.split('\t')[:2]
                    if chrom != current_contig:
                        self.contigs[current_contig] = [int(current_length) + 10, 0]
                    current_contig = chrom
                    current_length = pos
                    line = f.readline()
                self.contigs[current_contig] = [int(current_length) + 10, 0]
        for c in self.contigs:
            self.contigs_tot += self.contigs[c][0]
        self.bin_size = self.contigs_tot//2000
        return 0

    def meta_inf_parse(self, ):
        """
        parse out arbirtary meta-inf lines
        """   
    
    def summary(self):
        """
        This is where all the work happens
        """
        with open(self.vcf_file, 'r') as f:
            # Jump over Header
            line = f.readline()
            while line[:2] == '##':
                line = f.readline()
            # Column names
            self.samples = [s.strip() for s in line.split('\t')[9:]]
            line = f.readline()
            # Body
            # here goes some variables to count bins
            bin_global = 1
            bin_local = 1
            current_contig = line.split('\t')[0] # first contig
            
            next_lines = [[bin_global, sample, current_contig, bin_local, 0, 0, 0, 0, 0] for sample in self.samples]
            afd_lines = [[0,0] for sample in self.samples]
            self.variant_summary = [[[0, 0], [0, 0], [0, 0], [0, 0]] for sample in self.samples]
            #[sample, snv, ins, ins_bp, del, del_bp, oth, oth_bp]
            self.variant_overlap = [0 for sample in self.samples]
            self.base_change = [[0 for sample in self.samples] for i in range(12)]
            self.ti_tv = [[0,0] for sample in self.samples]
            while line:
                self.variants += 1
                tab = line[:-1].split('\t')
                contig, pos, ID, ref, alt, qual, filt, info, form = tab[:9]
                if current_contig != contig: # New contig represents a new bin always and resets bin_local
                    self.variant_density.extend(next_lines)
                    self.allele_frequency_density.extend([[bin_global, sample, current_contig, bin_local, afd_lines[i][0]/(afd_lines[i][1] + 1)] for i, sample in enumerate(self.samples)])
                    bin_global += 1
                    bin_local = 1
                    current_contig = contig
                    next_lines = [[bin_global, sample, current_contig, bin_local, 0, 0, 0, 0, 0] for sample in self.samples]
                    afd_lines = [[0,0] for sample in self.samples]
                if int(pos) > (bin_local * self.bin_size): # this is not elif because the first bin can be empty
                    self.variant_density.extend(next_lines)
                    self.allele_frequency_density.extend([[bin_global, sample, current_contig, bin_local, afd_lines[i][0]/(afd_lines[i][1] + 1)] for i, sample in enumerate(self.samples)])
                    bin_global += 1
                    bin_local += 1
                    next_lines = [[bin_global, sample, current_contig, bin_local, 0, 0, 0, 0, 0] for sample in self.samples]
                    afd_lines = [[0,0] for sample in self.samples]
                for i, a in enumerate(alt.split(',')):
                    variant_present_count = 0
                    all_freq_tot = 0
                    all_freq_zero = 0
                    for n, s in enumerate(self.samples):
                        genotypes = re.split('[\|/]', tab[9+n][:tab[9+n].find(':')])
                        all_freq_tot += len(genotypes)
                        afd_lines[n][1] += len(genotypes)
                        all_freq_zero += genotypes.count('0')
                        afd_lines[n][0] += genotypes.count('0')
                        if str(i+1) in genotypes:
                            variant_present_count += 1
                            next_lines[n][-1] += 1 #running total
                            # TODO: test if this actually works
                            # replace this with one that looks only for pre-defined non-standard alts
                            # fukkin also - make this its own function again - AND 
                            if a[0] == '<': # non standard
                                if a[1:4] == 'INS':
                                    next_lines[n][5] += 1
                                elif a[1:4] == 'DEL':
                                    next_lines[n][6] += 1
                                else:
                                    next_lines[n][7] += 1
                            elif a == '.': # deletion
                                next_lines[n][6] += 1
                                self.variant_summary[n][2][0] += 1
                                self.variant_summary[n][2][1] += 1
                            elif len(ref) == 1 and len(a) == 1: # SNV
                                next_lines[n][4] += 1
                                self.variant_summary[n][0][0] += 1
                                self.variant_summary[n][0][1] += 1
                                self.base_changer(ref, a, n)
                            elif len(ref) == 1 and len(a) > 1: # INS
                                next_lines[n][5] += 1
                                self.variant_summary[n][1][0] += 1
                                self.variant_summary[n][1][1] += (len(a) - 1)
                            elif len(ref) > 1 and len(a) == 1: # DEL
                                next_lines[n][6] += 1
                                self.variant_summary[n][2][0] += 1
                                self.variant_summary[n][2][1] += (len(ref) - 1)
                            elif len(ref) > 1 and len(a) > 1: # other
                                next_lines[n][7] += 1
                                self.variant_summary[n][3][0] += 1
                                self.variant_summary[n][3][1] += abs(len(ref) - len(a))
                    self.variant_overlap[variant_present_count-1] += 1
                    all_freq = (all_freq_tot - all_freq_zero)/all_freq_tot
                    if all_freq in self.allele_frequency:
                        self.allele_frequency[all_freq] += 1
                    else:
                        self.allele_frequency[all_freq] = 1
                line = f.readline()

    def base_changer(self, ref, alt, sample_index):
        """
        Base change index
        [A->C,T,G C->A,T,G T->A,C,G G->A, C, T
            0,1,2    3,4,5    6,7,8    9,10,11
            1,1,0    1,0,1    1,0,1    0, 1, 1
        """
        if ref == 'A':
            if alt == 'C':
                self.base_change[0][sample_index] += 1
                self.ti_tv[sample_index][1] += 1
            elif alt =='T':
                self.base_change[1][sample_index] += 1
                self.ti_tv[sample_index][1] += 1
            elif alt == 'G':
                self.base_change[2][sample_index] += 1
                self.ti_tv[sample_index][0] += 1
        elif ref == 'C':
            if alt == 'A':
                self.base_change[3][sample_index] += 1
                self.ti_tv[sample_index][1] += 1
            if alt == 'T':
                self.base_change[4][sample_index] += 1
                self.ti_tv[sample_index][0] += 1
            if alt == 'G':
                self.base_change[5][sample_index] += 1
                self.ti_tv[sample_index][1] += 1
        elif ref == 'T':
            if alt == 'A':
                self.base_change[6][sample_index] += 1
                self.ti_tv[sample_index][1] += 1
            if alt == 'C':
                self.base_change[7][sample_index] += 1
                self.ti_tv[sample_index][0] += 1
            if alt == 'G':
                self.base_change[8][sample_index] += 1
                self.ti_tv[sample_index][1] += 1
        elif ref == 'G':
            if alt == 'A':
                self.base_change[9][sample_index] += 1
                self.ti_tv[sample_index][0] += 1
            if alt == 'C':
                self.base_change[10][sample_index] += 1
                self.ti_tv[sample_index][1] += 1
            if alt == 'T':
                self.base_change[11][sample_index] += 1
                self.ti_tv[sample_index][1] += 1

    def rmd_writer(self):
        """
        Oh geez...
        Actually this worked straight away
        """
        with open('VCF_whisper.Rmd', 'w') as f:
            f.write(Rmd_template.HEADER)
            f.write(Rmd_template.BASIC.format(len(self.samples), self.variants, len(self.contigs), self.contigs_tot))
            for c in self.contigs:
                f.write(c + ' | ' + str(self.contigs[c][0]) + '\n')
            f.write(Rmd_template.VAR_TYPE)
            f.write(Rmd_template.VAR_DENS.format(self.bin_size))
            if len(self.samples) > 1:
                f.write(Rmd_template.VAR_DENS_BY_SAMPLE_HEAD)
                for s in self.samples:
                    f.write(Rmd_template.VAR_DENS_BY_SAMPLE.format(s, self.bin_size))
            f.write(Rmd_template.VAR_OVERLAP)
            f.write(Rmd_template.ALL_FREQ)
            f.write(Rmd_template.AFD)
            if len(self.samples) > 1:
                f.write(Rmd_template.AFD_BY_SAMPLE_HEAD)
                for s in self.samples:
                    f.write(Rmd_template.AFD_BY_SAMPLE.format(s))
            f.write(Rmd_template.BASE_OVERLAP)

    def printer(self):
        """
        Prints results for input into R script
        """
        with open('whisper_results_variant_density.csv', 'w') as f:
            csvwriter = csv.writer(f, delimiter = ',', quotechar = '"')
            csvwriter.writerow(['Gbin', 'sample', 'contig', 'lbin', 'SNV', 'INS', 'DEL', 'OTH', 'TOT'])
            for n in self.variant_density:
                csvwriter.writerow(n)

        with open('whisper_results_variant_sum.csv', 'w') as f:
            csvwriter = csv.writer(f, delimiter = ',', quotechar = '"')
            csvwriter.writerow(['sample', 'type', 'count', 'bp'])
            variants = ['snv', 'ins', 'del', 'oth']
            for i, v in enumerate(variants):
                rows = []
                for j, s in enumerate(self.samples):
                    rows.append([s, v, self.variant_summary[j][i][0], self.variant_summary[j][i][1]])
                for r in rows:
                    csvwriter.writerow(r)

        with open('whisper_results_variant_overlap.csv', 'w') as f:
            csvwriter = csv.writer(f, delimiter = ',', quotechar = '"')
            csvwriter.writerow(['N', 'count'])
            for i,n in enumerate(self.variant_overlap):
                csvwriter.writerow([str(i+1), n])

        with open('whisper_results_allele_frequency.csv', 'w') as f:
            csvwriter = csv.writer(f, delimiter = ',', quotechar = '"')
            csvwriter.writerow(['zygosity_fraction', 'count'])
            jk = list(self.allele_frequency.keys())
            jk.sort()
            for i in jk:
                csvwriter.writerow([i, self.allele_frequency[i]])

        with open('whisper_results_allele_frequency_density.csv', 'w') as f:
            csvwriter = csv.writer(f, delimiter = ',', quotechar = '"')
            csvwriter.writerow(['Gbin', 'sample', 'contig', 'lbin', 'afd'])
            for n in self.allele_frequency_density:
                csvwriter.writerow(n)

        with open('whisper_results_ti_tv.csv', 'w') as f:
            csvwriter = csv.writer(f, delimiter = ',', quotechar = '"')
            csvwriter.writerow(['sample', 'type', 'count'])
            types = ['Ti', 'Tv']
            for i,v in enumerate(types):
                rows = []
                for j, s in enumerate(self.samples):
                    rows.append([s,v,self.ti_tv[j][i]])
                for r in rows:
                    csvwriter.writerow(r)

        with open('whisper_results_base_changes.csv', 'w') as f:
            csvwriter = csv.writer(f, delimiter = ',', quotechar = '"')
            row = ['sample', 'ref','alt', 'count']
            csvwriter.writerow(row)
            base_ref = ['A','A','A','C','C','C','T','T','T','G','G','G']
            base_alt = ['C','T','G','A','T','G','A','C','G','A','C','T']
            for i in range(len(base_ref)):
                rows = [[sample, base_ref[i], base_alt[i], 0] for sample in self.samples]
                for j in range(len(self.samples)):
                    rows[j][3] = self.base_change[i][j]
                for r in rows:
                    csvwriter.writerow(r)

    def plotter(self):
        """
        All this does is call an Rscipt
        """
        subprocess.call(["/usr/bin/Rscript", RSCRIPT])
        subprocess.run(["mv", "VCF_whisper.html", self.vcf_file + ".whisper.html"])

# FUNCTIONS

# MAINBLOCK

if __name__ == "__main__":
    cli_parser = argparse.ArgumentParser(description = 'VCF Whisper ' + VERSION)
    cli_parser.add_argument('filename', type=str, help = 'path to VCF file')
    args = cli_parser.parse_args()
    vcf1 = VCF(args.filename)
    vcf1.controller()
