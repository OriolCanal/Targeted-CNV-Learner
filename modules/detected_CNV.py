import os
import re
import pybedtools
import matplotlib.pyplot as plt
from Bio import SeqIO
import pybedtools

class CNV():
    all_cnvs = list()
    def __init__(self, start, end, chr, type, sample, algorithm=None, qual=None):
        self.start = int(start)
        self.end = int(end)
        self.chr = chr
        self.type = type
        self.algorithm = algorithm
        self.qual = qual
        self.sample = sample
        CNV.all_cnvs.append(self)


    def get_vcf_line(self):
        line = f"{self.chr}\t{self.start}\tN\t{self.type}\t{self.qual}\t.\tEND={self.end},SAMPLE={self.sample},ALGORITHM={self.algorithm}\t.\t.\n"
        return line
    
    @classmethod
    def create_cnvs_vcf(cls, vcf_path):
        with open(vcf_path, "w") as f:
            for cnv in cls.all_cnvs:
                vcf_line = cnv.get_vcf_line()
                f.write(vcf_line)

    def get_cnv_length(self):
        return(self.end - self.start)

    def calculate_overlap_percentage(self, cnv2):
        overlap_start = max(self.start, cnv2.start)
        overlap_end = min(self.end, cnv2.end)

        overlap_length = max(0, overlap_end - overlap_start)

        length_cnv1 = self.get_cnv_length()
        length_cnv2 = cnv2.get_cnv_length()

        overlap_percentage1 = (overlap_length / length_cnv1) * 100 if length_cnv1 != 0 else 0
        overlap_percentage2 = (overlap_length / length_cnv2) * 100 if length_cnv2 != 0 else 0

        return(overlap_percentage1, overlap_percentage2) 

    def is_overlap_acceptable(self, cnv2, threshold=30):
        """
        the overlap will be considered acceptable only if both CNVs overlap with 
        each other by at least the threshold percentage. This means that both the
        in silico CNV and the algorithm-detected CNV must overlap each other by the specified percentage.
        """
        overlap_percentage1, overlap_percentage2 = self.calculate_overlap_percentage(cnv2)
        return overlap_percentage1 >= threshold and overlap_percentage2 >= threshold
    
class Detected_CNV(CNV):
    cnvs = dict()
    # % of grapes CNV overlapping a in silico CNV
    grapes_overlaps = list()
    def __init__(self, start, end, chr, type, sample, numb_exons, gene, algorithm, qual=None):
        super().__init__(start, end, chr, type, sample, algorithm, qual)
        self.numb_exons = numb_exons
        self.gene = gene
        self.algorithm = algorithm.upper()
        if self.algorithm not in Detected_CNV.cnvs:
            Detected_CNV.cnvs[self.algorithm] = list()
        Detected_CNV.cnvs[self.algorithm].append(self)
        self.true_positive_cnv = False
        self.overlap_percentage = 0
        self.sequence = None

    
    def __str__(self):
        return(f"CNV detected by {self.algorithm} in sample {self.sample}, {self.chr}:{self.start}-{self.end}")

    def __repr__(self) -> str:
        return(f"CNV detected by {self.algorithm} in sample {self.sample}, {self.chr}:{self.start}-{self.end}")
    
    def set_in_silico_cnv(self, in_silico_cnv):
        if self.sample != in_silico_cnv.sample:
            raise ValueError(
                f"Samples for in silico and detected cnv are not the same"
            )
        self.true_positive_cnv = True
        self.in_silico_cnv = in_silico_cnv

    def extract_sequence(self, fasta_file):
        """
        Extract a sequence from a FASTA file
        """

        for record in SeqIO.parse(fasta_file, "fasta"):
            if record.id == self.chr:
                self.sequence = record.seq[self.start-1:self.end]
                return(str(self.sequence))
            raise ValueError(
                f"Chromosome {self.chr} not found in {fasta_file}"
            )
    
    def get_gc_content(self):
        if not self.sequence:
            self.extract_sequence()
        g_count = self.sequence.count("G")
        c_count = self.sequence.count("C")

        self.gc_content = float(g_count + c_count) / len(self.sequence) * 100
        return(self.gc_content)


    def get_mappability(self, mappability_file):
        region = f"{self.chr}:{self.start}-{self.end}"
        bedtool = pybedtools.BedTool(mappability_file)
        mappability_scores = list()

        for interval in bedtool.tabix_intervals(region):
            mappability_scores.append(float(interval[3]))

        if mappability_scores:
            self.mappability = sum(mappability_scores) / len(mappability_scores)
            return self.mappability
        else:
            raise ValueError(f"Mappability score could not be obtained for CNV: {self}")





    def get_numb_exons(self, Bed):

        exon_bed = pybedtools.BedTool(Bed.path)

        range_bed = pybedtools.BedTool(f"{self.chr}\t{self.start}\t{self.end}", from_string=True)

        overlapping_exons = exon_bed.intersect(range_bed)

        self.numb_exons = len(overlapping_exons)
        return(self.numb_exons)


    def get_gene_name(self, Bed):
        exon_bed = pybedtools.BedTool(Bed.path)

        range_bed = pybedtools.BedTool(f"{self.chr}\t{self.start}\t{self.end}", from_string=True)

        intersected = exon_bed.intersect(range_bed)

        gene_name = [interval[3] for interval in intersected]

        self.gene_name = gene_name[0]
        return(self.gene_name)


class In_Silico_CNV(CNV):
    all_in_silico_cnvs = list()
    # % of insilico cnv covered by detected grapes cnv
    grapes_overlaps = list()
    def __init__(self, start, end, chr, type, sample, numb_exons, gene, algorithm="in_silico", qual=None):
        super().__init__(start, end, chr, type, sample, algorithm, qual)
        self.numb_exons = numb_exons
        In_Silico_CNV.all_in_silico_cnvs.append(self)
        self.gene = gene
        # store algorihtms that have detected this cnv
        self.algorithms_detected = set()
        # stores algorithms as keys and % of overlapping as value
        self.algorithms_overlap = {
            "cnvkit": 0,
            "grapes2": 0,
            "decon": 0,
            "gatk": 0
        }
    

    @staticmethod
    def parse_cnvs_config(config_path, sample_id_sample_obj: dict =None):
        with open(config_path, "r") as f:
            for line in f:
                line = line.strip().split("\t")
                bam_path = line[0]
                sample = os.path.basename(bam_path).split(".")[0]
                chr, start, end, cnv_type, gene_exons, single_multipe_exons = line[1], line[2], line[3], line[4], line[5], line[6]
                if "_" not in gene_exons:
                    gene = gene_exons
                    exons = 1
                else:
                    match = re.search(r"(\w+)_(\d+)_to_(\d+)", gene_exons)  # Updated pattern
                    if match:
                        gene = match.group(1)
                        exon_start = int(match.group(2))
                        exon_end = int(match.group(3))
                        exons = exon_end - exon_start + 1

                cnv = In_Silico_CNV(start, end, chr, cnv_type, sample, exons, gene)
                if sample_id_sample_obj:
                    sample_obj = sample_id_sample_obj[sample]
                    sample_obj.cnvs["in_silico"].append(cnv)
        return(In_Silico_CNV.all_in_silico_cnvs)
    
    def plot_overlap_histograms(self, algorithm_overlaps, algorithm_name="Algorithm", ref_conf=None):

        plt.figure(figsize=(10, 6))

        plt.hist(self.in_silico_overlaps, bins=20, alpha=0.5, label='In Silico Overlap')
        plt.hist(algorithm_overlaps, bins=20, alpha=0.5, label=f'{algorithm_name} Overlap')

        plt.xlabel('Overlap Percentage')
        plt.ylabel('Frequency')
        plt.title('Histogram of Overlap Percentages')
        plt.legend(loc='upper right')
        plt.grid(True)
        if ref_conf:
            plots_dir =  os.path.join(ref_conf.main_dir, "Plots", "overlap_cnv")
            if not os.path.exists(plots_dir):
                os.mkdir(plots_dir)
            plot_path = os.path.join(plots_dir, f"In_silico_vs_{algorithm_name}_overlap.png")
            plt.savefig(plot_path)
        plt.show()
