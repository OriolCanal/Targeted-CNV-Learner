import os
import re
import pybedtools
import matplotlib.pyplot as plt
from Bio import SeqIO
import pybedtools
from pyfaidx import Fasta
import pandas as pd
import pyBigWig
from modules.log import logger
from pybedtools import BedTool


class CNV():
    all_cnvs = list()
    def __init__(self, start, end, chr, type, sample, algorithm=None, qual=None):
        self.start = int(start)
        self.end = int(end)
        self.chr = chr
        if type.upper() not in ["DUP", "DEL"]:
            raise ValueError(
                f"CNV type should be either DEL or DUP and {type} was specified"
            )
        self.type = type
        self.algorithm = algorithm
        self.qual = qual
        self.sample = sample
        CNV.all_cnvs.append(self)
    def __str__(self):
        return(f"CNV in sample {self.sample}, {self.type}: {self.chr}:{self.start}-{self.end}")

    def __repr__(self) -> str:
        return(f"CNV in sample {self.sample}, {self.type}: {self.chr}:{self.start}-{self.end}")


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
        if self.chr != cnv2.chr:
            return 0

        if self.type != cnv2.type:
            return 0

        overlap_start = max(self.start, cnv2.start)
        overlap_end = min(self.end, cnv2.end)

        overlap_length = max(0, overlap_end - overlap_start)

        length_cnv1 = self.get_cnv_length()
        length_cnv2 = cnv2.get_cnv_length()

        overlap_percentage1 = (overlap_length / length_cnv1) * 100 if length_cnv1 != 0 else 0
        overlap_percentage2 = (overlap_length / length_cnv2) * 100 if length_cnv2 != 0 else 0
        # print(overlap_percentage1, overlap_percentage2, "\n\n\noverlap_percentages")
        return(overlap_percentage1, overlap_percentage2) 

    def is_overlap_acceptable(self, cnv2, threshold=30):
        """
        the overlap will be considered acceptable only if both CNVs overlap with 
        each other by at least the threshold percentage. This means that both the
        in silico CNV and the algorithm-detected CNV must overlap each other by the specified percentage.
        """
        result = self.calculate_overlap_percentage(cnv2)
        logger.info(f"overlap percentage: {result}")
        if result == 0:
            return(False)
        else:
            overlap_percentage1, overlap_percentage2 = result
            return overlap_percentage1 >= threshold and overlap_percentage2 >= threshold

    def extract_all_sequence(self, fasta_file):
        """
        Extract a sequence from a FASTA file
        """
        fasta = Fasta(fasta_file)
        self.sequence = str(fasta[self.chr][self.start-1:self.end]).upper()

        if not self.sequence:
            ValueError(f"Sequence not foud for {self} in {fasta_file}")

        fasta.close()

        return self.sequence


    def extract_bed_sequence(self, fasta_file, bed_file):
        bed = BedTool(bed_file)
        bed.sequence(fi=fasta_file, fo="")
    def get_gc_content(self, fasta_file):
        if not self.sequence:
            self.extract_all_sequence(fasta_file)
        g_count = self.sequence.count("G")
        c_count = self.sequence.count("C")
        a_count = self.sequence.count("A")
        t_count = self.sequence.count("T")
        self.at_content = float(a_count + t_count) / len(self.sequence) * 100
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

        exon_bed = pybedtools.BedTool(Bed.roi_bed)

        range_bed = pybedtools.BedTool(f"{self.chr}\t{self.start}\t{self.end}", from_string=True)

        overlapping_exons = exon_bed.intersect(range_bed)

        self.numb_exons = len(overlapping_exons)
        return(self.numb_exons)


    def get_gene_name(self, Bed):
        exon_bed = pybedtools.BedTool(Bed.roi_bed)

        range_bed = pybedtools.BedTool(f"{self.chr}\t{self.start}\t{self.end}", from_string=True)

        intersected = exon_bed.intersect(range_bed)
        # print(intersected, "intesected")
        if not intersected:
            self.gene_name = None
            return(self.gene_name)
        gene_name = [interval[3] for interval in intersected]

        self.gene = gene_name[0]
        return(self.gene)
    
    def get_exon_number_in_gene(self, Bed):
        if not hasattr(Bed, "sorted_roi"):
            Bed.sort_roi_bed()
        exon_bed = pybedtools.BedTool(Bed.sorted_roi)

class Detected_CNV(CNV):
    cnvs = dict()
    sample_cnvs = dict()
    # % of grapes CNV overlapping a in silico CNV
    grapes_overlaps = list()
    def __init__(self, start, end, chr, type, sample, numb_exons, gene, algorithm, qual=None):
        super().__init__(start, end, chr, type, sample, algorithm, qual)
        self.numb_exons = numb_exons
        self.gene = gene
        
        self.true_positive_cnv = False
        self.overlap_percentage = 0
        self.sequence = None
        self.gc_content = None
        self.at_content = None
        self.algorithm = algorithm.upper()
        if self.algorithm not in Detected_CNV.cnvs:
            Detected_CNV.cnvs[self.algorithm] = list()
        Detected_CNV.cnvs[self.algorithm].append(self)

        if self.sample not in Detected_CNV.sample_cnvs:
            Detected_CNV.sample_cnvs[self.sample] = list()
        Detected_CNV.sample_cnvs[self.sample].append(self)

    def __str__(self) -> str:
        return(f"CNV detected by {self.algorithm} in sample {self.sample}, {self.type}: {self.chr}:{self.start}-{self.end}")

    def __repr__(self) -> str:
        return(f"CNV detected by {self.algorithm} in sample {self.sample}, {self.type}: {self.chr}:{self.start}-{self.end}")


    @classmethod
    def get_overlapping_cnvs(cls):
        overlap_results = list()
        for sample, sample_cnvs in cls.sample_cnvs.items():
            for i in range(len(sample_cnvs)):
                for j in range(i+1, len(sample_cnvs)):
                    cnv1 = sample_cnvs[i]
                    cnv2 = sample_cnvs[j]
                    if cnv1.algorithm != cnv2.algorithm:
                        overlap_percentage = cls.calculate_overlap(cnv1, cnv2)
                        if overlap_percentage > 0:
                            overlap_results.append((cnv1, cnv2, overlap_percentage))
        return overlap_results
            # sample_cnvs = Detected_CNV.sample_cnvs[sample]
            # for cnv1 in range(len(sample_cnvs)):
            #     for cnv2 in range(i+1, len(sample_cnvs)):
            #         cnv1 = 

    def set_in_silico_cnv(self, in_silico_cnv):
        if self.sample != in_silico_cnv.sample:
            raise ValueError(
                f"Samples for in silico and detected cnv are not the same"
            )
        self.true_positive_cnv = True
        self.in_silico_cnv = in_silico_cnv

class Overlapping_CNV(CNV):
    available_algs = ["GRAPES2", "CNVKIT", "DECON", "GATK_GCNV"]
    sample_cnvs = dict()
    all_joined_cnvs = list()
    def __init__(self, start, end, chr, type, sample, algorithm, numb_exons=None, gene=None, qual=None):
        super().__init__(start, end, chr, type, sample, algorithm, qual)
        self.numb_exons = numb_exons
        self.gene = gene
        
        self.true_positive_cnv = False
        self.overlap_percentage = 0
        self.sequence = None
        self.gc_content = None
        self.at_content = None
        self.algorihtm = algorithm
        self.qual = qual
        self.origin = None
        for algorithm in self.algorihtm:
            if algorithm.upper() not in Overlapping_CNV.available_algs:
                raise ValueError(
                    f"Algorithm {algorithm} not in available algorihms: {Overlapping_CNV.available_algs}"
                )

        if self.sample not in Overlapping_CNV.sample_cnvs:
            Overlapping_CNV.sample_cnvs[self.sample] = list()
        Overlapping_CNV.sample_cnvs[self.sample].append(self)
        # self.sample.joined_cnv.append(self)
        Overlapping_CNV.all_joined_cnvs.append(self)
    def __str__(self) -> str:
        return(f"CNV detected by {self.algorithm} in sample {self.sample}, {self.type}: {self.chr}:{self.start}-{self.end}")

    def __repr__(self) -> str:
        return(f"CNV detected by {self.algorithm} in sample {self.sample}, {self.type}: {self.chr}:{self.start}-{self.end}")
    
    @classmethod
    def get_df_all_cnvs(cls, Joined_Mosdepth_df=None):
        data = list()
        for cnv in cls.all_joined_cnvs:
            # if isinstance(cnv, In_Silico_CNV):
            #     cnv_origin = "in_silico"
            # elif isinstance(cnv, Real_CNV):
            #     cnv_origin = "real_cnv"
            # elif isinstance(cnv, Overlapping_CNV):
            #     cnv_origin = "overlapping"
            # else:
            #     raise ValueError(
            #         f"CNV: {cnv} is not an insilico or real CNV"
            #     )
            
            if "DECON" in cnv.algorithm:
                decon = 1
            else:
                decon = 0

            if "GRAPES2" in cnv.algorithm:
                grapes = 1
            else:
                grapes = 0
            
            if "CNVKIT" in cnv.algorithm:
                cnvkit = 1
            else:
                cnvkit = 0
            
            if "GATK_GCNV" in cnv.algorithm:
                gatk = 1
            else:
                gatk = 0
            cnv_dict = {
                "start": cnv.start,
                "end": cnv.end,
                "chr": cnv.chr,
                "type": cnv.type,
                "sample": cnv.sample.sample_id,
                "numb_exons": cnv.numb_exons,
                "decon": decon,
                "gatk": gatk,
                "cnvkit": cnvkit,
                "grapes": grapes,
                "gc_content": cnv.gc_content,
                "mappability": cnv.mappability,
                "cnv_length": cnv.get_cnv_length(),
                "gene": cnv.gene,
                "qual": cnv.qual,
                "true_positive": cnv.true_positive_cnv,
                "decon_qual": cnv.qual["DECON"],
                "gatk_qual": cnv.qual["GATK_GCNV"],
                "cnvkit_qual": cnv.qual["CNVKIT"],
                "grapes_qual": cnv.qual["GRAPES2"],
                "cnv_origin": cnv.origin

            }
            if Joined_Mosdepth_df:
                sample_corr = Joined_Mosdepth_df.get_sample_correlation(cnv.sample.sample_id)
                cnv_dict["sample_correlation"] = sample_corr
                cnv_dict["is_outlier"] = cnv.sample.is_outlier

            data.append(cnv_dict)
        df = pd.DataFrame(data)
        return(df)


    def set_in_silico_cnv(self, in_silico_cnv):
        if self.sample.sample_id != in_silico_cnv.sample:
            raise ValueError(
                f"Samples for in silico and detected cnv are not the same"
            )
        self.true_positive_cnv = True
        self.in_silico_cnv = in_silico_cnv

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



class Real_CNV(CNV):
    all_real_cnvs = list()

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

    @classmethod
    def get_detected_cnvs(cls, detected_cnvs_csv, sample_id_sample_obj):
        cnvs_df = pd.read_csv(detected_cnvs_csv)
        all_real_cnvs = list()
        for index, row in cnvs_df.iterrows():
            sample_id = row.LAB_ID
            if sample_id not in sample_id_sample_obj:
                print(sample_id_sample_obj)
                logger.warning(
                    f"Sample id: {sample_id} specified in detected CNVS csv file: {detected_cnvs_csv} not match with any sample_id"
                )
                continue

            pos = row.hg19
            chr, start_end = pos.split(":")
            chr = chr.replace("chr", "")
            start, end = map(int, start_end.split("-"))

            cnv_type = row.TYPE
            gene = row.GENE
            numb_exons = None
            algorithm = "real_cnv"

            real_cnv = cls(start, end, chr, cnv_type, sample_id, numb_exons, gene, algorithm)
            sample_obj = sample_id_sample_obj[sample_id]
            logger.info(
                f"sample in {sample_id} found, {sample_obj}, assigning real cnv to this cnv"
            )
            sample_obj.cnvs["real_cnv"].append(real_cnv)
            logger.info(
                f"real CNV: {real_cnv} parsed correctly"
            )
            all_real_cnvs.append(real_cnv)
        return(all_real_cnvs)