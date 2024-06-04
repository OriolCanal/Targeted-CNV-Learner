import os
import re

class CNV():
    # all_cnvs = list()
    def __init__(self, start, end, chr, type, sample, algorithm=None, qual=None):
        self.start = start
        self.end = end
        self.chr = chr
        self.type = type
        self.algorithm = algorithm
        self.qual = qual
        self.sample = sample
        # CNV.all_cnvs.append(self)


    def get_vcf_line(self):
        line = f"{self.chr}\t{self.start}\tN\t{self.type}\t{self.qual}\t.\tEND={self.end},SAMPLE={self.sample},ALGORITHM={self.algorithm}\t.\t.\n"
        return line
    
    @classmethod
    def create_cnvs_vcf(cls, vcf_path):
        with open(vcf_path, "w") as f:
            for cnv in cls.all_cnvs:
                vcf_line = cnv.get_vcf_line()
                f.write(vcf_line)

class Detected_CNV(CNV):
    all_detected_cnvs = list()
    def __init__(self, start, end, chr, type, sample, numb_exons, gene, algorithm="in_silico", qual=None):
        super().__init__(start, end, chr, type, sample, algorithm, qual)
        Detected_CNV.all_detected_cnvs.append(self)


class In_Silico_CNV(CNV):
    all_in_silico_cnvs = list()
    def __init__(self, start, end, chr, type, sample, numb_exons, gene, algorithm="in_silico", qual=None):
        super().__init__(start, end, chr, type, sample, algorithm, qual)
        self.numb_exons = numb_exons
        In_Silico_CNV.all_in_silico_cnvs.append(self)
        self.gene = gene

    @staticmethod
    def parse_cnvs_config(config_path, sample_id_sample_obj=None):
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

