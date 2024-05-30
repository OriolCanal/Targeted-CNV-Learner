

class CNV():
    all_cnvs = list()
    def __init__(self, start, end, chr, type, algorithm, qual, sample):
        self.start = start
        self.end = end
        self.chr = chr
        self.type = type
        self.algorithm = algorithm
        self.qual = qual
        self.sample = sample
        CNV.all_cnvs.append(self)


    def get_vcf_line(self):
        line = f"{self.chr}\t{self.start}\tN\t{self.type}\t{self.qual}\t.\tEND={self.end},SAMPLE={self.sample},ALGORITHM={self.algorithm}\t.\t.\n"
        return line
    
    def create_cnvs_vcf(self, vcf_path):
        with open(vcf_path, "w") as f:
            for cnv in CNV.all_cnvs:
                vcf_line = cnv.get_vcf_line()
                f.write(vcf_line)
        