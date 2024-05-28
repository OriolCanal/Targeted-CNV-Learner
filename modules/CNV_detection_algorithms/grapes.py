import os
import json

from modules.log import logger
from modules.CNV_detection_algorithms.CNV_algorithm import CNV_Algorithm


class Grapes(CNV_Algorithm):
    def __init__(self, docker_conf, reference_conf, Bed, force_run=False):
        super().__init__(docker_conf, reference_conf, Bed, force_run)
        self.grapes_image = docker_conf.grapes["image"]
        self.grapes_version = docker_conf.grapes["version"]
        self.grapes_dir = os.path.join(self.main_dir, "Grapes")
        self.grapes_results_dir = os.path.join(self.results_dir, "Grapes")
        if not os.path.exists(self.grapes_results_dir):
            os.mkdir(self.grapes_results_dir)
        
        
    def get_grapes_bed(self):
        bed_path = self.Bed.path

        grapes_bed_path = bed_path.replace(".bed", ".grapes.bed")
        if os.path.exists(grapes_bed_path):
            self.Bed.set_grapes_bed(grapes_bed)
            logger.info(
                f"Grapes bed already exists: {grapes_bed_path}"
            )
            return True

        logger.info(
            f"Creating Grapes bed in {grapes_bed_path}"
        )
        with open(bed_path, "r") as bed:
            with open(grapes_bed_path, "w") as grapes_bed:
                for line in bed:
                    line = line.rstrip("\n")
                    tmp = line.split("\t")
                    if len(tmp) > 3:
                        info = tmp[3].replace("'", '"')
                        try:
                            info_dict = json.loads(info)
                        except:
                            info = "ROI"
                        else:
                            gene_name = "."
                            if "gene_name" in info_dict:
                                gene_name = info_dict["gene_name"]
                            isoform = "."
                            if "Parent" in info_dict:
                                isoform = info_dict["Parent"]
                            exon = "."
                            if "exon_number" in info_dict:
                                exon = info_dict["exon_number"]
                            exon_str = str(int(exon) - 1) + "_" + exon
                            info = f"{isoform}_{exon_str};{gene_name}"
                        tmp[3] = info
                        line = "\t".join(tmp)
                        grapes_bed.write(line + "\n")
        
        self.Bed.set_grapes_bed(grapes_bed)
    
    def run_grapes(self, control_samples, analysis_samples):
        cohort_symbolic_link = control_samples[0].bam.symbolic_link
        analysis_symbolic_link = analysis_samples[0].bam.symbolic_link
        cmd = [
            self.docker_path, "run",
            "-v", f"{self.fasta_dir}:/fasta_dir",
            "-v", f"{self.Bed.dir}:{self.Bed.volume}",
            "-v", f"{cohort_symbolic_link}:/cohort_symbolic_link",
            "-v", f"{analysis_symbolic_link}:/analysis_symbolic_link",
            "-v", f"{self.grapes_results_dir}:/grapes_results",
            f"{self.grapes_image}:{self.grapes_version}", "grapes2.py",
            "--cases", "/analysis_symbolic_link/", 
            "--control", "/cohort_symbolic_link/",
            "--output", "/grapes_results",
            "--bed", f"{self.Bed.volume}/{self.Bed.grapes_bed_filename}",
            "--fasta", f"/fasta_dir/{self.fasta_filename}",
            "--genome", "hg19",
            "--all"
        ]

        self.run_cmd(cmd, "Running Grapes:")