import os
import subprocess

from modules.log import logger
from modules.detected_CNV import Detected_CNV
from modules.CNV_detection_algorithms.CNV_algorithm import CNV_Algorithm

class Decon(CNV_Algorithm):

    def __init__(self, docker_conf, reference_conf, Bed, Run_class, force_run=False):
        super().__init__(docker_conf, reference_conf, Bed, force_run)
        self.decon_image = docker_conf.decon["image"]
        self.decon_version = docker_conf.decon["version"]

        self.Run = Run_class
        self.decon_folder = os.path.join(self.main_dir, "DECON")
        if not os.path.exists(self.decon_folder):
            os.mkdir(self.decon_folder)
        self.bams_dir = self.Run.run_path
        self.run_id = self.Run.run_id
        self.input_filename = f"{self.run_id}_decon_input.txt"
        self.input_path = os.path.join(self.decon_folder, self.input_filename)
        self.bam_volume = "bam_vol"
        self.decon_results_dir = os.path.join(self.results_dir, "DECON")
        if not os.path.exists(self.decon_results_dir):
            os.mkdir(self.decon_results_dir)
        
    
    def get_input_file(self):
        with open(self.input_path, "w") as f:
            for sample in self.Run.samples_147:
                # take of outliers
                if sample.is_outlier is not False:
                    continue
                sample_path = os.path.join(
                    "/bam_vol",
                    os.path.basename(sample.bam.path)
                )

                f.write(f"{sample_path}\n")
    
    def run_read_in_bams(self):
        output_path = os.path.join(self.decon_folder, f"{self.run_id}.RData")
        if os.path.exists(output_path):
            logger.info(
                f"Output of Decon ReadInBams for runID: {self.run_id} already exists: {output_path}"
            )
            return(True)
        fasta_dir = os.path.dirname(self.reference_fasta)
        fasta_filename = os.path.basename(self.reference_fasta)

        cmd = [
            self.docker_path, "run",
            "-v", f"{self.decon_folder}:/decon_folder",
            "-v", f"{self.bams_dir}:/bam_vol",
            "-v", f"{self.Bed.dir}:{self.Bed.volume}",
            "-v", f"{fasta_dir}:/fasta_dir",
            f"{self.decon_image}:{self.decon_version}",
            "ReadInBams.R", 
            "--bams", f"/decon_folder/{self.input_filename}",
            "--bed", f"{self.Bed.volume}/{self.Bed.sorted_merged_bed_filename}",
            "--fasta", f"/fasta_dir/{fasta_filename}",
            "--out", f"/decon_folder/{self.run_id}"
        ]

        self.run_cmd(cmd, "DECON ReadInBams")
        
    def run_identify_failed_rois(self):

        output_path = os.path.join(self.decon_folder, f"{self.run_id}_failed_rois_Failures.txt")
        if os.path.exists(output_path):
            logger.info(
                f"Output of Decon IdentifyFailures for runID: {self.run_id} already exists: {output_path}"
            )
            return True
        cmd = [
            self.docker_path, "run",
            "-v", f"{self.decon_folder}:/decon_folder",
            f"{self.decon_image}:{self.decon_version}",
            "IdentifyFailures.R", 
            "--Rdata", f"/decon_folder/{self.run_id}.RData",
            "--mincorr", "0.98",
            "--mincov", "100",
            "--out", f"/decon_folder/{self.run_id}_failed_rois"
        ]

        self.run_cmd(cmd, "DECON IdentifyFailures")

    def run_make_CNVcalls(self):
        output_path = os.path.join(self.decon_results_dir, f"decon_{self.run_id}_cnvs.RData")
        if os.path.exists(output_path):
            logger.info(
                f"Output of Decon makeCNVcalls for runID: {self.run_id} already exists: {output_path}"
            )
            return True
        cmd = [
            self.docker_path, "run",
            "-v", f"{self.decon_folder}:/decon_folder",
            "-v", f"{self.decon_results_dir}:/decon_results",
            f"{self.decon_image}:{self.decon_version}",
            "makeCNVcalls.R",
            "--Rdata", f"/decon_folder/{self.run_id}.RData",
            "--transprob", "0.01",
            "--out", f"/decon_results/decon_{self.run_id}_cnvs"
        ]

        self.run_cmd(cmd, "DECON makeCNVcalls")
    

    def parse_DECON_result(self, sample_id_sample_obj: dict =None):
        output_files = os.listdir(self.decon_results_dir)
        for file in output_files:
            if not file.endswith("all.txt"):
                continue
            output_path = os.path.join(self.decon_results_dir, file)
            with open(output_path, "r") as f:
                for i, line in enumerate(f):
                    if i == 0:
                        continue
                    line = line.strip().split("\t")
                    # print(line[1], "hey")
                    sample =  line[1].split(".")[0]
                    cnv_type, n_exons, start, end, chr, bayes_factor, gene = line[6], line[7], line[8], line[9], line[10], line[12], line[16]
                    if cnv_type == "deletion":
                        cnv_type = "DEL"
                    elif cnv_type == "duplication":
                        cnv_type = "DUP"
                    else:
                        raise ValueError(
                            f"CNV type detected by DECON is neither a duplication or deletion. Type is : {cnv_type}"
                        )

                    cnv = Detected_CNV(start, end, chr, cnv_type, sample, n_exons, gene, algorithm="DECON", qual=bayes_factor)
                    if sample_id_sample_obj:
                        # print(Sample_class.sample_id_sample_obj)
                        sample_obj = sample_id_sample_obj[sample]
                        sample_obj.cnvs["decon"].append(cnv)


