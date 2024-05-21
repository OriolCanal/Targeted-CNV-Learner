import os
import subprocess

from modules.log import logger

from modules.CNV_detection_algorithms.CNV_algorithm import CNV_Algorithm

class Decon(CNV_Algorithm):

    def __init__(self, docker_conf, reference_conf, Bed, control_samples, force_run=False):
        super().__init__(docker_conf, reference_conf, Bed, force_run)
        self.decon_image = docker_conf.decon["image"]
        self.decon_version = docker_conf.decon["version"]

        self.control_samples = control_samples
        self.decon_folder = os.path.join(self.main_dir, "DECON")
        if not os.path.exists(self.decon_folder):
            os.mkdir(self.decon_folder)
        self.bams_dir = os.path.dirname(self.control_samples[0].bam.path)
        self.run_id = self.control_samples[0].run_id
        self.input_filename = f"{self.run_id}_decon_input.txt"
        self.input_path = os.path.join(self.decon_folder, self.input_filename)
        self.bam_volume = "bam_vol"
        self.decon_results_dir = os.path.join(self.results_dir, "DECON")
        if not os.path.exists(self.decon_results_dir):
            os.mkdir(self.decon_results_dir)
        
    
    def get_input_file(self):
        with open(self.input_path, "w") as f:
            for sample in self.control_samples:
                sample_path = os.path.join(
                    "/bam_vol",
                    sample.run_id,
                    os.path.basename(sample.bam.path)
                )

                f.write(f"{sample_path}\n")
    
    def run_read_in_bams(self):
        output_path = os.path.join(self.decon_folder, f"{self.run_id}.RData")
        if os.path.exists(output_path):
            logger.info(
                f"Output of Decon ReadInBams for runID: {self.run_id} already exists: {output_path}"
            )
        fasta_dir = os.path.dirname(self.reference_fasta)
        fasta_filename = os.path.basename(self.reference_fasta)

        cmd = [
            self.docker_path, "run",
            "-v", f"{self.decon_folder}:/decon_folder",
            "-v", f"{self.bams_dir}:/bam_dir",
            "-v", f"{self.Bed.dir}:/{self.Bed.volume}",
            "-v", f"{fasta_dir}:/fasta_dir",
            f"{self.decon_image}:{self.decon_version}",
            "ReadInBams.R", 
            "--bams", f"/decon_folder/{self.input_filename}",
            "--bed", f"{self.Bed.volume}/{self.Bed.filename}",
            "--fasta", f"/fasta_dir/{fasta_filename}",
            "--out", f"/decon_folder/{self.run_id}"
        ]

        self.run_cmd(cmd, "DECON ReadInBams")
        
    def run_identify_failed_rois(self):

        output_path = os.path.join(self.decon_folder, f"{self.run_id}_failed_rois.txt")
        if os.path.exists(output_path):
            logger.info(
                f"Output of Decon IdentifyFailures for runID: {self.run_id} already exists: {output_path}"
            )
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
            

