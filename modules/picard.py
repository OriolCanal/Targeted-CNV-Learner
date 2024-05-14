import subprocess
import os
import pandas as pd
from modules.log import logger

class Picard():
    def __init__(self, docker_conf, ref_conf):
        self.docker_path = "/usr/bin/docker"
        self.docker_conf = docker_conf    
        self.reference_dir = ref_conf.hg19.dir_path
        self.ref_dict_file = ref_conf.hg19.dict
        self.ref_fasta_file = ref_conf.hg19.fasta

    def create_fasta_dict(self):
        self.ref_dict_file = self.ref_fasta_file.replace("fasta", "dict")
        picard_cmd = [
            "docker", "run",
            "-v", f"{self.reference_dir}:/ref_dir",
            self.docker_conf.picard["image"],
            "java", "-Xmx60g", "-jar",
            "/usr/picard/picard.jar",
            "CreateSequenceDictionary",
            "-R", f"/ref_dir/{self.ref_fasta_file}",
            "-O", f"/ref_dir/{self.ref_dict_file}"
        ]

        cmd_str = " ".join(picard_cmd)
        logger.info(
            f"Creating dict file from fasta, if already exists it won't be created:\n{cmd_str}"
        )
        subprocess.run(picard_cmd, encoding="utf-8", capture_output=True)
        return(self.ref_dict_file)

    def run_bed_to_interval_list(self, Bed):
        
        interval_list_filename = f"{Bed.filename}.interval_list"
        interval_list_path = os.path.join(Bed.dir, interval_list_filename)
        if os.path.isfile(Bed.interval_list_path):
            logger.info(
                f"Listfile for bed {interval_list_path} already created in {interval_list_path}"
            )
            Bed.set_interval_list_path(interval_list_path)
            return(interval_list_path)
        picard_cmd = [
            "docker", "run",
            "-v", f"{Bed.dir}:/bed_dir",
            "-v", f"{self.reference_dir}:/ref_dir",
            self.docker_conf.picard["image"],
            "java", "-Xmx60g", "-jar",
            "/usr/picard/picard.jar",
            "BedToIntervalList",
            "-I", f"/bed_dir/{Bed.filename}",
            "-O", f"/bed_dir/{interval_list_filename}",
            "-SD", f"/ref_dir/{self.ref_dict_file}"
        ]

        command_str = " ".join(picard_cmd)
        logger.info(
            f"Creating picard interval list:\n{command_str}"
        )
        subprocess.run(picard_cmd, encoding="utf-8", capture_output=True)
        Bed.set_interval_list_path(interval_list_path)

        return(interval_list_path)
    
    def run_collectHsMetrics(self, Bam, Bed):
        if not hasattr(Bed, "interval_list_path"):
            self.run_bed_to_interval_list(Bed)

        hsmetrics_filename = f"{Bam.filename}_hs_metrics.txt"
        picard_dirname = "Picard"
        picard_dir = os.path.join(Bam.dir, picard_dirname)
        hsmetrics_path = os.path.join(picard_dir, hsmetrics_filename)
        if os.path.isfile(hsmetrics_path):
            logger.info(
                f"HsMetrics output file already available: {hsmetrics_path}"
            )
            Bam.set_hs_metrics(hsmetrics_path)
            return(hsmetrics_path)

        if not os.path.isdir(picard_dir):
            os.mkdir(picard_dir)
        picard_cmd = [
            "docker", "run",
            "-v", f"{Bed.dir}:/bed_dir",
            "-v", f"{self.reference_dir}:/ref_dir",
            "-v", f"{Bam.dir}:/bam_dir",
            self.docker_conf.picard["image"],
            "java", "-Xmx60g", "-jar",
            "/usr/picard/picard.jar",
            "CollectHsMetrics",
            "--INPUT", f"/bam_dir/{Bam.filename}",
            "--OUTPUT", f"/bam_dir/{picard_dirname}/{hsmetrics_filename}",
            "--BAIT_INTERVALS", f"/bed_dir/{Bed.interval_list_filename}",
            "--TARGET_INTERVALS", f"/bed_dir/{Bed.interval_list_filename}"
        ]

        command_str = " ".join(picard_cmd)
        logger.info(
            f"Running picard CollectHsMetrics:\n{command_str}"
        )
        subprocess.run(picard_cmd, encoding="utf-8", capture_output=True)
        Bam.set_hs_metrics(hsmetrics_path)

        return(hsmetrics_path)

    def get_picard_metrics(self, Bam):

        with open(Bam.hs_metrics_path) as f:
            for line in f:
                if not line.startswith("#"):
                    continue
                if line.startswith("## METRICS CLASS"):
                    metrics_info_line = next(f)
                    
                    metrics_value_line = next(f)
                
                    return(metrics_info_line, metrics_value_line)


class Metrics_Df():
    def __init__(self):
        self.metrics_df = pd.DataFrame()
        self.df_has_header = False

    def has_df_header(self):
        return not self.metrics_df.columns.empty

    def add_metrics_header(self, line):
        header = line.strip().split("\t")
        # this items are not specified and should be removed
        items_to_remove = ["LIBRARY", "SAMPLE", "READ_GROUP"]
        for item in items_to_remove:
            if item in header:
                header.remove(item)
        
            
        self.metrics_df = pd.DataFrame(columns=header)  # Initialize DataFrame with header
        self.df_has_header = True

    def add_metrics_line(self, values_line):
        metrics_values = values_line.strip().split("\t")

        if len(metrics_values) != len(self.metrics_df.columns):
            raise ValueError("Length of metrics values does not match the length of columns.")
        
        self.metrics_df = self.metrics_df._append(dict(zip(self.metrics_df.columns, metrics_values)), ignore_index=True)
    