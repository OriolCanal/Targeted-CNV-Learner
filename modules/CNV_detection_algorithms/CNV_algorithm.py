import os
import subprocess
from modules.log import logger

class CNV_Algorithm():
    def __init__(self, docker_conf, reference_conf, Bed, force_run=False):
        self.docker_path = "/usr/bin/docker"
        self.main_dir = reference_conf.main_dir
        self.Bed = Bed
        self.reference_fasta = reference_conf.hg19.fasta_path
        self.dict_filename = reference_conf.hg19.dict
        self.results_dir = os.path.join(self.main_dir, "Results")
        if not os.path.exists(self.results_dir):
            os.mkdir(self.results_dir)
        self.force_run = force_run
    
    def run_cmd(self, cmd, cmd_description=""):

        cmd_str = " ".join(cmd)
        logger.info(
            f"Running {cmd_description}:\n{cmd_str}"
        )
        result = subprocess.run(cmd, encoding="utf-8", capture_output=True)

        if result.returncode == 0:
            logger.info("Command ran successfully!")
        else:
            raise(ValueError(
                f"Command Failed:\n{cmd_str}"
            ))