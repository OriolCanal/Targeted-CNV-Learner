import os
import subprocess
import shutil

from modules.log import logger
from modules.mongo_classes import Sample
from modules.bam import Bam
class Analysis_Run():
    def __init__(self, run_id):
        self.run_id = run_id
        self.samples_ids = list()
        self.samples_147 = list()
    
    def get_samples(self, Run_mongo, Sample_mongo):
        try:
            run_doc = Run_mongo.objects.get(run_id=self.run_id)
        except:
            logger.warning(
                f"Run ID: {self.run_id} not found in MongoDb!"
            )
        
        if hasattr(run_doc, "samples"):
            samples = run_doc.samples
            for sample_id in samples:
                Sample = Analysis_Sample(run_id=self.run_id, sample_id=sample_id)
                self.samples_ids.append(Sample)
                if Sample.is_panel_147(Sample_mongo):
                    self.samples_147.append(Sample)

        else:
            raise (ValueError(
                f"Run ID: {self.run_id} does not contain any sample"
            ))
        

class Analysis_Sample():
    compendi_bai_path = "http://172.16.83.24:8001/download_sample_bai/"
    compendi_bam_path = "http://172.16.83.24:8001/download_sample_bam/"

    def __init__(self, run_id, sample_id):
        self.run_id = run_id
        self.sample_id = sample_id
        self.downloaded_bam = False
        self.bam = None # Bam class containing info about bam files
        self.cluster = None
    def is_panel_147(self, Sample_mongo):
        try:
            sample_doc = Sample_mongo.objects.get(run_id=self.run_id, lab_id=self.sample_id)
        except:
            logger.warning(
                f"Sample ID: {self.sample_id} not found in MongoDb Sample collection!"
            )
            return False
        
        if hasattr(sample_doc, "panel"):
            if "SUDD_147" in sample_doc.panel:
                self.panel_147 = True
                return(True)
            elif sample_doc.panel == "147":
                self.panel_147 = True
                return(True)
            self.panel_147 = False
            return False
        
        raise ValueError(
            f"Sample document with ID {self.sample_id} does not contain a panel attribute!"
        )
    
            
    def get_bam_bai_from_compendi(self, ref_conf):
        if not self.panel_147:
            logger.warning(f"Bam file won't be downloaded as it's not panel SUDD_147")
            return(0)
        output_dir = os.path.join(ref_conf.main_dir, "runs")

        if not os.path.isdir(output_dir):
            logger.info(f"Creating directory {output_dir} to store bam bai files.")
            os.mkdir(output_dir)

        run_dir = os.path.join(output_dir, self.run_id)
        if not os.path.isdir(run_dir):
            os.mkdir(run_dir)
        sample_bam = f"{self.sample_id}.bam"
        sample_bai = f"{self.sample_id}.bai"
        bam_path = os.path.join(run_dir, sample_bam)
        bai_path = os.path.join(run_dir, sample_bai)
        if not os.path.isfile(bam_path):
            download_bam_compendi = os.path.join(Analysis_Sample.compendi_bam_path, self.sample_id, self.run_id)

            result = subprocess.run(["wget", "-O", bam_path, download_bam_compendi])

            if result.returncode == 0:
                logger.info(f"{bam_path} downloaded successfully!")
                self.bam = Bam(bam_path)
            else:
                logger.error(
                    f"Download from {download_bam_compendi} failed with return code: {result.returncode}"
                )
                if os.path.isfile(bam_path):
                    shutil.rmtree(run_dir)
                return(False)
        else:
            self.bam = Bam(bam_path)
            logger.info(
                f"Bam file: {bam_path} has previously been downloaded!"
            )

        if not os.path.isfile(bai_path):
            download_bai_compendi = os.path.join(Analysis_Sample.compendi_bai_path, self.sample_id, self.run_id)

            result_bai = subprocess.run(["wget", "-O", bai_path, download_bai_compendi])

            if result_bai.returncode == 0:
                logger.info(f"{bai_path} downloaded successfully!")
                self.bam = Bam(bam_path)
            else:
                logger.error(
                    f"Download from {download_bai_compendi} failed with return code: {result.returncode}"
                )
                if os.path.isfile(bai_path):
                    shutil.rmtree(run_dir)
                return(False)
        else:
            logger.info(
                f"Bai file: {bai_path} has previously been downloaded!"
            )

        return(self.bam)
        
