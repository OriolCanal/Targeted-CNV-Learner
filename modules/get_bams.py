import os
import subprocess
from ftplib import FTP
from modules.mongo_classes import Run
from modules.log import logger


def get_samples_from_run_name(Run, run_id):
    try:
        run_doc = Run.objects.get(run_id=run_id)
    except:
        raise (ValueError(
            f"Run ID: {run_id} not found in MongoDb!"
        )
        )
    
    if hasattr(run_doc, "samples"):
        return(run_doc.samples)
    else:
        raise (ValueError(
            f"Run ID: {run_id} does not contain any sample"
        ))

def get_bam_bai_from_compendi(sample_id, run_id, ref_conf):
    compendi_bai_path = "http://172.16.83.24:8001/download_sample_bai/"
    compendi_bam_path = "http://172.16.83.24:8001/download_sample_bam/"
    output_dir = os.path.join(ref_conf.main_dir, "runs")

    if not os.path.isdir(output_dir):
        logger.info(f"Creating directory {output_dir} to store bam bai files.")
        os.mkdir(output_dir)

    run_dir = os.path.join(output_dir, run_id)
    if not os.path.isdir(run_dir):
        os.mkdir(run_dir)
    sample_bam = f"{sample_id}.bam"
    sample_bai = f"{sample_id}.bai"
    sample_bam_path = os.path.join(run_dir, sample_bam)
    sample_bai_path = os.path.join(run_dir, sample_bai)
    if not os.path.isfile(sample_bam_path):
        download_bam_compendi = os.path.join(compendi_bam_path, sample_id, run_id)

        result = subprocess.run(["wget", "-O", sample_bam_path, download_bam_compendi])

        if result.returncode == 0:
            logger.info(f"{sample_bam_path} downloaded successfully!")
        else:
            raise RuntimeError(
                f"Download from {download_bam_compendi} failes with return code: {result.returncode}"
            )
    else:
        logger.info(
            f"Bam file: {sample_bam_path} has previously been downloaded!"
        )
    
    if not os.path.isfile(sample_bai_path):
        download_bai_compendi = os.path.join(compendi_bai_path, sample_id, run_id)

        result_bai = subprocess.run(["wget", "-O", sample_bai_path, download_bai_compendi])

        if result_bai.returncode == 0:
            logger.info(f"{sample_bam_path} downloaded successfully!")
        else:
            raise RuntimeError(
                f"Download from {download_bam_compendi} failes with return code: {result.returncode}"
            )
    else:
        logger.info(
            f"Bai file: {sample_bai_path} has previously been downloaded!"
        )

    return(sample_bam_path, sample_bai_path)

def get_bams_from_run_name(Run, run_id, ref_conf):

    samples = get_samples_from_run_name(Run, run_id)
    downloaded_run_samples = set()
    for sample in samples:
        sample_bam, sample_bai = get_bam_bai_from_compendi(sample, run_id, ref_conf)
        downloaded_run_samples.add
        



# get_bams_from_run_name(Run, "RUN20240415-CGC6400")