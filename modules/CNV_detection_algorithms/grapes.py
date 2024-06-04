import os
import json

from modules.log import logger
from modules.CNV_detection_algorithms.CNV_algorithm import CNV_Algorithm


class Grapes(CNV_Algorithm):
    def __init__(self, docker_conf, reference_conf, Bed, force_run=False):
        super().__init__(docker_conf, reference_conf, Bed, force_run)
        self.grapes_image = docker_conf.grapes2["image"]
        self.grapes_version = docker_conf.grapes2["version"]
        self.grapes_dir = os.path.join(self.main_dir, "Grapes")
        self.grapes_results_dir = os.path.join(self.results_dir, "Grapes")
        if not os.path.exists(self.grapes_results_dir):
            os.mkdir(self.grapes_results_dir)
        
        
    def get_grapes_bed(self):
        bed_path = self.Bed.path

        grapes_bed_path = bed_path.replace(".bed", ".grapes.bed")
        if os.path.exists(grapes_bed_path):
            self.Bed.set_grapes_bed(grapes_bed_path)
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
        
        self.Bed.set_grapes_bed(grapes_bed_path)
    

    def create_baseline(self, control_samples):
        for sample in control_samples:
            cohort_dir = sample.bam.cohort_bam_dir
            break
        cmd = [
            self.docker_path, "run",
            "-v", f"{self.fasta_dir}:/fasta_dir",
            "-v", f"{self.Bed.dir}:{self.Bed.volume}",
            "-v", f"{cohort_dir}:/cohort_dir",
            "-v", f"{self.grapes_results_dir}:/grapes_results",
            f"{self.grapes_image}:{self.grapes_version}", "grapes2.py",
            "--bam_dir", "/cohort_dir",
            "--output_dir", "/grapes_results",
            "--bed", f"{self.Bed.volume}/{self.Bed.grapes_bed_filename}",
            "--fasta", f"/fasta_dir/{self.fasta_filename}",
            "--genome", "hg19",
            "--baseline_db", "/grapes_results/grapes2_baseline.db"
        ]

        # runs_volumes = self.get_runs_volumes(control_samples, analysis_samples)
        # docker_control_samples = self.get_docker_samples(control_samples)
        # docker
        #     f"{self.grapes_image}:{self.grapes_version}", "grapes2.py",
        #     "--cases", "/analysis_symbolic_link/", 
        #     "--control", "/cohort_symbolic_link/",
        #     "--output", "/grapes_results",
        #     "--bed", f"{self.Bed.volume}/{self.Bed.grapes_bed_filename}",
        #     "--fasta", f"/fasta_dir/{self.fasta_filename}",
        #     "--genome", "hg19",
        #     "--all"
        # ]

        self.run_cmd(cmd, "Grapes to create a baseline")


    def get_runs_volumes(self, control_samples, analysis_samples):
        """
        It gets all the volumes of all the runs analysed by the pipeline
        """
        runs_analysed = list()
        run_volumes = list()
        volume = "-v"
        for sample in control_samples:
            run_path = os.path.dirname(sample.bam.path)
            run_id = sample.run_id
            if run_path not in runs_analysed:
                run_volume = f"{run_path}:/{run_id}"
                runs_analysed.append(run_path)
                run_volumes.append(volume)
                run_volumes.append(run_volume)
        
        for sample in analysis_samples:
            run_path = os.path.dirname(sample.bam.path)
            run_id = sample.run_id
            if run_path not in runs_analysed:
                run_volume = f"{run_path}:/{run_id}"
                runs_analysed.append(run_path)
                run_volumes.append(volume)
                run_volumes.append(run_volume)
        
        return run_volumes
    
        
    def get_docker_samples(self, samples):
        """
        Get the relative docker path of the samples based on the volumes created by runs.
        
        E.g. RUN20230405 will create a volume in docker named RUN20230405,
        this function will obtain the relative path of each sample in the docker volume:
        /RUN20230405/RB35342.bam
        """
        sample_docker_path = list()
        for sample in samples:
            run_id = sample.run_id
            print(f"runid: {run_id}")
            sample_docker_path.append(
                f"/{run_id}/{sample.bam.filename}")
        
        return(sample_docker_path)
    

