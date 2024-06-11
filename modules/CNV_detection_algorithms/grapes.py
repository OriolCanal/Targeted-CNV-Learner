import os
import json

from modules.log import logger
from modules.CNV_detection_algorithms.CNV_algorithm import CNV_Algorithm
from modules.detected_CNV import Detected_CNV

class Grapes(CNV_Algorithm):
    def __init__(self, docker_conf, reference_conf, Run_class, Bed, force_run=False):
        super().__init__(docker_conf, reference_conf, Bed, force_run)
        self.grapes_image = docker_conf.grapes2["image"]
        self.grapes_version = docker_conf.grapes2["version"]
        self.Run = Run_class
        self.run_id = Run_class.run_id
        self.bams_dir = self.Run.run_path
        self.grapes_dir = os.path.join(self.main_dir, "Grapes")
        self.grapes_results_dir = os.path.join(self.results_dir, "Grapes")
        self.run_results_dir = os.path.join(self.grapes_results_dir, self.run_id)
        if not os.path.exists(self.grapes_results_dir):
            os.mkdir(self.grapes_results_dir)
            os.chmod(self.run_results_dir, 0o777)
        if not os.path.exists(self.run_results_dir):
            os.mkdir(self.run_results_dir)
            os.chmod(self.run_results_dir, 0o777)
        if not os.path.exists(self.grapes_dir):
            os.mkdir(self.grapes_dir)
            os.chmod(self.grapes_dir, 0o777)
        self.input_filename = f"{self.run_id}_grapes_input.txt"
        self.input_path = os.path.join(self.grapes_dir, self.input_filename)
        
        self.output_file = os.path.join(self.run_results_dir, "output_dir.all.calls.bed")

    def create_input_file(self):
        samples_analysed = 0
        with open(self.input_path, "w") as f:
            for sample in self.Run.samples_147:
                if sample.is_outlier is not False:
                    continue
                samples_analysed += 1
                sample_path = os.path.join(
                    "/bam_vol",
                    os.path.basename(sample.bam.path)
                )

                f.write(f"{sample_path}\n")
        
        if samples_analysed > 1:
            return True
        else:
            return False

    def run_grapes_run_mode(self):
        if os.path.exists(self.output_file):
            logger.info(f"GRAPES2 already run for RUN: {self.run_id} as output: {self.output_file} exists")
            return(True)
        fasta_dir = os.path.dirname(self.reference_fasta)
        fasta_filename = os.path.basename(self.reference_fasta)
        cmd = [
            self.docker_path, "run",
            "-v", f"{fasta_dir}:/fasta_dir",
            "-v", f"{self.run_results_dir}:/output_dir",
            "-v", f"{self.grapes_dir}:/input_dir",
            "-v", f"{self.bams_dir}:/bam_vol",
            "-v", f"{self.Bed.dir}:{self.Bed.volume}",
            f"{self.grapes_image}:{self.grapes_version}",
            "grapes2.py", 
            "--bam_dir", f"/input_dir/{self.input_filename}",
            "--output_dir", "/output_dir",
            "--bed", f"{self.Bed.volume}/{self.Bed.grapes_bed_filename}",
            "--fasta", f"/fasta_dir/{fasta_filename}",
            "--genome_version", "hg19"
        ]

        self.run_cmd(cmd, "GRAPES in run mode:")
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
            "--baseline_db", "/grapes_results/grapes2_baseline.db",
            "--use_baseline"
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
    
    def parse_dectected_cnvs(self, Bed_obj, sample_id_sample_obj: dict =None):
        self.output_file = os.path.join(self.run_results_dir, "output_dir.all.calls.bed")

        with open(self.output_file, "r") as f:
            for line in f:
                if line.startswith("sample"):
                    continue
                line = line.strip().split("\t")
                sample = line[0].split(".")[0]
                chr, start, end = line[1], line[2], line[3]
                info = line[4]
                infos = info.split(";")
                for info in infos:
                    if "=" not in info:
                        continue
                    key, value = info.split("=")
                    if key == "SVTYPE":
                        svtype = value
                Detected_Cnv = Detected_CNV(start, end, chr, svtype, sample, None, None, "Grapes2")
                Detected_Cnv.get_gene_name(Bed_obj)
                Detected_Cnv.get_numb_exons(Bed_obj)

                if sample_id_sample_obj:
                    sample_obj = sample_id_sample_obj[sample]
                    sample_obj.cnvs["grapes2"].append(Detected_Cnv)

    

