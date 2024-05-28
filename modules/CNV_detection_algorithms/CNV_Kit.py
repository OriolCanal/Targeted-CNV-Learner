import os
import subprocess

from modules.log import logger

from modules.CNV_detection_algorithms.CNV_algorithm import CNV_Algorithm

class CNV_Kit(CNV_Algorithm):
    def __init__(self, docker_conf, reference_conf, Bed, force_run=False):
        super().__init__(docker_conf, reference_conf, Bed, force_run)
        self.cnv_kit_image = docker_conf.cnvkit["image"]
        self.cnv_kit_version = docker_conf.cnvkit["version"]
        self.decon_dir = os.path.join(self.main_dir, "DECON")
        self.decon_results_dir = os.path.join(self.results_dir, "Decon")
        self.normal_cnn_filename = "reference.cnn" 
        self.normal_cnn = os.path.join(self.decon_dir, self.normal_cnn_filename)
        self.results_filename = "decon_results"
        self.resutls_dir = os.path.join(self.decon_results_dir, self.results_filename)


        if not os.path.exists(self.decon_results_dir):
            os.mkdir(self.decon_results_dir)
        
    def run_batch_germline_pipeline(self, control_samples, analysis_samples):
        """
        Run batch CNV-Kit pipeline
        
        Params:
            control_samples(list): list of control samples objects that will be used to create the reference cnn
            analysis_samples(list): list of analysis samples objects that will be analyse for CNVs.
        """
        fasta_dir = os.path.dirname(self.reference_fasta)
        fasta_filename = os.path.basename(self.reference_fasta)

        fasta_volume = "/fasta_dir"
        cmd = [
            self.docker_path, "run",
            "-v", f"{self.Bed.dir}:{self.Bed.volume}",
            "-v", f"{fasta_dir}:{fasta_volume}",
            "-v", f"{self.decon_dir}:/decon_dir",
            "-v", f"{self.decon_results_dir}:/decon_results",

        ]
        runs_volumes = self.get_runs_volumes(control_samples, analysis_samples)
        
        docker_control_samples = self.get_docker_samples(control_samples)
        docker_control_samples.insert(0, "--normal")
        docker_analysis_samples = self.get_docker_samples(analysis_samples)
        cmd.extend(runs_volumes)
        cmd.extend(
            [
                f"{self.cnv_kit_image}:{self.cnv_kit_version}",
                "cnvkit.py", "batch",
            ]
        )
        cmd.extend(docker_analysis_samples)
        cmd.extend(docker_control_samples)
        cmd.extend(
            [
                f"--targets", f"{self.Bed.volume}/{self.Bed.filename}",
                "--fasta", f"{fasta_volume}/{fasta_filename}",
                "--output-reference", f"/decon_dir/{self.normal_cnn_filename}",
                "--output-dir", f"/decon_results/{self.results_filename}"
            ]
        )

        self.run_cmd(cmd, "CNVKit copy number calling pipeline")
    

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
                run_volume = f"{run_path}:/{run_id.lower()}"
                runs_analysed.append(run_path)
                run_volumes.append(volume)
                run_volumes.append(run_volume)
        
        for sample in analysis_samples:
            run_path = os.path.dirname(sample.bam.path)
            run_id = sample.run_id
            if run_path not in runs_analysed:
                run_volume = f"{run_path}:/{run_id.lower()}"
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
                f"/{run_id.lower()}/{sample.bam.filename}")
        
        return(sample_docker_path)
    

