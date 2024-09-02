import os
import subprocess

from modules.log import logger

from modules.CNV_detection_algorithms.CNV_algorithm import CNV_Algorithm
from modules.detected_CNV import Detected_CNV
class CNV_Kit(CNV_Algorithm):
    def __init__(self, docker_conf, reference_conf, Bed, force_run=False):
        super().__init__(docker_conf, reference_conf, Bed, force_run)
        self.cnv_kit_image = docker_conf.cnvkit["image"]
        self.cnv_kit_version = docker_conf.cnvkit["version"]
        self.cnv_kit_dir = os.path.join(self.main_dir, "CNVKit")
        self.cnv_kit_results_dir = os.path.join(self.results_dir, "CNVKit")
        self.normal_cnn_filename = "reference.cnn" 
        self.normal_cnn = os.path.join(self.cnv_kit_dir, self.normal_cnn_filename)
        self.results_folder = "CNVKit_results"
        self.resutls_dir = os.path.join(self.cnv_kit_results_dir, self.results_folder)

        if not os.path.exists(self.cnv_kit_dir):
            os.mkdir(self.cnv_kit_dir)

        if not os.path.exists(self.cnv_kit_results_dir):
            os.mkdir(self.cnv_kit_results_dir)
        
    def run_batch_germline_pipeline(self, control_samples, analysis_samples):
        """
        Run batch CNV-Kit pipeline
        
        Params:
            control_samples(list): list of control samples objects that will be used to create the reference cnn
            analysis_samples(list): list of analysis samples objects that will be analyse for CNVs.
        """
        fasta_dir = os.path.dirname(self.reference_fasta)
        fasta_filename = os.path.basename(self.reference_fasta)

        # Getting samples already analysed
        samples_to_analyse = list()
        for sample in analysis_samples:
            output_filename = sample.bam.filename.replace(".bam", ".call.cns")
            output_path = os.path.join(self.resutls_dir, output_filename)
            # print(output_path)
            if os.path.exists(output_path):
                logger.info(
                    f"CNVKit won't be run over sample {sample.sample_id} as output already exists: {output_path}"
                )
                continue
            print(output_path)
            samples_to_analyse.append(sample)


        fasta_volume = "/fasta_dir"
        cmd = [
            self.docker_path, "run",
            "-v", f"{self.Bed.dir}:{self.Bed.volume}",
            "-v", f"{fasta_dir}:{fasta_volume}",
            "-v", f"{self.cnv_kit_dir}:/cnv_kit_dir",
            "-v", f"{self.cnv_kit_results_dir}:/cnv_kit_results",

        ]
        runs_volumes = self.get_runs_volumes(control_samples, samples_to_analyse)
        
        docker_control_samples = self.get_docker_samples(control_samples)
        docker_control_samples.insert(0, "--normal")
        docker_analysis_samples = self.get_docker_samples(samples_to_analyse)
        if not docker_analysis_samples or len(docker_analysis_samples) < 2:
            return True
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
                f"--targets", f"{self.Bed.volume}/{self.Bed.roi_bed_filename}",
                "--fasta", f"{fasta_volume}/{fasta_filename}",
                "--output-reference", f"/cnv_kit_dir/{self.normal_cnn_filename}",
                "--output-dir", f"/cnv_kit_results/{self.results_folder}"
            ]
        )

        self.run_cmd(cmd, "CNVKit copy number calling pipeline")
    
    def run_batch_germline_pipeline_with_existing_reference(self, analysis_samples, cpus=4):

        cmd = [
            self.docker_path, "run",

            "-v", f"{self.cnv_kit_dir}:/cnv_kit_dir",
            "-v", f"{self.cnv_kit_results_dir}:/cnv_kit_results",
            

        ]
        # no control samples in this case
        runs_volumes = self.get_runs_volumes([], analysis_samples)

        cmd.extend(runs_volumes)
        cmd.extend(
            [
                f"{self.cnv_kit_image}:{self.cnv_kit_version}",
                "cnvkit.py", "batch",
            ]
        )
        docker_analysis_samples = self.get_docker_samples(analysis_samples)

        cmd.extend(docker_analysis_samples)

        cmd.extend(
            [
                "-r", f"/cnv_kit_dir/{self.normal_cnn_filename}"
                "-p", cpus
            ]
        )


    def get_vcf(self, analysis_samples):
        for sample in analysis_samples:
            call_cns_filename = sample.bam.filename.replace(".bam", ".call.cns")
            output_path = os.path.join(self.results_folder, call_cns_filename)
            vcf_filename = f"{sample.sample_id}_cnvkit.vcf"
            vcf_path = os.path.join(self.results_folder, vcf_filename)
            if not os.path.exists(output_path):
                raise FileNotFoundError(
                    f"CNVKit cns file does not exists: {output_path}"
                )
            
            output_dir = os.path.join(self.cnv_kit_results_dir, self.results_folder)
            cmd = [
                self.docker_path, "run",
                "-v", f"{output_dir}:/cnv_kit_results",
                f"{self.cnv_kit_image}:{self.cnv_kit_version}",
                "cnvkit.py",
                "export", "vcf",
                f"/cnv_kit_results/{call_cns_filename}",
                "-o", f"/cnv_kit_results/{vcf_filename}"

            ]

            self.run_cmd(cmd, "CNVKit cns export to vcf:")
            sample.set_cnv_kit_vcf(vcf_path)
    def get_bed(self, analysis_samples):

        for sample in analysis_samples:
            call_cns_filename = sample.bam.filename.replace(".bam", ".call.cns")
            output_path = os.path.join(self.results_folder, call_cns_filename)
            bed_filename = f"{sample.sample_id}_cnvkit.bed"
            if not os.path.exists(output_path):
                raise FileNotFoundError(
                    f"CNVKit cns file does not exists: {output_path}"
                )
            
            output_dir = os.path.join(self.cnv_kit_results_dir, self.results_folder)
            cmd = [
                self.docker_path, "run",
                "-v", f"{output_dir}:/cnv_kit_results",
                f"{self.cnv_kit_image}:{self.cnv_kit_version}",
                "cnvkit.py",
                "export", "vcf",
                f"/cnv_kit_results/{call_cns_filename}",
                "-o", f"/cnv_kit_results/{bed_filename}"

            ]

            self.run_cmd(cmd, "CNVKit cns export to bed:")
         
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
            sample_docker_path.append(
                f"/{run_id}/{sample.bam.filename}")
        
        return(sample_docker_path)
    

    def parse_detected_cnvs(self, Bed_obj, sample_id_sample_obj: dict =None):
        out_files = os.listdir(self.resutls_dir)
        calls_files = [os.path.join(self.resutls_dir, out_file) for out_file in out_files if out_file.endswith(".call.cns")]

        for calls_file in calls_files:
            sample = os.path.basename(calls_file).split(".")[0]
            # print(calls_file)
            with open(calls_file, "r") as f:
                for line in f:
                    if line.startswith("chromosome"):
                        continue
                    line = line.strip().split("\t")

                    cn = line[5]
                    if int(cn) == 2:
                        continue
                    chromosome = line[0]
                    if chromosome == "chrX" or chromosome == "chrY":
                        continue

                    if int(cn) > 2:
                        sv_type = "DUP"
                    elif int(cn) < 2: 
                        sv_type = "DEL"
                    
                    start, end, gene, = line[1], line[2], line[3]
                    # print("cnv found")
                    cnv =Detected_CNV(start, end, chromosome, sv_type, sample, None, gene, "CNVKit")
                    cnv.get_numb_exons(Bed_obj)
                    if sample_id_sample_obj:
                        sample_obj = sample_id_sample_obj[sample]
                        sample_obj.cnvs["cnvkit"].append(cnv)