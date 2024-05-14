import os
import subprocess

from modules.log import logger

class Gatk_gCNV():
    """
    Run all the steps of GATK gCNV in:
    https://gatk.broadinstitute.org/hc/en-us/articles/360035531152--How-to-Call-rare-germline-copy-number-variants
    """
    def __init__(self, docker_conf, reference_conf, Bed):

        self.docker_path = "/usr/bin/docker"
        self.main_dir = reference_conf.main_dir
        self.Bed = Bed
        self.gatk_image = docker_conf.gatk["image"]
        self.gatk_version = docker_conf.gatk["version"]
        self.reference_fasta = reference_conf.hg19.fasta_path
        self.dict_filename = reference_conf.hg19.dict
        self.gatk_folder = os.path.join(reference_conf.main_dir, "GATK")
        self.gatk_read_counts_dirname = "Read_counts"
        self.gatk_read_counts_folder = os.path.join(self.gatk_folder, self.gatk_read_counts_dirname)
        self.mappability_folder = os.path.join(self.gatk_folder, "mappability_track")
        self.mappability_track_path = os.path.join(self.mappability_folder, "k50.umap.bed")
        self.gatk_volume = "/gatk_vol"
        self.read_counts_volume = "/read_counts"
        self.contig_ploidy_filename = "contig_ploidy_prior.tsv"
        self.contig_ploidy = os.path.join(self.gatk_folder, self.contig_ploidy_filename)
        results_dir = os.path.join(self.main_dir, "Results")
        if not os.path.exists(results_dir):
            os.mkdir(results_dir)
        self.gatk_results_dir = os.path.join(results_dir, "GATK")
        if not os.path.exists(self.gatk_results_dir):
            os.mkdir(self.gatk_results_dir)
    def run_preprocess_intervals(self, Bed):
        
        fasta_dir = os.path.dirname(self.reference_fasta)
        fasta_filename = os.path.basename(self.reference_fasta)

        fasta_volume = "/fasta_dir"

        preprocessed_interval_list_filename = f"preprocessed_{Bed.interval_list_filename}"
        preprocessed_interval_list = os.path.join(Bed.dir, preprocessed_interval_list_filename)

        if os.path.exists(preprocessed_interval_list):
            logger.info(f"{preprocessed_interval_list} already exists, Gatk PreprocessIntervals won't be run")
            Bed.set_preprocessed_intervals_list_path(preprocessed_interval_list)
        cmd = [
            self.docker_path, "run",
            "-v", f"{fasta_dir}:{fasta_volume}",
            "-v", f"{Bed.dir}:{Bed.volume}",
            f"{self.gatk_image}:{self.gatk_version}",
            "gatk", "PreprocessIntervals",
            "-R", f"{fasta_volume}/{fasta_filename}",
            "-L", f"{Bed.volume}/{Bed.interval_list_filename}",
            "--bin-length", "0",
            "-imr", "OVERLAPPING_ONLY",
            "--padding", "250",
            "-O", f"{Bed.volume}/{preprocessed_interval_list_filename}"
        ]
        cmd_str = " ".join(cmd)
        logger.info(f"Running GATK Preprocess Intervals:\n{cmd_str}")
        subprocess.run(cmd, encoding="utf-8", capture_output=True)
        
        Bed.set_preprocessed_intervals_list_path(preprocessed_interval_list)
        return (preprocessed_interval_list)
    
    def create_gatk_read_counts_folder(self, Bam):


        if not os.path.isdir(self.gatk_read_counts_folder):
            logger.info(
                f"Creating new folder to store GATK gCNV files: {self.gatk_read_counts_folder}"
            )
            os.mkdir(self.gatk_read_counts_folder)
        return(self.gatk_read_counts_folder)
    
    
    def run_collect_read_counts(self, Bam, Bed):
        """
        Run GATK CollectReadCounts.
        """

        self.create_gatk_read_counts_folder(Bam)

        if not Bed.preprocessed_intervals_path:
            raise ValueError(
                f"You need to preprocess the interval list file first! Run the funciton Gatk_gCNV.run_preprocess_intervals()"
            )


        fasta_dir = os.path.dirname(self.reference_fasta)
        fasta_filename = os.path.basename(self.reference_fasta)


        hdf5_filename = f"{Bam.filename}.hdf5"
        hdf5_path = os.path.join(self.gatk_read_counts_folder, hdf5_filename)
        
        if os.path.exists(hdf5_path):
            logger.info(
                f"GATK CollectReadCounts won't be executed as {hdf5_path} file already exists."
            )
            Bam.set_hdf5_path_and_filename(hdf5_path)
            return(Bam.hdf5_path)

        fasta_volume = "/fasta_dir"

        cmd = [
            self.docker_path, "run",
            "-v", f"{self.gatk_read_counts_folder}:/read_counts",
            "-v", f"{fasta_dir}:{fasta_volume}",
            "-v", f"{Bed.dir}:{Bed.volume}",
            "-v", f"{Bam.dir}:{Bam.volume}",
            f"{self.gatk_image}:{self.gatk_version}",
            "gatk", "CollectReadCounts",
            "-L", f"{Bed.volume}/{Bed.preprocessed_intervals_filename}",
            "-R", f"{fasta_volume}/{fasta_filename}",
            "-imr", "OVERLAPPING_ONLY",
            "-I", f"{Bam.volume}/{Bam.filename}",
            "-O", f"/read_counts/{hdf5_filename}"
        ]

        cmd_str = " ".join(cmd)
        logger.info(f"Running GATK CollectReadCounts:\n{cmd_str}")
        subprocess.run(cmd, encoding="utf-8", capture_output=True)
        Bam.set_hdf5_path_and_filename(hdf5_path)
        return(Bam.hdf5_path)
    
    def run_index_feature_file(self, Bed):
        """
        Indexing ./bed/mappability_track/k36.umap.bed.gz required to run annotate intervals
        """

        mappability_folder = os.path.basename(self.mappability_folder)
        mappability_filename = os.path.basename(self.mappability_track_path)

        index_map_track_path = f"{self.mappability_track_path}.idx"
        if os.path.isfile(index_map_track_path):
            return(index_map_track_path)
        cmd = [
            self.docker_path, "run",
            "-v", f"{Bed.dir}:{Bed.volume}",
            f"{self.gatk_image}:{self.gatk_version}",
            "gatk", "IndexFeatureFile",
            "-I", f"{Bed.volume}/{mappability_folder}/{mappability_filename}"
        ]
        cmd_str = " ".join(cmd)
        logger.info(
            f"Indexing mappability track:\n{cmd_str}"
        )
        subprocess.run(cmd, encoding="utf-8", capture_output=True)
        return(index_map_track_path)


    def run_annotate_intervals(self, Bed):
        """
        Exclude problematic reagions based on GC content and mappability of intervals
        """

        if not Bed.preprocessed_intervals_path:
            raise ValueError(
                f"You need to preprocess the interval list file first! Run the funciton Gatk_gCNV.run_preprocess_intervals()"
            )
        
        fasta_dir = os.path.dirname(self.reference_fasta)
        fasta_filename = os.path.basename(self.reference_fasta)


        mappability_dirname = os.path.basename(self.mappability_folder)
        mappability_filename = os.path.basename(self.mappability_track_path)


        gc_annotated_bed_filename = f"annotated_{Bed.preprocessed_intervals_filename}"
        gc_annotated_bed_path = os.path.join(Bed.dir, gc_annotated_bed_filename)
        if os.path.exists(gc_annotated_bed_path):
            logger.info(
                f"GATK Annotate intevals won't be executed as {gc_annotated_bed_path} file already exists."
            )
            Bed.set_annotated_intervals_path(gc_annotated_bed_path)
            return(gc_annotated_bed_path)
        
        fasta_volume = "/fasta_dir"

        cmd = [
            self.docker_path, "run",
            "-v", f"{fasta_dir}:{fasta_volume}",
            "-v", f"{Bed.dir}:{Bed.volume}",
            "-v", f"{self.gatk_folder}:{self.gatk_volume}",
            f"{self.gatk_image}:{self.gatk_version}",
            "gatk", "AnnotateIntervals",
            "-L", f"{Bed.volume}/{Bed.preprocessed_intervals_filename}",
            "-R", f"{fasta_volume}/{fasta_filename}",
            "-imr", "OVERLAPPING_ONLY",
            "--mappability-track", f"{self.gatk_volume}/{mappability_dirname}/{mappability_filename}",
            "-O", f"{Bed.volume}/{gc_annotated_bed_filename}"
        ]
    
        cmd_str = " ".join(cmd)
        logger.info(f"Running GATK AnnotateIntervals:\n{cmd_str}")

        subprocess.run(cmd, encoding="utf-8", capture_output=True)
        Bed.set_annotated_intervals_path(gc_annotated_bed_path)

        return(gc_annotated_bed_path)
    
    def get_input_read_count_files(self, analysed_sample_id, cohort_samples):
        """
        Get a list of input read counts files found in runs/GATK-gCNV7 to be given to gatk docker.
        E.g. ["-I", "gatk_vol/RB35645.hdf5", "-I", "gatk_vol/RB34532.hdf5]"""
        sample_hdf5_counts_filenames = [f"{sample.sample_id}.bam.hdf5" for sample in cohort_samples if sample.sample_id != analysed_sample_id]
        read_counts_input = list()
        for sample_hdf5_filename in sample_hdf5_counts_filenames:

            input_files_cmd = ["-I", f"{self.read_counts_volume}/{sample_hdf5_filename}"]
            read_counts_input.extend(input_files_cmd)
        
        return(read_counts_input)
    
    def run_filter_intervals(self, Bed, Sample, cohort_samples):
        """
        Given specific intervals (annotated intervals), and counts output by CollectReadCoutns,
        outputs a filtered Picard interval list
        """

        filtered_intervals_filename = f"filtered_{Bed.preprocessed_intervals_filename}"
        filtered_intervals_path = os.path.join(Bed.dir, filtered_intervals_filename)

        if os.path.isfile(filtered_intervals_path):
            logger.info(
                f"Filtered intervals file already exists: {filtered_intervals_path}, GATK FilteredIntervals won't be run"
            )
            Bed.set_filtered_intervals_path(filtered_intervals_path)
            return(filtered_intervals_path)

        cmd = [
            self.docker_path, "run",
            "-v", f"{Bed.dir}:{Bed.volume}",
            "-v", f"{self.gatk_read_counts_folder}:{self.read_counts_volume}",
            f"{self.gatk_image}:{self.gatk_version}",
            "gatk", "FilterIntervals",
            "-L", f"{Bed.volume}/{Bed.preprocessed_intervals_filename}",
            "--annotated-intervals", f"{Bed.volume}/{Bed.c}",
        ]

        read_counts_cmd = self.get_input_read_count_files(Sample.id, cohort_samples)
        
        cmd.extend(read_counts_cmd)

        cmd2 = [
            "-imr", "OVERLAPPING_ONLY",
            "-O", f"{Bed.volume}/{filtered_intervals_filename}",
        ]

        cmd.extend(cmd2)

        cmd_str = " ".join(cmd)
        logger.info(
            f"Running GATK FilterIntervals:\n {cmd_str}"
        )
        subprocess.run(cmd, encoding="utf-8", capture_output=True)

        Bed.set_filtered_intervals_path(filtered_intervals_path)

        return(filtered_intervals_path)
  
    def run_interval_list_tools(self, Bed):
        scatter_dirname = "scatter"
        cmd = [
            self.docker_path, "run",
            "-v", f"{Bed.dir}:{Bed.volume}",
            "-v", f"{self.model_dir}:/model_dir",
            f"{self.gatk_image}:{self.gatk_version}",
            "gatk", "IntervalListTools",
            "--INPUT", f"{Bed.volume}/{Bed.filtered_intervals_filename}",
            "--SUBDIVISION_MODE", "INTERVAL_COUNT",
            "--SCATTER_CONTENT", "5000",
            "--OUTPUT", F"{Bed.volume}/{scatter_dirname}"
        ]

        cmd_str = " ".join(cmd)
        logger.info(
            f"Running GATK IntervalListTools:\n{cmd_str}"
        )
        subprocess.run(cmd, encoding="utf-8", capture_output=True)
        scatter_path = os.path.join(Bed.dir, scatter_dirname)
        Bed.set_scatter_path(scatter_path)



class Case_Gatk_gCNV(Gatk_gCNV):
    def __init__(self, sample, docker_conf, reference_conf, Bed):
        self.sample = sample
        # Initialize attributes that belong to the parent class
        super().__init__(docker_conf, reference_conf, Bed)
    

    def create_model(self, Bed, Sample, cohort_samples):

        self.model_dir = os.path.join(self.gatk_folder, self.sample.sample_id)
        if not os.path.exists(self.model_dir):
            os.mkdir(self.model_dir)
        

        self.ploidy_prefix = "ploidy"
        output = os.path.join(self.model_dir, f"{self.ploidy_prefix}-model")
        if os.path.exists(output):
            logger.info(
                f"DetermineGermlineContigPloidy already run in cohort mode for sample {self.sample.sample_id}"
            )
            return(True)
        
        
        input_read_files = self.get_input_read_count_files(Sample.sample_id, cohort_samples)
        cmd = [
            self.docker_path, "run",
            "-v", f"{self.model_dir}:/model_dir", 
            "-v", f"{Bed.dir}:{Bed.volume}",
            "-v", f"{self.gatk_folder}:{self.gatk_volume}",
            "-v", f"{self.gatk_read_counts_folder}:{self.read_counts_volume}",
            f"{self.gatk_image}:{self.gatk_version}",
            "gatk", "DetermineGermlineContigPloidy",
            "-L", f"{Bed.volume}/{Bed.filtered_intervals_filename}",
            "--interval-merging-rule", "OVERLAPPING_ONLY",
        ]

        cmd.extend(input_read_files)

        cmd2 = [
            "--contig-ploidy-priors", f"{self.gatk_volume}/{self.contig_ploidy_filename}",
            "--output", f"/model_dir",
            "--output-prefix", self.ploidy_prefix,
            "--verbosity", "DEBUG"
        ]


        cmd.extend(cmd2)

        cmd_str = " ".join(cmd)

        logger.info(
            f"Running GATK DetermineGermlineContigPloidy:\n{cmd_str}"
        )

        subprocess.run(cmd, encoding="utf-8", capture_output=True)

    def run_determine_germline_contig_ploidy(self):
        self.ploidy_case_prefix = "ploidy-case"
        if not hasattr(self, "model_dir"):
            self.model_dir = os.path.join(self.gatk_folder, self.sample.sample_id)
        if not os.path.exists(self.model_dir):
            os.mkdir(self.model_dir)

        output = os.path.join(self.model_dir, f"{self.ploidy_case_prefix}-model")
        if os.path.exists(output):
            logger.info(
                f"DetermineGermlineContigPloidy already run in case mode for sample {self.sample.sample_id}"
            )
            return(True)
        cmd = [
            self.docker_path, "run",
            "-v", f"{self.gatk_folder}:{self.gatk_volume}",
            "-v", f"{self.model_dir}:/model_dir",
            "-v", f"{self.gatk_read_counts_folder}:{self.read_counts_volume}",
            f"{self.gatk_image}:{self.gatk_version}",
            "gatk", "DetermineGermlineContigPloidy",
            "--model", f"/model_dir/{self.ploidy_prefix}-model",
            "-I", f"{self.read_counts_volume}/{self.sample.sample_id}.bam.hdf5",
            "-O", "/model_dir",
            "--output-prefix", self.ploidy_case_prefix,
            "--verbosity", "DEBUG"
        ]

        cmd_str = " ".join(cmd)
        logger.info(
            f"Running GATK DetermineGermlineContigPloidy over sample {self.sample.sample_id}:\n{cmd_str}"
        )
        subprocess.run(cmd, encoding="utf-8", capture_output=True)
  

    def run_germline_cnv_caller_cohort_mode(self, Bed, Sample, cohort_samples):
        self.cohort_caller_prefix = f"cohort"
        if not hasattr(self, "model_dir"):
            self.model_dir = os.path.join(self.gatk_folder, self.sample.sample_id)
        cmd = [
            self.docker_path, "run",
            "-v", f"{Bed.dir}:{Bed.volume}",
            "-v", f"{self.model_dir}:/model_dir",
            "-v", f"{self.gatk_read_counts_folder}:{self.read_counts_volume}",
            f"{self.gatk_image}:{self.gatk_version}",
            "gatk", "GermlineCNVCaller",
            "-L", f"{Bed.volume}/{Bed.scatter_dirname}"
        ]
        input_read_files = self.get_input_read_count_files(Sample.sample_id, cohort_samples)
        cmd.extend(input_read_files)
        self.cohort_cnv_caller_filename = "cohort_caller"
        cmd2 = [
            "--contig-ploidy-calls", f"/model_dir/{self.ploidy_case_prefix}-calls"
            "--annotated_intervals", f"{Bed.volume}/{Bed.annotated_intevals_filename}",
            "--intervals-merging-rule", "OVERLAPPING_ONLY",
            "--output", f"/model_dir/",
            "--output-prefix", f"{self.cohort_cnv_caller_filename}",
            "--verbosity", "DEBUG"
        ]

        cmd.extend(cmd2)

        cmd_str = " ".join(cmd)
        logger.info(
            f"Running GermlineCNVCaller in cohort mode:\n{cmd_str}"
        )
        subprocess.run(cmd, encoding="utf-8", capture_output=True)

    def run_germline_cnv_caller_case_mode(self):
        self.caller_prefix = f"{self.sample.sample_id}_vs_cluster_cohort"
        if not hasattr(self, "model_dir"):
            self.model_dir = os.path.join(self.gatk_folder, self.sample.sample_id)
        cmd = [
            self.docker_path, "run",
            "-v", f"{self.model_dir}:/model_dir",
            "-v", f"{self.gatk_read_counts_folder}:{self.read_counts_volume}",
            f"{self.gatk_image}:{self.gatk_version}",
            "gatk", "GermlineCNVCaller",
            "--run-mode", "CASE",
            "-I", f"{self.read_counts_volume}/{self.sample.sample_name}.bam.hdf5",
            "--contig-ploidy-calls", f"/model_dir/{self.ploidy_case_prefix}",
            "--model", f"/model_dir/{self.cohort_cnv_caller_filename}-model",
            "--output", f"/model_dir/",
            "--output-prefix", f"{self.caller_prefix}",
            "--verbosity", "DEBUG"
        ]
        cmd_str = " ".join(cmd)
        logger.info(f"Running GermlineCNVCaller:\n{cmd_str}")
        subprocess.run(cmd, encoding="utf-8", capture_output=True)

    def run_postprocess_germline_calls(self):
        fasta_dir = os.path.dirname(self.reference_fasta)

        cmd = [
            self.docker_path, "run",
            "-v", f"{self.model_dir}:/model_dir",
            "-v", f"{self.gatk_results_dir}:/results_dir"
            "-v", f"{fasta_dir}:/ref_dir",
            f"{self.gatk_image}:{self.gatk_version}",
            "PostprocessGermlineCNVCalls",
            "--model-shard-path", f"/model_dir/{self.ploidy_prefix}-model",
            "--calls-shard-path", f"/model_dir/{self.ploidy_case_prefix}-calls",
            "--allosomal-contig", "chrX",
            "--allosomal-contig", "chrY",
            "--contig-ploidy-calls", f"/model_dir/{self.ploidy_case_prefix}-calls",
            "--sample-index", self.sample.sample_id,
            "--output-genotyped-intervals", f"/results_dir/{self.sample.sample_id}_intervals_cluster{self.sample.cluster}.vcf.gz",
            "--output-genotyped-segments", f"/results_dir/{self.sample.sample_id}_segments_cluster{self.sample.cluster}.vcf.gz",
            "--sequence-dictionary", f"/ref_dir/{self.dict_filename}"
        ]
        cmd_str = " ".join(cmd)
        logger.info(
            f"Running GATK PostprocessGermlineCNVCalls over sample {self.sample.sample_id}\n{cmd_str}"
        )

        subprocess.run(cmd, encoding="utf-8", capture_output=True)

        
# class Cohort_Gatk_gCNV(Gatk_gCNV)dict     def __init__(self, sample, docker_conf, reference_conf, Bed):
#         super().__init__(sample, docker_conf, reference_conf, Bed)
    
#     def run_determine_germline_contig_ploidy(self, Bed):

#         input_read_files = self.get_input_read_count_files(Sample.id, cohort_samples)
#         cmd = [
#             self.docker_path, "run",
#             "-v", f"{Bed.dir}:{Bed.volume}",
#             "-v", f"{self.gatk_folder}:{self.gatk_volume}",
#             "-v", f"{self.gatk_read_counts_folder}:{self.read_counts_volume}",
#             f"{self.gatk_image}:{self.gatk_version}",
#             "gatk", "DetermineGermlineContigPloidy",
#             "-L", f"{Bed.volume}/{Bed.annotated_intevals_filename}",
#             "--interval-merging-rule", "OVERLAPPING_ONLY",
#         ]
#         cmd.extend(input_read_files)

#         cmd2 = [
#             "--contig-ploidy-priors", f"{self.gatk_volume}/{self.contig_ploidy_filename}",
#             "--output", f"{self.gatk_volume}",
#             "--output-prefix", "ploidy_prefix",
#             "--verbosity", "DEBUG"
#         ]

#         cmd.extend(cmd2)

#         cmd_str = " ".join(cmd)

#         logger.info(
#             f"Running GATK DetermineGermlineContigPloidy:\n{cmd_str}"
#         )

#         subprocess.run(cmd, encoding="utf-8", capture_output=True)
