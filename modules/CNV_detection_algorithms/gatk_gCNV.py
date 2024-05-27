import os
import gzip
from modules.log import logger
from modules.CNV_detection_algorithms.CNV_algorithm import CNV_Algorithm
from modules.detected_CNV import CNV
class Gatk_gCNV(CNV_Algorithm):
    """
    Run all the steps of GATK gCNV in:
    https://gatk.broadinstitute.org/hc/en-us/articles/360035531152--How-to-Call-rare-germline-copy-number-variants
    """
    def __init__(self, docker_conf, reference_conf, Bed, force_run=False):
        super().__init__(docker_conf, reference_conf, Bed, force_run)

        self.gatk_image = docker_conf.gatk["image"]
        self.gatk_version = docker_conf.gatk["version"]
        self.gatk_folder = os.path.join(reference_conf.main_dir, "GATK")
        self.gatk_read_counts_dirname = "Read_counts"
        self.gatk_read_counts_folder = os.path.join(self.gatk_folder, self.gatk_read_counts_dirname)
        self.mappability_folder = os.path.join(self.gatk_folder, "mappability_track")
        self.mappability_track_path = os.path.join(self.mappability_folder, "k50.umap.bed")
        self.gatk_volume = "/gatk_vol"
        self.read_counts_volume = "/read_counts"
        self.contig_ploidy_filename = "contig_ploidy_prior.tsv"
        self.contig_ploidy = os.path.join(self.gatk_folder, self.contig_ploidy_filename)
        self.scatter_dirname = "scatter"
        self.scatter_path = os.path.join(Bed.dir, self.scatter_dirname)
        self.gatk_results_dir = os.path.join(self.results_dir, "GATK")
        if not os.path.exists(self.gatk_results_dir):
            os.mkdir(self.gatk_results_dir)
    
    def run_preprocess_intervals(self):
        
        fasta_dir = os.path.dirname(self.reference_fasta)
        fasta_filename = os.path.basename(self.reference_fasta)

        fasta_volume = "/fasta_dir"

        preprocessed_interval_list_filename = f"preprocessed_{self.Bed.interval_list_filename}"
        preprocessed_interval_list = os.path.join(self.Bed.dir, preprocessed_interval_list_filename)

        if os.path.exists(preprocessed_interval_list) and not self.force_run:
            logger.info(f"{preprocessed_interval_list} already exists, Gatk PreprocessIntervals won't be run")
            self.Bed.set_preprocessed_intervals_list_path(preprocessed_interval_list)
        cmd = [
            self.docker_path, "run",
            "-v", f"{fasta_dir}:{fasta_volume}",
            "-v", f"{self.Bed.dir}:{self.Bed.volume}",
            f"{self.gatk_image}:{self.gatk_version}",
            "gatk", "PreprocessIntervals",
            "-R", f"{fasta_volume}/{fasta_filename}",
            "-L", f"{self.Bed.volume}/{self.Bed.interval_list_filename}",
            "--padding", "25",
            "--bin-length", "0",
            "-imr", "OVERLAPPING_ONLY",
            "-O", f"{self.Bed.volume}/{preprocessed_interval_list_filename}"
        ]
        self.run_cmd(cmd, "GATK PreprocessIntervals")
        
        self.Bed.set_preprocessed_intervals_list_path(preprocessed_interval_list)
        return (preprocessed_interval_list)
    
    def create_gatk_read_counts_folder(self):


        if not os.path.isdir(self.gatk_read_counts_folder):
            logger.info(
                f"Creating new folder to store GATK gCNV files: {self.gatk_read_counts_folder}"
            )
            os.mkdir(self.gatk_read_counts_folder)
        return(self.gatk_read_counts_folder)
    
    
    def run_collect_read_counts(self, Sample):
        """
        Run GATK CollectReadCounts.
        """
        self.create_gatk_read_counts_folder()

        if not self.Bed.preprocessed_intervals_path:
            raise ValueError(
                f"You need to preprocess the interval list file first! Run the funciton Gatk_gCNV.run_preprocess_intervals()"
            )


        fasta_dir = os.path.dirname(self.reference_fasta)
        fasta_filename = os.path.basename(self.reference_fasta)


        hdf5_filename = f"{Sample.bam.filename}.hdf5"
        hdf5_path = os.path.join(self.gatk_read_counts_folder, hdf5_filename)
        
        if os.path.exists(hdf5_path) and not self.force_run:
            logger.info(
                f"GATK CollectReadCounts won't be executed as {hdf5_path} file already exists."
            )
            Sample.bam.set_hdf5_path_and_filename(hdf5_path)
            return(Sample.bam.hdf5_path)

        fasta_volume = "/fasta_dir"

        cmd = [
            self.docker_path, "run",
            "-v", f"{self.gatk_read_counts_folder}:/read_counts",
            "-v", f"{fasta_dir}:{fasta_volume}",
            "-v", f"{self.Bed.dir}:{self.Bed.volume}",
            "-v", f"{Sample.bam.dir}:{Sample.bam.volume}",
            f"{self.gatk_image}:{self.gatk_version}",
            "gatk", "CollectReadCounts",
            "-L", f"{self.Bed.volume}/{self.Bed.preprocessed_intervals_filename}",
            "-R", f"{fasta_volume}/{fasta_filename}",
            "-imr", "OVERLAPPING_ONLY",
            "-I", f"{Sample.bam.volume}/{Sample.bam.filename}",
            "--format", "HDF5",
            "-O", f"/read_counts/{hdf5_filename}"
        ]

        self.run_cmd(cmd, "GATK CollectReadCounts")

        Sample.bam.set_hdf5_path_and_filename(hdf5_path)
        return(Sample.bam.hdf5_path)
    
    def run_index_feature_file(self):
        """
        Indexing ./bed/mappability_track/k36.umap.bed.gz required to run annotate intervals
        """

        mappability_folder = os.path.basename(self.mappability_folder)
        mappability_filename = os.path.basename(self.mappability_track_path)

        index_map_track_path = f"{self.mappability_track_path}.idx"
        if os.path.isfile(index_map_track_path) and not self.force_run:
            return(index_map_track_path)
        cmd = [
            self.docker_path, "run",
            "-v", f"{self.Bed.dir}:{self.Bed.volume}",
            f"{self.gatk_image}:{self.gatk_version}",
            "gatk", "IndexFeatureFile",
            "-I", f"{self.Bed.volume}/{mappability_folder}/{mappability_filename}"
        ]
        self.run_cmd(cmd, "GATK IndexFeatureFile")
        return(index_map_track_path)


    def run_annotate_intervals(self):
        """
        Exclude problematic reagions based on GC content and mappability of intervals
        """

        if not self.Bed.preprocessed_intervals_path:
            raise ValueError(
                f"You need to preprocess the interval list file first! Run the funciton Gatk_gCNV.run_preprocess_intervals()"
            )
        
        fasta_dir = os.path.dirname(self.reference_fasta)
        fasta_filename = os.path.basename(self.reference_fasta)


        mappability_dirname = os.path.basename(self.mappability_folder)
        mappability_filename = os.path.basename(self.mappability_track_path)


        gc_annotated_bed_filename = f"annotated_{self.Bed.preprocessed_intervals_filename}"
        gc_annotated_bed_path = os.path.join(self.Bed.dir, gc_annotated_bed_filename)
        if os.path.exists(gc_annotated_bed_path) and not self.force_run:
            logger.info(
                f"GATK Annotate intevals won't be executed as {gc_annotated_bed_path} file already exists."
            )
            self.Bed.set_annotated_intervals_path(gc_annotated_bed_path)
            return(gc_annotated_bed_path)
        
        fasta_volume = "/fasta_dir"

        cmd = [
            self.docker_path, "run",
            "-v", f"{fasta_dir}:{fasta_volume}",
            "-v", f"{self.Bed.dir}:{self.Bed.volume}",
            "-v", f"{self.gatk_folder}:{self.gatk_volume}",
            f"{self.gatk_image}:{self.gatk_version}",
            "gatk", "AnnotateIntervals",
            "-L", f"{self.Bed.volume}/{self.Bed.preprocessed_intervals_filename}",
            "-R", f"{fasta_volume}/{fasta_filename}",
            "-imr", "OVERLAPPING_ONLY",
            "--mappability-track", f"{self.Bed.volume}/{mappability_dirname}/{mappability_filename}",
            "-O", f"{self.Bed.volume}/{gc_annotated_bed_filename}"
        ]
    
        self.run_cmd(cmd, "GATK AnnotateIntervals")
        self.Bed.set_annotated_intervals_path(gc_annotated_bed_path)

        return(gc_annotated_bed_path)
    
    def get_input_read_count_files(self, cohort_samples):
        """
        Get a list of input read counts files found in runs/GATK-gCNV7 to be given to gatk docker.
        E.g. ["-I", "gatk_vol/RB35645.hdf5", "-I", "gatk_vol/RB34532.hdf5]"""


        sample_hdf5_counts_filenames = [sample.bam.hdf5_filename for sample in cohort_samples]
        read_counts_input = list()
        for sample_hdf5_filename in sample_hdf5_counts_filenames:

            input_files_cmd = ["-I", f"{self.read_counts_volume}/{sample_hdf5_filename}"]
            read_counts_input.extend(input_files_cmd)
        
        return(read_counts_input)
    
    def run_filter_intervals(self, cohort_samples):
        """
        Given specific intervals (annotated intervals), and counts output by CollectReadCoutns,
        outputs a filtered Picard interval list
        """

        filtered_intervals_filename = f"filtered_{self.Bed.preprocessed_intervals_filename}"
        filtered_intervals_path = os.path.join(self.Bed.dir, filtered_intervals_filename)

        # it allways should be calculated as it depends on the sample analysed
        if os.path.isfile(filtered_intervals_path) and not self.force_run:
            logger.info(
                f"Filtered intervals file already exists: {filtered_intervals_path}, GATK FilteredIntervals won't be run"
            )
            self.Bed.set_filtered_intervals_path(filtered_intervals_path)
            print(self.Bed.filtered_intervals_filename)
            return(filtered_intervals_path)

        cmd = [
            self.docker_path, "run",
            "-v", f"{self.Bed.dir}:{self.Bed.volume}",
            "-v", f"{self.gatk_read_counts_folder}:{self.read_counts_volume}",
            f"{self.gatk_image}:{self.gatk_version}",
            "gatk", "FilterIntervals",
            "-L", f"{self.Bed.volume}/{self.Bed.preprocessed_intervals_filename}",
            "--annotated-intervals", f"{self.Bed.volume}/{self.Bed.annotated_intevals_filename}",
        ]

        read_counts_cmd = self.get_input_read_count_files(cohort_samples)
        
        cmd.extend(read_counts_cmd)

        cmd2 = [
            "-imr", "OVERLAPPING_ONLY",
            "-O", f"{self.Bed.volume}/{filtered_intervals_filename}",
        ]

        cmd.extend(cmd2)

        self.run_cmd(cmd, "GATK FilterIntervals")

        self.Bed.set_filtered_intervals_path(filtered_intervals_path)
        return(filtered_intervals_path)
  
    def run_interval_list_tools(self):
        cmd = [
            self.docker_path, "run",
            "-v", f"{self.Bed.dir}:{self.Bed.volume}",
            f"{self.gatk_image}:{self.gatk_version}",
            "gatk", "IntervalListTools",
            "--INPUT", f"{self.Bed.volume}/{self.Bed.filtered_intervals_filename}",
            "--SUBDIVISION_MODE", "INTERVAL_COUNT",
            "--SCATTER_CONTENT", "1000",
            "--OUTPUT", F"{self.Bed.volume}/{self.scatter_dirname}"
        ]

        self.run_cmd(cmd, "GATK IntervalListTools")

    def get_scatters(self):
        scatter_dirnames = os.listdir(self.scatter_path)
        scatter_paths = [os.path.join(self.scatter_path, scatter_dirname, "scattered.interval_list") for scatter_dirname in scatter_dirnames]
        return sorted(scatter_paths)


class Cohort_Gatk_gCNV(Gatk_gCNV):
    def __init__(self, docker_conf, reference_conf, Bed, cohort_samples, force_run=False):
        # Initialize attributes that belong to the parent class
        super().__init__(docker_conf, reference_conf, Bed, force_run)
        self.ploidy_prefix = "ploidy"
        self.model_dir = os.path.join(self.gatk_folder, "cohort")
        self.cohort_cnv_caller_filename = "cohort_caller"
        self.cohort_samples = cohort_samples

    
    def run_determine_germline_contig_ploidy(self):


        if not os.path.exists(self.model_dir):
            os.mkdir(self.model_dir)

        self.ploidy_model = os.path.join(self.model_dir, f"{self.ploidy_prefix}-model")
        if os.path.exists(self.ploidy_model) and not self.force_run:
            logger.info(
                f"DetermineGermlineContigPloidy already run in cohort mode for cohort samples"
            )
            return(True)
        
        
        input_read_files = self.get_input_read_count_files(self.cohort_samples)
        cmd = [
            self.docker_path, "run",
            "-v", f"{self.model_dir}:/model_dir", 
            "-v", f"{self.Bed.dir}:{self.Bed.volume}",
            "-v", f"{self.gatk_folder}:{self.gatk_volume}",
            "-v", f"{self.gatk_read_counts_folder}:{self.read_counts_volume}",
            f"{self.gatk_image}:{self.gatk_version}",
            "gatk", "DetermineGermlineContigPloidy",
            "-L", f"{self.Bed.volume}/{self.Bed.filtered_intervals_filename}",
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

        self.run_cmd(cmd, "GATK DetermineGermlineContigPloidy in cohort mode")



    def run_germline_cnv_caller(self):

        self.model = list()

        
        scatter_files = self.get_scatters()

        for i, scatter in enumerate(scatter_files):
            self.model.append(os.path.join(self.model_dir, f"{i}_{self.cohort_cnv_caller_filename}-model"))
            if os.path.exists(os.path.join(self.model_dir, f"{i}_{self.cohort_cnv_caller_filename}-model")) and not self.force_run:
                logger.info(
                    f"GATK GermlineCNVCaller already run in cohort mode: {self.model[i]}")
                continue
            scatter_dir = os.path.dirname(scatter)
            scatter_filename = os.path.basename(scatter)
            cmd = [
                self.docker_path, "run",
                "-v", f"{scatter_dir}:/scatter_dir",
                "-v", f"{self.Bed.dir}:{self.Bed.volume}",
                "-v", f"{self.model_dir}:/model_dir",
                "-v", f"{self.gatk_read_counts_folder}:{self.read_counts_volume}",
                f"{self.gatk_image}:{self.gatk_version}",
                "gatk", "GermlineCNVCaller",
                "--run-mode", "COHORT",
                "-L", f"/scatter_dir/{scatter_filename}"
            ]
            input_read_files = self.get_input_read_count_files(self.cohort_samples)
            cmd.extend(input_read_files)
            cmd2 = [
                "--contig-ploidy-calls", f"/model_dir/{self.ploidy_prefix}-calls",
                "--annotated-intervals", f"{self.Bed.volume}/{self.Bed.annotated_intevals_filename}",
                "--interval-merging-rule", "OVERLAPPING_ONLY",
                "--output", f"/model_dir/",
                "--output-prefix", f"{i}_{self.cohort_cnv_caller_filename}",
                "--verbosity", "DEBUG"
            ]

            cmd.extend(cmd2)

            self.run_cmd(cmd, "GATK GermlineCNVCaller in cohort mode")

class Case_Gatk_gCNV(Gatk_gCNV):
    def __init__(self, Cohort_gatk_gCNV, sample, docker_conf, reference_conf, Bed, force_run=False):
        self.sample = sample
        # Initialize attributes that belong to the parent class

        super().__init__(docker_conf, reference_conf, Bed, force_run)
        self.Cohort_Gatk = Cohort_gatk_gCNV
        self.ploidy_case_prefix = "ploidy-case"
        self.model_dir = os.path.join(self.gatk_folder, self.sample.sample_id)
        self.caller_prefix = f"{self.sample.sample_id}_vs_cluster_cohort"


    
    def run_determine_germline_contig_ploidy(self):

        if not os.path.exists(self.model_dir):
            os.mkdir(self.model_dir)

        output = os.path.join(self.model_dir, f"{self.ploidy_case_prefix}-calls")
        if os.path.exists(output) and not self.force_run:
            logger.info(
                f"DetermineGermlineContigPloidy already run in case mode for sample {self.sample.sample_id}"
            )
            return(True)
        cmd = [
            self.docker_path, "run",
            "-v", f"{self.gatk_folder}:{self.gatk_volume}",
            "-v", f"{self.model_dir}:/model_dir",
            "-v", f"{self.Cohort_Gatk.model_dir}:/cohort_model_dir",
            "-v", f"{self.gatk_read_counts_folder}:{self.read_counts_volume}",
            f"{self.gatk_image}:{self.gatk_version}",
            "gatk", "DetermineGermlineContigPloidy",
            "--model", f"/cohort_model_dir/{self.Cohort_Gatk.ploidy_prefix}-model",
            "-I", f"{self.read_counts_volume}/{self.sample.bam.hdf5_filename}",
            "-O", "/model_dir",
            "--output-prefix", self.ploidy_case_prefix,
            "--verbosity", "DEBUG"
        ]

        self.run_cmd(cmd, "GATK DetermineGermlineContigPloidy in case mode")


    def run_germline_cnv_caller(self):
        if not hasattr(self, "model_dir"):
            self.model_dir = os.path.join(self.gatk_folder, self.sample.sample_id)
        
        if os.path.exists(os.path.join(self.model_dir, f"{self.caller_prefix}-calls")) and not self.force_run:
            return True
        
        scatter_files = self.get_scatters()
        for i, scatter in enumerate(scatter_files):

            cmd = [
                self.docker_path, "run",
                "-v", f"{self.model_dir}:/model_dir",
                "-v", f"{self.Cohort_Gatk.model_dir}:/cohort_model_dir",
                "-v", f"{self.gatk_read_counts_folder}:{self.read_counts_volume}",
                f"{self.gatk_image}:{self.gatk_version}",
                "gatk", "GermlineCNVCaller",
                "--run-mode", "CASE",
                "-I", f"{self.read_counts_volume}/{self.sample.bam.hdf5_filename}",
                "--contig-ploidy-calls", f"/model_dir/{self.ploidy_case_prefix}-calls",
                "--model", f"/cohort_model_dir/{i}_{self.Cohort_Gatk.cohort_cnv_caller_filename}-model",
                "--output", f"/model_dir/",
                "--output-prefix", f"{i}_{self.caller_prefix}",
                "--verbosity", "DEBUG"
            ]
            self.run_cmd(cmd, "GATK GermlineCNVCaller in case mode")

    def run_postprocess_germline_calls(self):
        fasta_dir = os.path.dirname(self.reference_fasta)
        # sample index should be an integer
        sample_int = self.sample.sample_id.replace("RB", "")

        scatter_files = self.get_scatters()



        cmd = [
            self.docker_path, "run",
            "-v", f"{self.model_dir}:/model_dir",
            "-v", f"{self.Cohort_Gatk.model_dir}:/cohort_model_dir",
            "-v", f"{self.gatk_results_dir}:/results_dir",
            "-v", f"{fasta_dir}:/ref_dir",
            f"{self.gatk_image}:{self.gatk_version}",
            "gatk", "PostprocessGermlineCNVCalls",
        ] 
        model_shard = list()
        calls_shard = list()
        for i, scatter in enumerate(scatter_files):
            scatter_dir = os.path.dirname(scatter)
            scatter_filename = os.path.basename(scatter)
            model_shard.extend([
                "--model-shard-path", f"/cohort_model_dir/{i}_{self.Cohort_Gatk.cohort_cnv_caller_filename}-model"

            ])
            calls_shard.extend([
                "--calls-shard-path", f"/model_dir/{i}_{self.caller_prefix}-calls",

            ])
        cmd.extend(model_shard)
        cmd.extend(calls_shard)

        cmd2 = [
            # "--model-shard-path", f"/model_dir/{i}_{cohort_obj.cohort_cnv_caller_filename}-model",
            # "--calls-shard-path", f"/model_dir/{self.caller_prefix}-calls",
            "--allosomal-contig", "chrX",
            "--allosomal-contig", "chrY",
            "--contig-ploidy-calls", f"/model_dir/{self.ploidy_case_prefix}-calls",
            "--sample-index", "0", # in case mode is allways 0
            "--output-genotyped-intervals", f"/results_dir/{self.sample.sample_id}_intervals_cluster.vcf.gz",
            "--output-genotyped-segments", f"/results_dir/{self.sample.sample_id}_segments_cluster.vcf.gz",
            "--output-denoised-copy-ratios", f"/results_dir/{self.sample.sample_id}_denoised_copy_rations.hdf5",
            "--sequence-dictionary", f"/ref_dir/{self.dict_filename}"
        ]
        cmd.extend(cmd2)
        self.run_cmd(cmd, "GATK PostprocessGermlineCNVCalls in case mode")

    def process_cnvs(self):
        segments_file = os.path.join(self.gatk_results_dir, f"{self.sample.sample_id}_segments_cluster{self.sample.cluster}.vcf.gz")
        cnvs = list()
        with gzip.open(segments_file, "rt") as f:
            for line in f:
                if line.startswith("#"):
                    continue
                fields = line.strip().split("\t")
                cnv_type = fields[4]
                if cnv_type == ".":
                    # not a CNV segment
                    continue
                chr = fields[0]
                pos = fields[1]
                id = fields[2]
                ref = fields[3]
                qual = fields[5]
                filter = fields[6]
                end = fields[7].replace("END=", "")
                format = fields[8]
                other = fields[9]
                cnv = CNV(chr=chr, start=pos, end=end, type=cnv_type, algorithm="GATK_gCNV", qual=qual, sample=self.sample.sample_id)
                yield(cnv)
        
        return None

            


