# USED WHEN BAMS WERE IMPORTED FROM THE DB
import os
import subprocess
import shutil
import matplotlib.pyplot as plt

from modules.log import logger
from modules.mongo_classes import Sample
from modules.bam import Bam
class Analysis_Run():
    samples_analysed = set()
    duplicated_samples = set()
    def __init__(self, run_path, is_cohort=False):
        self.run_path = run_path
        self.run_id = os.path.basename(run_path)
        self.samples_147 = list()
        self.is_cohort = is_cohort
    
    def put_bams_in_cohort_dir(self, ref_conf):
        """
        Create a cohort directory where there will be stored all the cohort bam files
        """
        if self.is_cohort is True:
            cohort_dir = os.path.join(ref_conf.main_dir, "cohort_bams")
        else:
            logger.warning(
                f"Trying to put bams of run {self.run_id} into cohort_dir, but this run has not \
                 been specified to be a cohort run!")
            return(False)

        if not os.path.exists(cohort_dir):
            os.mkdir(cohort_dir)
        
        logger.info(
            f"Copying bams of RUN: {self.run_id} to {cohort_dir}"
            )
        for sample in self.samples_147:
            out_path = os.path.join(cohort_dir, sample.bam.filename)
            if os.path.exists(out_path) and os.path.getsize(out_path) > 0:
                continue
            cmd = [
                "cp", sample.bam.path, out_path
            ]
            subprocess.run(cmd)
            out_path_bai = os.path.join(cohort_dir, sample.bam.bai_filename)
            cmd_bai = [
                "cp", sample.bam.bai_path, out_path_bai
            ]
            subprocess.run(cmd_bai)

        
        Bam.set_cohort_bam_dir(cohort_dir)
        # self.run_path = 
    def put_bams_in_analysis_dir(self, ref_conf):

        if self.is_cohort is False:
            analysis_dir = os.path.join(ref_conf.main_dir, "analysis_bams")
        else:
            logger.warning(
                f"You can't put bam files of a cohort run into analysis directory,"
            )
            return(False)
        if not os.path.exists(analysis_dir):
            os.mkdir(analysis_dir)

        logger.info(
            f"Copying bams of RUN: {self.run_id} to {analysis_dir}"
        )

        for sample in self.samples_147:
            out_path = os.path.join(analysis_dir, sample.bam.filename)
            if os.path.exists(out_path) and os.path.getsize(out_path) > 0:
                continue
            cmd = [
                "cp", sample.bam.path, out_path
            ]
            subprocess.run(cmd)
            out_path_bai = os.path.join(analysis_dir, sample.bam.bai_filename)
            cmd_bai = [
                "cp", sample.bam.bai_path, out_path_bai
            ]
            subprocess.run(cmd_bai)

        
        Bam.set_analysis_bam_dir(analysis_dir)
    
    def put_bam_in_run_directory(self):

        for sample in self.samples_147:
            bam_obj = sample.bam
            bam_path = bam_obj.path
            bam_dir = bam_obj.dir
            # if already in run directory, don't change it
            # if bam_obj.dir.endswith(self.run_id):
            #     continue
            # if run id already is the last directory,
            if os.path.basename(bam_dir) == self.run_id:
                # bam_obj.set_bam_new_path(new_bam_path)
                self.run_path = bam_dir
                return 0
            
            new_dir = os.path.join(bam_obj.dir, self.run_id)
            if not os.path.exists(new_dir):
                os.mkdir(new_dir)
            
            new_bam_path = os.path.join(new_dir, bam_obj.filename)
            if os.path.exists(new_bam_path):
                continue
            cmd = [
                "mv", bam_path, new_bam_path
            ]
            cmd2 = [
                "mv", f"{bam_path}.bai", f"{new_bam_path}.bai"
            ]
            cmd_str = " ".join(cmd)
            logger.info(
                f"Moving bam to run directory:\n{cmd_str}"
            )
            subprocess.run(cmd)
            subprocess.run(cmd2)
            bam_obj.set_bam_new_path(new_bam_path)
            self.run_path = new_dir

    def create_symbolic_link(self, ref_conf):
        if self.is_cohort:
            self.symbolic_link_dir = os.path.join(ref_conf.main_dir, "cohort_symbolic_link")
            if not os.path.exists(self.symbolic_link_dir):
                os.mkdir(self.symbolic_link_dir)
        if not self.is_cohort:
            self.symbolic_link_dir = os.path.join(ref_conf.main_dir, "analysis_symbolic_link")
            if not os.path.exists(self.symbolic_link_dir):
                os.mkdir(self.symbolic_link_dir)
        
        for sample in self.samples_147:

            destination_bam_path = os.path.join(self.symbolic_link_dir, sample.bam.filename)
            destination_bai_path = os.path.join(self.symbolic_link_dir, sample.bam.bai_filename)
            if not os.path.exists(destination_bam_path):
                os.symlink(sample.bam.path, destination_bam_path)
                os.symlink(sample.bam.bai_path, destination_bai_path)

            bam = sample.bam
            bam.set_symbolic_link(destination_bam_path)

    # def get_samples_from_db(self, Run_mongo):
    #     try:
    #         run_doc = Run_mongo.objects.get(run_id=self.run_id)
    #     except:
    #         logger.warning(
    #             f"Run ID: {self.run_id} not found in MongoDb!"
    #         )
        
    #     if hasattr(run_doc, "samples"):
    #         samples = run_doc.samples
    #         for sample_id in samples:
    #             Sample = Analysis_Sample(run_id=self.run_id, sample_id=sample_id)
    #             self.samples_ids.append(Sample)
    #             if Sample.is_panel_147(Sample_mongo):
    #                 self.samples_147.append(Sample)

    #     else:
    #         raise (ValueError(
    #             f"Run ID: {self.run_id} does not contain any sample"
    #         ))

    def remove_sample(self, Sample):
        if Sample in self.samples_147:
            self.samples_147.remove(Sample)
        else:
            raise ValueError(
                f" Sample: {Sample.sample_id} not in run {self.run_id}"
            )


    def get_Bams_and_Samples(self):
        bams_bais = os.listdir(self.run_path)
        bams = [bam for bam in bams_bais if bam.endswith(".bam")]

        for bam in bams:
            # avoiding duplicated samples
            if bam not in self.samples_analysed:
                self.samples_analysed.add(bam)
            else:
                logger.info(
                    f"Bam: {bam} already analysed"
                )
                continue
                
            bam_path = os.path.join(self.run_path, bam)
            self.duplicated_samples.add(bam_path)
            Bam_class = Bam(bam_path)
            sample_name = bam.split(".")[0]
            sample = Sample(Bam=Bam_class, sample_id=sample_name, run_id=self.run_id, is_cohort=self.is_cohort)
            self.samples_147.append(sample)

class Sample():
    compendi_bai_path = "http://172.16.83.24:8001/download_sample_bai/"
    compendi_bam_path = "http://172.16.83.24:8001/download_sample_bam/"

    sample_id_sample_obj = dict()
    def __init__(self, Bam, sample_id=None, run_id=None, is_cohort=False):
        self.run_id = run_id
        self.sample_id = sample_id
        self.downloaded_bam = False
        self.bam = Bam # Bam class containing info about bam files
        self.cluster = None
        self.is_outlier = False
        self.is_cohort = is_cohort
        self.cnv_kit_vcf_path = None
        self.cnv_kit_vcf_filename = None
        self.gatk_vcf_filename = None
        self.gatk_vcf_path = None
        self.cnvs = {
            "gatk": list(),
            "decon": list(),
            "cnvkit": list(),
            "grapes2": list(),
            "in_silico": list()
        }
        Sample.sample_id_sample_obj[self.sample_id] = self


    def get_all_callers_cnvs(self):
        all_cnvs = list()
        for key, value in self.cnvs.items():
            if key == "in_silico":
                continue
            all_cnvs.extend(value)
        
        return(all_cnvs)
    
    def check_sample_overlap_cnv(self):
        threshold = 15
        for in_silico_cnv in self.cnvs["in_silico"]:
            for decon_cnv in self.cnvs["decon"]:
                if in_silico_cnv.chr != decon_cnv.chr:
                    continue
                
                if in_silico_cnv.start > decon_cnv.end or decon_cnv.start > in_silico_cnv.end:
                    continue

                if in_silico_cnv.end < decon_cnv.start or decon_cnv.end < in_silico_cnv.start:
                    continue

                in_silico_overlap, decon_overlap = in_silico_cnv.calculate_overlap_percentage(decon_cnv)
                in_silico_cnv.algorithms_overlap["decon"] = in_silico_overlap
                decon_cnv.overlap_percentage = decon_overlap
                if float(in_silico_overlap) > threshold and float(decon_overlap) > threshold:
                    in_silico_cnv.algorithms_detected.add("decon")
                    decon_cnv.set_in_silico_cnv(in_silico_cnv)

            for grapes_cnv in self.cnvs["grapes2"]:
                if in_silico_cnv.chr != grapes_cnv.chr:
                    continue
                
                if in_silico_cnv.start > grapes_cnv.end or grapes_cnv.start > in_silico_cnv.end:
                    continue

                if in_silico_cnv.end < grapes_cnv.start or grapes_cnv.end < in_silico_cnv.start:
                    continue

                in_silico_overlap, grapes_overlap = in_silico_cnv.calculate_overlap_percentage(grapes_cnv)
                in_silico_cnv.algorithms_overlap["grapes2"] = in_silico_overlap
                grapes_cnv.overlap_percentage = grapes_overlap
                if float(in_silico_overlap) > threshold and float(grapes_overlap) > threshold:
                    in_silico_cnv.algorithms_detected.add("grapes2")
                    grapes_cnv.overlap = grapes_overlap
                    grapes_cnv.set_in_silico_cnv(in_silico_cnv)
                    
                

            for cnvkit_cnv in self.cnvs["cnvkit"]:
                if in_silico_cnv.chr != cnvkit_cnv.chr:
                    continue
                
                if in_silico_cnv.start > cnvkit_cnv.end or cnvkit_cnv.start > in_silico_cnv.end:
                    continue

                if in_silico_cnv.end < cnvkit_cnv.start or cnvkit_cnv.end < in_silico_cnv.start:
                    continue

                in_silico_overlap, cnvkit_overlap = in_silico_cnv.calculate_overlap_percentage(cnvkit_cnv)
                in_silico_cnv.algorithms_overlap["cnvkit"] = in_silico_overlap
                cnvkit_cnv.overlap_percentage = cnvkit_overlap
                if float(in_silico_overlap) > threshold and float(cnvkit_overlap) > threshold:
                    in_silico_cnv.algorithms_detected.add("cnvkit")
                    cnvkit_cnv.set_in_silico_cnv(in_silico_cnv)
    
    
    # def plot_overlap_histograms(self, in_silico_overlaps, algorithm_overlaps, algorithm_name="Algorithm", ref_conf=None):

    #     plt.figure(figsize=(10, 6))

    #     plt.hist(in_silico_overlaps, bins=20, alpha=0.5, label='In Silico Overlap')
    #     plt.hist(algorithm_overlaps, bins=20, alpha=0.5, label=f'{algorithm_name} Overlap')

    #     plt.xlabel('Overlap Percentage')
    #     plt.ylabel('Frequency')
    #     plt.title('Histogram of Overlap Percentages')
    #     plt.legend(loc='upper right')
    #     plt.grid(True)
    #     if ref_conf:
    #         plots_dir =  os.path.join(ref_conf.main_dir, "Plots", "overlap_cnv")
    #         if not os.path.exists(plots_dir):
    #             os.mkdir(plots_dir)
    #         plot_path = os.path.join(plots_dir, f"In_silico_vs_{algorithm_name}_overlap.png")
    #         plt.savefig(plot_path)
    #     plt.show()

        
    def set_cnv_kit_vcf(self, cnv_kit_vcf_path):
        if os.path.isfile(cnv_kit_vcf_path):
            self.cnv_kit_vcf_path = cnv_kit_vcf_path
        else:
            raise ValueError(
                f"Hdf5 path does not exist: {cnv_kit_vcf_path}"
            )
        self.cnv_kit_vcf_filename = os.path.basename(cnv_kit_vcf_path)

    def set_gatk_vcf(self, gatk_vcf_path):
        if os.path.isfile(gatk_vcf_path):
            self.gatk_vcf_path = gatk_vcf_path
        else:
            raise ValueError(
                f"Hdf5 path does not exist: {gatk_vcf_path}"
            )
        self.gatk_vcf_filename = os.path.basename(gatk_vcf_path)

# class Analysis_Sample():
#     compendi_bai_path = "http://172.16.83.24:8001/download_sample_bai/"
#     compendi_bam_path = "http://172.16.83.24:8001/download_sample_bam/"

#     def __init__(self, Bam, sample_id=None, run_id=None):
#         self.run_id = run_id
#         self.sample_id = sample_id
#         self.downloaded_bam = False
#         self.bam = Bam # Bam class containing info about bam files
#         self.assigned_cluster = None
#         self.is_outlier = False


# class Analysis_Sample():
#     compendi_bai_path = "http://172.16.83.24:8001/download_sample_bai/"
#     compendi_bam_path = "http://172.16.83.24:8001/download_sample_bam/"

#     def __init__(self, sample_id, Bam=None, run_id=None):
#         self.run_id = run_id
#         self.sample_id = sample_id
#         self.downloaded_bam = False
#         self.bam = Bam # Bam class containing info about bam files
#         self.cluster = None
#     def is_panel_147(self, Sample_mongo):
#         try:
#             sample_doc = Sample_mongo.objects.get(run_id=self.run_id, lab_id=self.sample_id)
#         except:
#             logger.warning(
#                 f"Sample ID: {self.sample_id} not found in MongoDb Sample collection!"
#             )
#             return False
        
#         if hasattr(sample_doc, "panel"):
#             if "SUDD_147" in sample_doc.panel:
#                 self.panel_147 = True
#                 return(True)
#             elif sample_doc.panel == "147":
#                 self.panel_147 = True
#                 return(True)
#             self.panel_147 = False
#             return False
        
#         raise ValueError(
#             f"Sample document with ID {self.sample_id} does not contain a panel attribute!"
#         )
    
            
#     def get_bam_bai_from_compendi(self, ref_conf):
#         if not self.panel_147:
#             logger.warning(f"Bam file won't be downloaded as it's not panel SUDD_147")
#             return(0)
#         output_dir = os.path.join(ref_conf.main_dir, "runs")

#         if not os.path.isdir(output_dir):
#             logger.info(f"Creating directory {output_dir} to store bam bai files.")
#             os.mkdir(output_dir)

#         run_dir = os.path.join(output_dir, self.run_id)
#         if not os.path.isdir(run_dir):
#             os.mkdir(run_dir)
#         sample_bam = f"{self.sample_id}.bam"
#         sample_bai = f"{self.sample_id}.bai"
#         bam_path = os.path.join(run_dir, sample_bam)
#         bai_path = os.path.join(run_dir, sample_bai)
#         if not os.path.isfile(bam_path):
#             download_bam_compendi = os.path.join(Analysis_Sample.compendi_bam_path, self.sample_id, self.run_id)

#             result = subprocess.run(["wget", "-O", bam_path, download_bam_compendi])

#             if result.returncode == 0:
#                 logger.info(f"{bam_path} downloaded successfully!")
#                 self.bam = Bam(bam_path)
#             else:
#                 logger.error(
#                     f"Download from {download_bam_compendi} failed with return code: {result.returncode}"
#                 )
#                 if os.path.isfile(bam_path):
#                     shutil.rmtree(run_dir)
#                 return(False)
#         else:
#             self.bam = Bam(bam_path)
#             logger.info(
#                 f"Bam file: {bam_path} has previously been downloaded!"
#             )

#         if not os.path.isfile(bai_path):
#             download_bai_compendi = os.path.join(Analysis_Sample.compendi_bai_path, self.sample_id, self.run_id)

#             result_bairef_config.bam_dir) = subprocess.run(["wget", "-O", bai_path, download_bai_compendi])

#             if result_bai.returncode == 0:
#                 logger.info(f"{bai_path} downloaded successfully!")
#                 self.bam = Bam(bam_path)
#             else:
#                 logger.error(
#                     f"Download from {download_bai_compendi} failed with return code: {result.returncode}"
#                 )
#                 if os.path.isfile(bai_path):
#                     shutil.rmtree(run_dir)
#                 return(False)
#         else:
#             logger.info(
#                 f"Bai file: {bai_path} has previously been downloaded!"
#             )

#         return(self.bam)
        
