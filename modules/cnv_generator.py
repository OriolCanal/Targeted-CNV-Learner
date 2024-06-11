import os
import subprocess
from modules.log import logger
from modules.bam import Bam
class CNV_Generator():

    def __init__(self, bam_dir, ref_conf, Bed_obj):
        self.analysis_bam_dir = bam_dir
        self.main_dir = ref_conf.main_dir
        self.cnvs_config_path = os.path.join(self.main_dir, "CNVs_config")
        if not os.path.exists(self.cnvs_config_path):
            os.mkdir(self.cnvs_config_path)
        self.fasta_path = ref_conf.hg19.fasta_path

        self.spike_in_bam_dirname = "SpikeInBAM"
        self.spike_in_bam_path = os.path.join(self.main_dir, self.spike_in_bam_dirname)
        self.filtered_config = os.path.join(self.cnvs_config_path, "filtered_cnvs.config")
        self.Bed_obj = Bed_obj

        # script that autogenerates CNVs config files
        self.generate_config_path = os.path.join(self.spike_in_bam_path, "scripts", "generate_config.py")

        self.configs_generated = list()

        self.bams_with_cnv_path = os.path.join(self.main_dir, "bams_with_cnvs")
        if not os.path.exists(self.bams_with_cnv_path):
            os.mkdir(self.bams_with_cnv_path)

    def create_multiple_exons_config(self, cnv_type, out_filename):
        allowed_types = ["del", "dup"]
        config_output_path = os.path.join(self.cnvs_config_path, out_filename)
        if cnv_type not in allowed_types:
            raise ValueError(
                f"Type {cnv_type} not in allowed types: {allowed_types}, "
            )

        if os.path.exists(config_output_path):
            logger.info(
                f"Multiple exons config file already created: {config_output_path}"
            )
            self.configs_generated.append(config_output_path)
            return(0)
        print(self.generate_config_path, self.analysis_bam_dir, cnv_type, self.Bed_obj.roi_bed)
        cmd = [
            "python", self.generate_config_path,
            "--indir", self.analysis_bam_dir,
            "-c", cnv_type,
            "-b", self.Bed_obj.roi_bed,
            "--mode", "multiple",
        ]
        str_cmd = " ".join(cmd)
        logger.info(
            f"Creating config for multiple exons CNVs: {config_output_path}\n{str_cmd}"
        )

        try:
            with open(config_output_path, 'w') as out_file:
                result = subprocess.run(cmd, stdout=out_file, stderr=subprocess.PIPE, text=True)
                
            if result.returncode != 0:
                logger.error(f"Error creating config for single exons CNVs: {result.stderr}")
                raise subprocess.CalledProcessError(result.returncode, cmd, output=result.stdout, stderr=result.stderr)
            
            self.configs_generated.append(config_output_path)
            logger.info(f"Config file created successfully: {config_output_path}")
        except Exception as e:
            logger.error(f"Failed to create config file: {e}")
            raise


    def create_single_exons_config(self, cnv_type, out_filename):
        allowed_types = ["del", "dup"]
        config_output_path = os.path.join(self.cnvs_config_path, out_filename)
        if cnv_type not in allowed_types:
            raise ValueError(
                f"Type {cnv_type} not in allowed types: {allowed_types}, "
            )
        if os.path.exists(config_output_path):
            logger.info(
                f"Single exons config file already created: {config_output_path}"
            )
            self.configs_generated.append(config_output_path)
            return(0)


        cmd = [
            "python", self.generate_config_path,
            "--indir", self.analysis_bam_dir,
            "-c", cnv_type,
            "-b", self.Bed_obj.roi_bed,
            "--mode", "single",
        ]

        str_cmd = " ".join(cmd)
        logger.info(
            f"Creating config for single exons CNVs: {config_output_path}:\n{str_cmd}"
        )
                
        try:
            with open(config_output_path, 'w') as out_file:
                result = subprocess.run(cmd, stdout=out_file, stderr=subprocess.PIPE, text=True)
                
            if result.returncode != 0:
                logger.error(f"Error creating config for single exons CNVs: {result.stderr}")
                raise subprocess.CalledProcessError(result.returncode, cmd, output=result.stdout, stderr=result.stderr)
            
            self.configs_generated.append(config_output_path)
            logger.info(f"Config file created successfully: {config_output_path}")
        except Exception as e:
            logger.error(f"Failed to create config file: {e}")
            raise


    def check_config_overlap(self):
        """
        Checks if 2 CNVs overlaps in the same gene in the same sample. If overlaps just includes 1 CNV in filtered.config
        
        Return:
        se_filtered_path: path of the config file without overlapping CNVs
        """
        sample_cnv = dict()
        logger.info(
            f"Creating filtered CNVs config with no overlapping CNV: {self.filtered_config}"
        )
        with open(self.filtered_config, "w") as outfile:
            for file in self.configs_generated:
                
                with open(file, "r") as f:
                    for str_line in f:
                        variant_overlap = False
                        line = str_line.strip().split("\t")
                        sample = os.path.basename(line[0])
                        cnv_gene = line[5].split("\t")[0]
                        if sample not in sample_cnv:
                            sample_cnv[sample] = [cnv_gene]
                            outfile.write(str_line)
                        else:
                            genes = sample_cnv[sample]
                            if cnv_gene not in genes:
                                outfile.write(str_line)
                                sample_cnv[sample].append(cnv_gene)
 
            
        return(self.filtered_config)
    
    def order_config(self, by_column=0):
        output_file = os.path.join(os.path.dirname(self.filtered_config), f"filtered_cnsv_ordered_by_col_{by_column}.config")
        # Read the lines of the input file into a list
        with open(self.filtered_config, 'r') as f:
            lines = f.readlines()
        
        # Sort the lines based on the values in the first column
        sorted_lines = sorted(lines, key=lambda line: line.split('\t')[by_column])
        
        # Write the sorted lines back to the output file
        with open(output_file, 'w') as f:
            f.writelines(sorted_lines)

    # def assign_in_silico_cnvs_to_sample(self, Sample_class):
    #     with open(self.filtered_config, "r") as f:
    #         for line in f:
    #             line = line.strip().split("\t")
    #             bam_path, chr, start, end, svtype, gene_exons = line[0], line[1], line[2], line[3], line[4], line[5]
    #             sample = os.path.basename(bam_path).split(".")[0]
    #             if "_" in gene_exons:
    #                 gene_exons = gene_exons.split("_")
    #                 gene = gene_exons[0]
    #                 exon_start = gene_exons[1]
    #                 exon_end = gene_exons[3]
    #                 numb_exons = exon_end - exon_start - 1
    #             sample_obj = Sample_class.sample_id_sample_obj[sample]
    #             sample_obj.cnvs["decon"].append(cnv)
                
    def generate_cnvs(self):
        
        script_path = os.path.join(self.spike_in_bam_path, "spikeinbam.py")
        cmd = [
            "python", script_path,
            "--variants", self.filtered_config,
            "--reference", self.fasta_path,
            "--threads", "4",
            "--suffix", ".simulated",
            "--output", self.bams_with_cnv_path
        ]
        str_cmd = " ".join(cmd)
        logger.info(
            f"Creating bams with in silico CNVs in {self.bams_with_cnv_path}:\n{str_cmd}"
        )
        if len(os.listdir(self.bams_with_cnv_path)) > 10:
            return (self.bams_with_cnv_path)
        subprocess.run(cmd)

        return self.bams_with_cnv_path
    
    
    # def sort_bams(self):
    #     bams = os.listdir(self.bams_with_cnv_path)
    #     bams_path = [os.path.join(self.bams_with_cnvs_path, bam) for bam in bams if bam.endswith(".bam")]
    #     for bam in bams_path:
    #         sorted_bam = bam.replace(".bam", ".sorted.bam")
    #         cmd = [
    #             "samtools", "sort", bam, "-o", sorted_bam
    #         ]
    #         cmd_str = " ".join(cmd)
    #         logger.info(
    #             f"sorting bam file:\n"
    #         )
    #         subprocess.run(cmd)
        
    def create_bam_objects(self):

        pass
                