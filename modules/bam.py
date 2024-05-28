import os

class Bam():
    def __init__(self, bam_path, hs_metrics_path="", hdf5_path=""):
        self.path = bam_path
        self.filename = os.path.basename(self.path)
        self.dir = os.path.dirname(self.path)

        self.bai_path = f"{self.path}.bai"
        self.bai_filename = os.path.basename(self.bai_path)
        self.validate_bam_bai()

        self.symbolic_link = None # symbolic link to have all cohort/analysis samples in corresponding directory

        self.volume = "/bam_dir"
        self.gatk_dirname = "GATK_gCNV"
        self.gatk_dir = os.path.join(self.dir, self.gatk_dirname)
        
        self.hdf5_path = hdf5_path
        self.hdf5_filename = os.path.basename(self.hdf5_path)

        self.hs_metrics_path = hs_metrics_path
        self.hs_metrics_filename = os.path.basename(self.hs_metrics_path)

    def set_symbolic_link(self, symbolic_link):
        self.symbolic_link = symbolic_link
    def validate_bam_bai(self):
        bam_size = os.path.getsize(self.path)
        bai_size = os.path.getsize(self.bai_path)
        if os.path.exists(self.path) and os.path.exists(self.bai_path):
            if bam_size > 0 and bai_size > 0:
                return True
        
        raise ValueError(
            f"Bam and bai files are not present or its size is 0: \n\t{self.path}\n\t{self.bai_path}"
        )
    

    def set_hdf5_path_and_filename(self, hdf5_path):
        if os.path.isfile(hdf5_path):
            self.hdf5_path = hdf5_path
        else:
            raise ValueError(
                f"Hdf5 path does not exist: {hdf5_path}"
            )
        self.hdf5_filename = os.path.basename(hdf5_path)

    
    def set_hs_metrics(self, hs_metrics_path):
        if os.path.isfile(hs_metrics_path):
            self.hs_metrics_path = hs_metrics_path
        else:
            raise ValueError(
                f"Hs_metrics path does not exist: {hs_metrics_path}"
            )
        self.hs_metrics_filename = os.path.basename(hs_metrics_path)

    