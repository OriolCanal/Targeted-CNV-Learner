import os
from modules.bam import Bam
from modules.run_class import Analysis_Sample
def get_bams_bais(bam_dir):

    bams_filename = os.listdir(bam_dir)
    for bam_filename in bams_filename:
        if not bam_filename.endswith(".bam"):
            continue
        bam_path = os.path.join(bam_dir, bam_filename)

        Bam(bam_path)
        sample_name = bam_filename.split(".")[0]


    bams = [Bam(os.path.join(bam_dir, bam_filename)) for bam_filename in bams_filename if bam_filename.endswith(".bam")]
    for bam in bams:
        Analysis_Sample()
    return(bams)