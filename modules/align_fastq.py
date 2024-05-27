import os
import subprocess

class BaseMapper:
    """
    Parent class to handle read alignment processes.
    This class checks that input fastq files are valid
    :param str sample_name: sample name extracted from fastq
    :param str fq1: raw fastq1 (R1)
    :param str fq2: raw fastq2 (R2)
    :param str ref: reference genome in fasta format
    :param str outdir: output directory
    :param bin_config: configuration object for binaries
    """

    def __init__(self, sample_name, fq1, fq2, ref, mem2, outdir, docker_config):
        self._sample_name = sample_name
        self._fq1 = fq1
        self._fq2 = fq2
        # Reference genome in FASTA
        self._ref = ref
        self._mem2 = mem2
        # Output BAM
        self._outdir = outdir
        self._docker_config = docker_config

        # Now check that both FASTQ files are valid
        try:
            Fastq(self._fq1, expect_paired=True)
        except Exception as e:
            logging.error(e.message)


class Bwa(BaseMapper):
    """
     BWA-MEM mapper class
    :param str sample_name: sample name extracted from fastq
    :param str fq1: raw fastq1 (R1)
    :param str fq2: raw fastq2 (R2)
    :param str ref: reference genome in fasta format
    :param str outdir: output directory
    :param bin_config: configuration object for binaries
    """

    def __init__(self, sample_name, fq1, fq2, ref, mem2, outdir, docker_config):
        # Inherits variant instances from Aligner class
        super().__init__(sample_name, fq1, fq2, ref, mem2, outdir, docker_config)
        # super().__init__(self)

    def align_mem2(self, num_cpu=4, sort=True):
        """
         Align with BWA-MEM2
        :param int num_cpu: number of CPU's
        :param bool sort: sort output BAM by coordinates
        """

        bam_name = f"{self._sample_name}.bam"
        bam_out = str(Path(self._outdir) / bam_name)

        samtools_cpus = 4
        bwamem_cpus = num_cpu-2


        cmd = ""
        if sort is True:
            cmd = (
                "{} run -v {}:/fastq_dir/ -v {}:/ref_dir/ -v {}:/bam_dir/ "
                " {} /bin/bash -c \"/opt/bwa-mem2-2.2.1_x64-linux/bwa-mem2.avx mem /ref_dir/{} -R '@RG\\tID:{}\\tSM:{}'"
                " -M -t {} /fastq_dir/{}"
                " /fastq_dir/{} | samtools view -Shu - "
                ' | samtools sort -@ {} -T /bam_dir/{} -o /bam_dir/{}" '
            ).format(
                self._docker_config.docker,
                os.path.dirname(self._fq1),
                os.path.dirname(self._ref),
                self._outdir,
                self._docker_config.ngs_binaries.image,
                os.path.basename(self._mem2),
                self._sample_name,
                self._sample_name,
                bwamem_cpus,
                os.path.basename(self._fq1),
                os.path.basename(self._fq2),
                samtools_cpus,
                self._sample_name,
                bam_name,
            )
        if not os.path.isfile(bam_out):
            msg = f" INFO: Mapping sample {self._sample_name}"
            logging.info(msg)
            logging.info(cmd)
            p1 = subprocess.run(
                cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
            )
            error = p1.stderr.decode("UTF-8")
            if p1.returncode != 0:
                msg = f" ERROR: Could map sample {self._sample_name} using bwa mem"
                logging.error(msg)
                logging.error(error)
                sys.exit()
        else:
            msg = (
                f" INFO: Skipping mapping for sample {self._sample_name}."
                " Bam is already available"
            )
            logging.info(msg)


    def align(self, num_cpu=4, sort=True):
        """
         Align with BWA-MEM
        :param int num_cpu: number of CPU's
        :param bool sort: sort output BAM by coordinates
        """

        bam_name = f"{self._sample_name}.bam"
        bam_out = str(Path(self._outdir) / bam_name)

        cmd = ""
        if sort is True:
            cmd = (
                "{} run -v {}:/fastq_dir/ -v {}:/ref_dir/ -v {}:/bam_dir/ "
                " {} /bin/bash -c \"bwa mem /ref_dir/{} -R '@RG\\tID:{}\\tSM:{}'"
                " -M -t {} /fastq_dir/{}"
                " /fastq_dir/{} | samtools view -Shu - "
                ' | samtools sort -T /bam_dir/{} -o /bam_dir/{}" '
            ).format(
                self._docker_config.docker,
                os.path.dirname(self._fq1),
                os.path.dirname(self._ref),
                self._outdir,
                self._docker_config.ngs_binaries.image,
                os.path.basename(self._ref),
                self._sample_name,
                self._sample_name,
                num_cpu,
                os.path.basename(self._fq1),
                os.path.basename(self._fq2),
                self._sample_name,
                bam_name,
            )
        if not os.path.isfile(bam_out):
            msg = f" INFO: Mapping sample {self._sample_name}"
            logging.info(msg)
            logging.info(cmd)
            p1 = subprocess.run(
                cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
            )
            error = p1.stderr.decode("UTF-8")
            if p1.returncode != 0:
                msg = f" ERROR: Could map sample {self._sample_name} using bwa mem"
                logging.error(msg)
                logging.error(error)
                sys.exit()
        else:
            msg = (
                f" INFO: Skipping mapping for sample {self._sample_name}."
                " Bam is already available"
            )
            logging.info(msg)


class Fastq():
    def __init__(self, path):
        self.path = path

    