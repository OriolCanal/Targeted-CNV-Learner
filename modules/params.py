
import os
import subprocess
import sys
import tempfile
from datetime import datetime
from pathlib import Path
import pandas as pd
import yaml

from modules.log import logger


class YamlConfig:
    """
    Base config class that sets default parameters
    """

    def __init__(self, config_file):
        # Setting defaults
        self.config_file = config_file
        config_yaml = self.yaml_to_dict(self.config_file)
        self.main_dir = config_yaml["main_dir"]
        self.mongo_restore_file = os.path.join( self.main_dir, config_yaml["mongo_restore_file"])
        self.docker_yaml = os.path.join(self.main_dir, config_yaml["docker_yaml_file"])
        # The following params can be overwritten with an yaml
        self.yaml_dir = os.path.join(self.main_dir, "yaml_files")
        self.reference_yaml = os.path.join(self.main_dir, config_yaml["reference_yaml_file"])
        self.annotations_yaml = os.path.join(self.main_dir, config_yaml["annotations_yaml_file"])
        self.isoform_path = os.path.join(self.main_dir, config_yaml["isoforms_to_annotate"])
        self.toy_vcf = os.path.join(self.main_dir, config_yaml["mock_vcf"])
        self.bed = os.path.join(self.main_dir, config_yaml["bed_file"])
    def yaml_to_dict(self, yaml_file) -> str:
        """
        Take a valid yaml file and return a dict. Sanity check
        """
        if yaml_file is None:
            msg = "missing yaml_file param"
            raise FileNotFoundError(msg)

        try:
            file_path = Path(yaml_file)
            file_path.is_file() is True
            # Path(yaml_file)
        except NameError:
            msg = f" ERROR: Invalid input yaml file {yaml_file}"
        else:
            try:
                Path(yaml_file).exists() is True
            except NameError:
                msg = f" ERROR: missing input yaml file {yaml_file}"

        with open(yaml_file, "r") as yaml_f:
            yaml_dict = yaml.safe_load(yaml_f)
        return yaml_dict

class AnnotationConfig(YamlConfig):
    """
    Configuration for annotation resources
    """

    def __init__(self, config_file):
        super().__init__(config_file)
        self._data = self.yaml_to_dict(self.annotations_yaml)
        if "ann_dir" in self._data:
            setattr(self, "_ann_dir", self._data["ann_dir"])

    def parse_data(self):
        # print(self._data)
        for resource, values in self._data.items():
            # print(values)
            if isinstance(values, (str, int, float, list, tuple)):
                setattr(
                    self,
                    resource,
                    values
                )
            elif isinstance(values, dict):
                # print(resource)
                setattr(
                    self,
                    resource,
                    AnnotationDatabase(self.config_file, values),
                )


class AnnotationDatabase(AnnotationConfig):
        def __init__(self, config_file, resource_data):
            super().__init__(config_file)
            self.resource_data = resource_data
            self.parse_resource_data()

        def parse_resource_data(self):
            if not isinstance(self.resource_data, dict):
                dict_example = {
                    'version': 1.68,
                    'resurce_name': 'yaml',
                    'dirname': 'yaml',
                    'file': 'hg38/annotation_resources_v1.68.yaml',
                    'last_review': '20220411',
                    'md5': '2d55cc751e662471fef7e5b0b2960c93'
                }
                msg = f"Input of AnnotationDatabase should be a dictionary: e.g.: {dict_example}\n{self.resource_data} was given"
                raise ValueError(msg)
            for field, item in self.resource_data.items():
                # print(field)

                setattr(
                    self,
                    str(field),
                    item
                )
                if field == "file":
                    # if no file specified, don't parse it
                    if not item:
                        full_path = os.path.join(self._ann_dir, self.dirname)
                        ann_relative_path = os.path.join(self.dirname)
                    else:
                        full_path = os.path.join(self._ann_dir, self.dirname, item)
                        ann_relative_path = os.path.join(self.dirname, item)
                    setattr(
                        self,
                        "file_path",
                        full_path
                    )
                    setattr(
                        self,
                        "ann_relative_path",
                        ann_relative_path
                    )
                if field == "dirname":
                    dir_path = os.path.join(self._ann_dir, item)
                    setattr(
                        self,
                        "dir_path",
                        dir_path
                    )


  
class ReferenceConfig(YamlConfig):
    def __init__(self, config_file):
        super().__init__(config_file)

        self._data = self.yaml_to_dict(self.reference_yaml)
        if "ref_dir" in self._data:
            setattr(self, "_ref_dir", self._data["ref_dir"])

    def parse_data(self):
        for resource, values in self._data.items():
            if isinstance(values, (str, int, float, list, tuple)):
                setattr(
                    self,
                    resource,
                    values
                )
            elif isinstance(values, dict):
                # print(resource)
                setattr(
                    self,
                    resource,
                    ReferenceGenome(self.config_file, values),
                )
    
    def get_genome_data(self, genome_version:str):
        if not hasattr(self, genome_version):
            raise ValueError(f"genome version {genome_version} not specified in {self.reference_yaml}")
        
        version_instance = getattr(self, genome_version)

        return (version_instance)

class ReferenceGenome(ReferenceConfig):
    def __init__(self, config_file, genome_data):
        super().__init__(config_file)
        self.genome_data = genome_data
        self.parse_genome_data()

    def parse_genome_data(self):
        if not isinstance(self.genome_data, dict):
            dict_example = {
                'driname': "",
                'fasta': 'hg38.analysisSet.fasta',
                "dict": "hg38.analysisSet.dict",
                "gene_bed": "genelist.hg38.bed.gz",
                "chrom_sizes": ""
            }
            msg = f"Input of ReferenceGenome should be a dictionary: e.g.: {dict_example}\n{self.genome_data} was given"
            raise ValueError(msg)
        for field, item in self.genome_data.items():
            setattr(
                self,
                str(field),
                str(item)
            )
            if field == "dirname":
                dir_path = os.path.join(self._ref_dir, item)
                setattr(
                    self,
                    "dir_path",
                    dir_path
                )
            if field == "fasta":
                full_path = os.path.join(self._ref_dir, self.dirname, item)
                
                setattr(
                    self,
                    "fasta_path",
                    full_path
                )
            if field == "dict":
                dict_full_path = os.path.join(self._ref_dir, self.dirname, item)
                
                setattr(
                    self,
                    "fasta_dict_path",
                    dict_full_path
                )

class DockerConfig(YamlConfig):
    """
    Class that sets docker image
    """

    def __init__(self,config_file):
        super().__init__(config_file)
        self._data = self.yaml_to_dict(self.docker_yaml)

        self.validate()
        if self._data is not None:
            for a, b in self._data.items():
                setattr(
                    self,
                    str(a),
                    b
                )
    @property
    def docker(self):
        try:
            self.get_bin_path("docker")
        except Exception:
            msg = " ERROR: docker was not found on PATH"
            logger.error(msg)
            raise Exception(msg)
        else:
            return self.get_bin_path("docker")

    def get_bin_path(self, program):
        """
        Get the PATH of a program
        """
        path = ""
        cmd = f"which {program}"
        p1 = subprocess.run(
            cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
        )
        output = p1.stdout.decode("UTF-8")
        error = p1.stderr.decode("UTF-8")
        if not error:
            if output:
                path = output.rstrip("\n")
                return path
            else:
                msg = f" ERROR: Unable to find the PATH of {program}"
                logger.error(msg)
        else:
            msg = f" ERROR: Unable to find the PATH of {program}"
            logger.error(msg)

    def validate(self, dump_messages=True):
        """Simple check for docker images installed"""
        for resource in self._data:
            image = self._data[resource]["image"]
            cmd = f"{self.docker} image ls {image}"
            print(cmd)
            p1 = subprocess.run(
                cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
            )
            output = p1.stdout.decode("UTF-8")
            error = p1.stderr.decode("UTF-8")
            print(output, error, "heeey")
            if p1.returncode != 0:
                msg = f"docker image {image} for {resource} was not found"
                logger.error(msg)
                logger.error(error)
                sys.exit()

            if output:
                c_lines = 0
                for line in output.split("\n"):
                    c_lines += 1
                if c_lines != 3:
                    msg = f"docker image {image} for {resource} was not found"
                    logger.error(msg)
                    sys.exit()
                else:
                    if dump_messages:
                        msg = f"Found docker image {image} for {resource}"
                        logger.info(msg)
        
class IsoformIds(YamlConfig):
    def __init__(self, config_file) -> None:
        super().__init__(config_file)
        self.transcripts_ids = set()
        self.get_transcripts()
    
    def get_transcripts(self):
        with open(self.isoform_path, "r") as f:
            for line in f:
                if line.startswith("#"):
                    continue
                line = line.strip()
                self.transcripts_ids.add(line)

def load_config():
    main_dir = Path(__file__).resolve().parents[1]
    config_file = os.path.join(main_dir, "yaml_files", "config.yaml")
    ann_config = AnnotationConfig(config_file)
    ann_config.parse_data()
    ref_config = ReferenceConfig(config_file)
    ref_config.parse_data()
    docker_config = DockerConfig(config_file)
    docker_config.validate()
    isoforms_conf = IsoformIds(config_file)

    return (ann_config, ref_config, docker_config, isoforms_conf)


# ann_conf, ref_conf, docker_conf, isoforms_conf = load_config()

if __name__ == "__main__":

    # Example on how to create and acess attributes of Reference Genome and Annotations

    ann_config  = AnnotationConfig()
    ann_config.parse_data()
    ref_config = ReferenceConfig()
    ref_config.parse_data()

    print(ref_config.hg38.fasta)
    print(ann_config.yaml.full_path)

    hg38_data = ref_config.get_genome_data("hg38")


    # we can also load it with the load_config function:
    ann_config, ref_config = load_config()