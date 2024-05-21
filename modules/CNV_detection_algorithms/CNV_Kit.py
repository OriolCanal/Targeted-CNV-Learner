import os
import subprocess

from modules.log import logger

from modules.CNV_detection_algorithms.CNV_algorithm import CNV_Algorithm

class CNV_Kit():
    def __init__(self, docker_conf, reference_conf, Bed, force_run=False):
        super().__init__(docker_conf, reference_conf, Bed, force_run)
