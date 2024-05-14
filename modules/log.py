import logging
from pathlib import Path
import os


# Define the path to the log file
log_file = os.path.join(str(Path(__file__).resolve().parents[1]), "my_package.log")

# Create a logger for the package
logger = logging.getLogger('my_package')
logger.setLevel(logging.DEBUG)  # Set the logger level to DEBUG

# Create a file handler for the log file
file_handler = logging.FileHandler(log_file)
file_handler.setLevel(logging.DEBUG)  # Set the file handler level to DEBUG

# Create a stream handler for printing to console
console_handler = logging.StreamHandler()
console_handler.setLevel(logging.INFO)  # Set the console handler level to INFO

# Create a formatter
formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')

# Set the formatter for both handlers
file_handler.setFormatter(formatter)
console_handler.setFormatter(formatter)

# Add the handlers to the logger
logger.addHandler(file_handler)
logger.addHandler(console_handler)