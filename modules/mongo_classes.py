import sys
import os
import pymongo
import mongoengine as me
from collections import defaultdict
from datetime import datetime
import pyliftover
import shlex
import subprocess
from pyfaidx import Fasta
# from modules.utils import db_name, main_dir
from modules.log import logger
from bson import ObjectId


class MongoDBManager:
    def __init__(self, db_name):
        self.db_name = db_name

    def init_mongodb(self, config=None):
        """Initialize MongoDB connection and create the database if it does not exist."""
        
        if config == None:
            config = {
                "username": "",
                "password": "",
                "host": "localhost:27017",
                "port": "",
                "authentication_source": "",
                "authentication_mechanism": "",
                "ssl": "",
                "ssl:cert_reqs": "",
            }

        # Check if the database exists
        db_exists = self.check_database_exists()

        if not db_exists:
            # If the database does not exist, create it from a restore file
            self.create_db_from_restore_file("your_restore_filename_here")

        try:
            me.connect(self.db_name, host="localhost:27017", alias="default")
        except ConnectionError:
            raise ConnectionError("Unable to connect to MongoDB")
        else:
            logger.info("Connected to MongoDB!")

    def check_database_exists(self):
        """Check if the MongoDB database exists."""
        client = pymongo.MongoClient("localhost", 27017)
        database_names = client.list_database_names()

        if self.db_name in database_names:
            logger.info(f"Database '{self.db_name}' exists.")
            return True
        else:
            logger.info(f"Database '{self.db_name}' does not exist.")
            return False

    def create_db_from_restore_file(self, restore_mongo_path):
        """Create MongoDB database from a restore file."""
        initial_db_name = "db"
        cmd = f"mongorestore --archive={restore_mongo_path} --nsFrom='{initial_db_name}.*' --nsTo='{self.db_name}.*'"
        logger.info(f"Creating mongodb database:\n{cmd}")
        cmd_parts = shlex.split(cmd)
        result = subprocess.run(cmd_parts, capture_output=True, text=True)

        if result.returncode == 0:
            logger.info("MongoDB created successfully!")
        else:
            logger.info("Error creating MongoDB:")
            logger.info(result.stderr.decode("utf-8"))




class HPO_Term(me.DynamicDocument):
    hpo_id = me.StringField()
    hpo_name = me.StringField()
    disease_id = me.StringField()
    disease_name = me.StringField()
    frequency = me.StringField()
    timestamp = me.DateTimeField(
        required=True, default=datetime.now().strftime("%d/%m/%Y")
    )
    hpo_genes = me.ListField(me.ReferenceField("HPO_Gene"))
    meta = {
        "strict": False,
        "collection": "hpo_term",
        "indexes": [
            {"fields": ["hpo_id", "hpo_name", "disease_id", "disease_name"]},
            "hpo_term.hpo_id",  # Add this line to index hpo_id within the hpo_terms
        ],
    }


class HPO_Gene(me.DynamicDocument):
    ncbi_gene_id = me.StringField()
    gene_symbol = me.StringField()
    association_type = me.StringField()
    disease_id = me.StringField()
    source = me.StringField()
    hpo_terms = me.ListField(me.ReferenceField("HPO_Term"))
    meta = {
        "strict": False,
        "collection": "hpo_gene",
        "indexes": [
            {"fields": ["gene_symbol", "association_type", "disease_id"]},
        ],
    }
    timestamp = me.DateTimeField(
        required=True, default=datetime.now().strftime("%d/%m/%Y")
    )


# class VirtualPanel(me.DynamicDocument):
#     name = me.StringField()
#     parent_panel = me.StringField()
#     hpo_terms = me.ListField(me.StringField())
#     genes = me.ListField(me.StringField())
#     timestamp = me.DateTimeField(
#         required=True, default=datetime.now().strftime("%d/%m/%Y")
#     )



class Diseases(me.DynamicDocument):
    timestamp = me.DateTimeField(
        required=True, default=datetime.now().strftime("%d/%m/%Y")
    )


class Pathologies(me.DynamicDocument):
    timestamp = me.DateTimeField(
        required=True, default=datetime.now().strftime("%d/%m/%Y")
    )


class Annotation(me.DynamicDocument):
    md5 = me.StringField(required=True)
    timestamp = me.DateTimeField(
        required=True, default=datetime.now().strftime("%d/%m/%Y")
    )
    variant = me.ReferenceField("Variant")
    meta = {
        "strict": False,
        "collection": "annotations",
        "indexes": [{"fields": ["md5"]}],
    }

    def get_synonyms(self, Gene_Synonyms):
        if hasattr(self, "Gene_Synonyms"):
            return(None)
        synonyms = set()
        ann_symbol = self.CSQ["SYMBOL"]
        if ann_symbol in Gene_Synonyms.gene_synonyms:
            synonyms = Gene_Synonyms.get_gene_synonyms(ann_symbol)
        synonyms.add(ann_symbol)
        # print(f"adding synonyms: {synonyms}")
        self.Gene_Synonyms = list(synonyms)
        self.save()
        return(synonyms)

    def get_timestamp(self):
        # If annotation is a DBRef object, return None
        if isinstance(self, ObjectId):
            return None
        
        # If annotation is an object with a 'timestamp' attribute
        if hasattr(self, "timestamp"):
            return self.timestamp
    
    def get_transcript_version(annotation):
        hgvsc = annotation.CSQ["HGVSc"]
        if ":" in hgvsc:
            trans_id = hgvsc.split(":")[0]
            if "." in trans_id:
                transcript_version = trans_id.split(".")[1]
        
                return transcript_version
        return(None)
    
    def compare_anns_values(self, ann_obj, vep_field, variant):


        if hasattr(ann_obj, "CSQ") and hasattr(self, "CSQ"):
            prev_transcript_csq = ann_obj.CSQ
            actual_transcript_csq = self.CSQ
        
        if vep_field not in prev_transcript_csq:
            msg = f"Vep field: {vep_field} not present in annotation id: {ann_obj.id}"
            # logger.warning(msg)
            return(None)
        if vep_field not in actual_transcript_csq:
            msg = f"Vep field: {vep_field} not present in annotation id: {self.id}"
            # logger.warning(msg)
            return(None)    
        # ensuring that we are comparing properties of same transcript ID
        if prev_transcript_csq["Feature"] != actual_transcript_csq["Feature"]:
            raise(ValueError("Comparing annotations with different transcript IDs!"))
        
        transcript_id = actual_transcript_csq["Feature"]
        prev_field = prev_transcript_csq[vep_field]
        actual_field = actual_transcript_csq[vep_field]
        if prev_field != actual_field:
            # HGVS change format, so we have to take into account in order to not log the change of format.
            if vep_field == "HGVSp" or vep_field == "HGVSc":
                if ":" in actual_field:
                    transcript_id, actual_field = actual_field.split(":")

                if ":" in prev_field:
                    # print("hey",prev_field)
                    prev_field = prev_field.split(":")[1]

                if prev_field == actual_field:
                    return None

            logger.info(
                f"Actual Annotation id: {self.id}\nPrevious annotation id:\
                {ann_obj.id}\n\thas  changed {vep_field}: \n\tprevious {prev_field}\n\tactual: {actual_field}"
            )
            annotation_log = AnnotationsLog(
                transcript_id = transcript_id,
                variant = variant.id,
                field = vep_field,
                previous_value = prev_field,
                actual_value = actual_field,
                original_transcript_version =  ann_obj.get_transcript_version(),
                new_transcript_version = self.get_transcript_version(),
                timestamp=datetime.now()
            )
            annotation_log.save()

            # Add the annotation_log to the changes list of the variant
            variant.update(push__changes=annotation_log)

            # Save the updated variant
            variant.save()
            return(vep_field, prev_field, actual_field)
        return(None)
    
class ClinicalClass(me.DynamicDocument):
    user = me.StringField()
    acmg_guideline = me.StringField()
    population_data = me.StringField()
    prevalence_data = me.StringField()
    computation_data = me.StringField()
    segregation_data = me.StringField()
    functional_data = me.StringField()
    de_novo_data = me.StringField()
    allelic_data = me.StringField()
    other_data = me.StringField()
    clinical_class = me.StringField()
    latest_acmg = me.StringField()
    lab_id = me.ReferenceField("Sample")
    timestamp = me.DateTimeField(
        required=True, default=datetime.now().strftime("%d/%m/%Y, %H:%M:%S")
    )
    observations = me.StringField()


class LostExon(me.DynamicDocument):
    timestamp = me.DateTimeField(
        required=True, default=datetime.now().strftime("%d/%m/%Y")
    )


class Sample(me.DynamicDocument):
    """ """

    lab_id = me.StringField(required=True)
    external_ids = me.ListField(me.StringField())
    run_id = me.StringField(required=True)
    panel = me.StringField(required=True)
    panel_version = me.StringField(required=True)
    subpanel = me.StringField(required=True)
    roi_file = me.StringField(required=True)
    lost_exons = me.ListField(me.ReferenceField("LostExon"))
    diseases = me.ListField(me.StringField())
    variant_call = me.StringField()
    bam_file = me.StringField()
    bai_file = me.StringField()
    vcf_file = me.StringField()
    virtual_panels = me.ListField(me.ReferenceField("VirtualPanel"))

    timestamp = me.DateTimeField(
        required=True, default=datetime.now().strftime("%d/%m/%Y")
    )

class VirtualPanel(me.DynamicDocument):
    """ """
    name = me.StringField(required=True)
    parent_panel = me.StringField(required=True)
    # parent_panel_version = me.StringField(required=True)    
    genes = me.ListField(me.StringField(required=True))
    hpo_terms = me.ListField(me.StringField())
    timestamp = me.DateTimeField(
        required=True, default=datetime.now().strftime("%d/%m/%Y")
    )


class Run(me.DynamicDocument):
    """ """

    run_id = me.StringField(required=True)
    panel = me.StringField(required=True)
    panel_version = me.StringField(required=True)
    subpanel = me.StringField(required=True)
    roi_file = me.StringField()
    samples = me.ListField(me.StringField())
    timestamp = me.DateTimeField(
        required=True, default=datetime.now().strftime("%d/%m/%Y")
    )
    meta = {
        "strict": False,
        "collection": "runs",
        "indexes": [{"fields": ["run_id", "panel", "panel_version"]}],
    }


class Variant(me.DynamicDocument):
    """ """

    # Required fields
    chromosome = me.StringField(required=True)
    pos = me.IntField(required=True)
    ref = me.StringField(required=True)
    alt = me.StringField(required=True)
    genome_version = me.StringField(required=True)
    rsid = me.StringField()
    invalid_format = me.BooleanField(default=False)
    invalid_reason = me.StringField()
    timestamp = me.DateTimeField(
        required=True, default=datetime.now().strftime("%d/%m/%Y")
    )
    calls = me.ListField(me.ReferenceField("Call"))
    changes = me.ListField(me.ReferenceField("AnnotationsLog"))
    annotations = me.ListField(me.ReferenceField(Annotation))
    clinical_class = me.ListField(me.ReferenceField(ClinicalClass))
    meta = {
        "strict": False,
        "collection": "variants",
        "indexes": [{"fields": ["chromosome", "pos", "ref", "alt"]}],
    }

    def __repr__(self):
        return f"Variant: {self.chromosome}:{str(self.pos)}{str(self.ref)}>{self.alt}"

    def __str__(self):
        return f"Variant: {self.chromosome}:{str(self.pos)}{str(self.ref)}>{self.alt}"

    def set_genome_v_grch37_to_hg19(self):
        if self.genome_version == "GRCh37/hg19":
            self.genome_version = "hg19"
            self.save()

    def is_sv(self):
        possible_variants = ["<DEL>", "<DUP>", "<BREAK>", "<BND>", "<INS>", "<DUP:TANDEM>"]
        if self.alt in possible_variants:
            return (True)
        
        # opening of breakpoint
        elif "[" in self.alt and "chr" in self.alt:
            return True
        
        # clousure of breakpoint junction:
        elif "]" in self.alt and "chr" in self.alt:
            return True

        return(False)

    def is_spliceai_indel(self):
        # to obtain if it's sv

        if self.is_sv():
            return(False)
        # there can be multiple alternative alleles
        alts = []
        if "," in self.alt:
            alts = self.alt.split(",")
        else:
            alts.append(self.alt)
        for alt in alts:
            # insertions of more than 1 base
            if len(alt) - len(self.ref) > 1:
                return(True)

            # deletions of more than 4 bases
            if len(self.ref) - len(alt) > 4:
                return(True)

        return(False)
    
    def check_variant(self, fasta_file):
        """
        Check if the variants are correctly formatted. If not, it adds a new feature to mongodb
        called invalid_format (boolean) that is True if the variant is incorrectly formatted and false
        if it is correctly formatted
        Params:
            fasta_file (str): reference fasta file to extract the reference nucleotide in case the variant doesn't have one.
        
        Returns:
            self.invalid_format (bool): variable that indicates if variant is correctly (False) or incorrectly (True) formatted
            Bool: 0 if incorrectly formatted, 1 if correclty formatted
        """
        if self.ref.upper() == self.alt.upper():
            self.invalid_format = True
            self.invalid_reason = "Same_Ref_and_Alt"
            self.save()
            return 0

        # if + in ref or alt (invalid format, variant called by sampileup)
        if any("+" in var for var in [self.ref, self.alt]):
            logger.info(f"Variant with invalid format: {self}")
            self.invalid_format = True
            self.invalid_reason = "+_in_Ref_or_Alt"
            self.save()
            return 0

        # if "." in self.alt:
        #     logger.info(f"Variant with invalid format: {self}")
        #     self.invalid_format = True
        #     self.save()
        #     return 0

        # opening of breakpoint junction:
        if "[" in self.alt and "chr" in self.alt:
            # correct format
            self.invalid_format = False
            return 1
        # clousure of breakpoint junction:
        if "]" in self.alt and "chr" in self.alt:
            # correct format
            self.invalid_format = False
            return 1

        # if digit in ref or alt
        if any(char.isdigit() for char in self.ref) or any(char.isdigit() for char in self.alt):
            logger.info(f"Variant with invalid format: {self}")
            self.invalid_format = True
            self.invalid_reason = "Digit_in_Ref_or_Alt"
            self.save()
            return 0

        if self.ref == "" or self.ref == " ":
            logger.info(f"Variant with no reference, extracting ref from reference fasta file: {fasta_file}")
            new_ref = self.extract_fasta_nucleotide_from_position(fasta_file)
            if new_ref.upper() != self.alt.upper():
                logger.info(f"Setting new reference {new_ref} in variant: {self}")
                self.ref = new_ref.upper()

            else:
                self.invalid_format = True
                self.invalid_reason = "Same_Ref_and_Alt"
                self.save()
        
        valid_nucleotides = ["A", "C", "G", "T", "N", ",", "*", "."]
        possible_variants = ["<DEL>", "<DUP>", "<BREAK>", "<BND>", "<INS>", "<DUP:TANDEM>"]
        if self.alt not in possible_variants:
            if not all(letter.upper() in valid_nucleotides for letter in self.ref):
                logger.info(f"Variant with no ACTGN, in ref. invalid format: {self}")
                self.invalid_format = True
                self.invalid_reason = "Not_valid_nucl_in_Ref"
                self.save()
                return 0
            if not all(letter.upper() in valid_nucleotides for letter in self.alt):
                logger.info(f"Variant with no ACTGN, in alt. invalid format: {self}")
                self.invalid_format = True
                self.invalid_reason = "Not_valid_nucl_in_Alt"
                self.save()
                return 0
        else:
            # valid_nucleotides.remove(".")
            if not all(letter.upper() in valid_nucleotides for letter in self.ref):
                logger.info(f"Structural variant with no ACTGN, in ref. invalid format: {self}")
                self.invalid_format = True
                self.invalid_reason = "SV_with_wrong_ref"
                self.save()
                return 0
        
        self.invalid_format = False

        self.save()


    def extract_fasta_nucleotide_from_position(self, fasta_file):
        
        
        fasta = Fasta(fasta_file)

        
        try:
            nucleotide = fasta[self.chromosome][self.pos-1].seq
            fasta.close()
            return nucleotide
        except KeyError:
            logger.error(f"Chromosome {self.chromosome} not found in the FASTA file: {fasta_file}")
        except IndexError:
            logger.error(f"Position {self.pos} out of range for chromosome {self.chromosome}")


        


    def get_most_recent_call(self):
        """Get the most recent Call associated with the Variant based on analysis date."""
        # calls = list()
        # for call in self.calls:
        #     call_obj = call.fetch()
        #     if call_obj:
        #         calls.append(call_obj)
        valid_calls = [call for call in self.calls if hasattr(call, 'analysis_date')]

        # Sort the calls based on analysis_date in descending order
        sorted_variant_calls = sorted(valid_calls, key=lambda x: x.analysis_date, reverse=True)

        # Retrieve the most recent (last) call
        if sorted_variant_calls:
            most_recent_call = sorted_variant_calls[0]
            # print(most_recent_call)
            return most_recent_call
        else:
            return None
    
    def get_ann(self) -> list:
        variant_ann = [ann for ann in self.annotations]

        if variant_ann:
            return variant_ann

    def get_position_and_alleles(self):
        return (f"{self.chromosome}:{self.pos}{self.ref}>{self.alt}")
    
    def get_position_and_ref(self):
        return (f"{self.chromosome}:{self.pos}{self.ref}")

    def register(self) -> bool:
        """ """
        curr_time = datetime.now().strftime("%d/%m/%Y")

        # if doesnt exist
        if not Variant.objects(
            chromosome=self.chromosome,
            pos=self.pos,
            ref=self.ref,
            alt=self.alt,
            timestamp=curr_time,
        ):
            self.save(
                chromosome=self.chromosome,
                pos=self.pos,
                ref=self.ref,
                alt=self.alt,
                timestamp=curr_time,
            )
            return True
        else:
            return False

    def get_lo_obj(self, from_genome, to_genome):
        """
        Creates and returns a liftover object from `pyliftover`
        library based on the provided `from_genome` and `to_genome` parameters.

        Parameters:
            from_genome (str): The source genome version, e.g., "hg19".
            to_genome (str): The target genome version, e.g., "hg38".

        Returns:
            list: A list containing the liftover object, source genome version,
                  and target genome version.

        Raises:
            ValueError: If the provided genome versions are not supported.
        """
        if from_genome == "hg19" and to_genome == "hg38":
            return ([pyliftover.LiftOver("hg19", "hg38"), from_genome, to_genome])
        elif from_genome == "hg38" and to_genome == "hg19":
            return ([pyliftover.LiftOver("hg38", "hg19"), from_genome, to_genome])
        else:
            msg = f"genome version {self.genome_version} not able to do liftover"
            raise ValueError(msg)
        
    def do_liftover(self, lo_obj):
        """
        Perform liftover for the variant using a pre-existing liftover object.

        This method takes a liftover object (`lo_obj`) from the self.get_lo_obj()
        library and performs liftover for the current variant.
        self.get_lo_obj() is not called inside this function as it's time consuming
        Parameters:
            lo_obj (list): A list containing the liftover object, source genome version,
                        and target genome version, obtained from `get_lo_obj` method.

        Returns:
            int or None: The lifted-over position in the target genome (1-based), or None
                        if the variant was not mapped.
        """
        lo, from_genome, to_genome = lo_obj

        # vcf positions are 1-based and pyliftover needs 0-based positions 
        zero_based_pos = int(self.pos) - 1 

        results = lo.convert_coordinate(self.chromosome, zero_based_pos)

        
        if len(results) == 0:
            logger.warning(f"Variant not mapped: {self} from genome {from_genome} to genome {to_genome}")
            one_based_mapped_pos = None
        else:
            mapped_chrom, mapped_pos, _, _ = results[0]
        one_based_mapped_pos = int(mapped_pos) + 1

        if from_genome == "hg19":
            self.hg19_pos = self.pos
            self.pos = one_based_mapped_pos
            logger.info(f"Setting liftovered coords: \n\thg19: {self.hg19_pos}\n\thg38: {self.pos}")
        elif from_genome == "hg38":
            self.hg19_pos = one_based_mapped_pos
            logger.info(f"Setting liftovered coords: \n\thg19: {self.hg19_pos}\n\thg38: {self.pos}")

        self.save()
        return (one_based_mapped_pos)



class AnnotationsLog(me.DynamicDocument):
    """ """

    # Required fields
    transcript_id = me.StringField(required=True)
    field = me.StringField(required=True)
    variant = me.ReferenceField("Variant")
    previous_value = me.StringField()
    actual_value = me.StringField()
    original_transcript_version = me.IntField()
    new_transcript_version = me.IntField()
    timestamp = me.DateTimeField(
        required=True, default=datetime.now().strftime("%d/%m/%Y")
    )








class Call(me.DynamicDocument):
    """ """

    chromosome = me.StringField(required=True)
    pos = me.IntField(required=True)
    ref = me.StringField(required=True)
    alt = me.StringField(required=True)
    lab_id = me.StringField(required=True)
    external_id = me.ListField(me.StringField())
    run_id = me.StringField(required=True)
    panel = me.StringField(required=True)
    panel_version = me.StringField(required=True)
    subpanel = me.StringField(required=True)
    roi_file = me.StringField(required=True)
    analysis_date = me.DateTimeField(
        required=True, default=datetime.now().strftime("%d/%m/%Y")
    )
    timestamp = me.DateTimeField(
        required=True, default=datetime.now().strftime("%d/%m/%Y")
    )
    hg19_pos = me.IntField()
    hg19_hgvsg = me.StringField()
    variant_type = me.StringField()
    hgvsg = me.StringField()
    hgvsp = me.StringField()
    hgvsc = me.StringField()
    gene = me.StringField()
    isoform = me.StringField()
    rsid = me.StringField()
    clinical_class = me.ListField(me.ReferenceField("ClinicalClass"))

    algorithms = me.ListField(me.StringField())
    variant_quality = me.StringField()
    variant_filter = me.StringField()
    genotype = me.StringField()
    read_support = me.StringField()
    depth = me.StringField()
    vaf = me.StringField()
    genotype_quality = me.StringField()
    phred_likelihood = me.StringField()
    copy_number = me.StringField()
    variant = me.ReferenceField(Variant)
    annotation = me.ReferenceField(Annotation)
    observation_variants = me.StringField()
    date_comment = me.StringField()

    meta = {
        "strict": False,
        "collection": "calls",
        "indexes": [
            {"fields": ["lab_id", "run_id", "chromosome", "pos", "ref", "alt"]}
        ],
    }

    def __str__(self):
        return (f"Call: {self.chromosome}:{self.pos}{self.ref}>{self.alt}")

    def __repr__(self):
        return (f"Call: {self.chromosome}:{self.pos}{self.ref}>{self.alt}")

    def get_lo_obj(self, from_genome, to_genome):
        """
        Creates and returns a liftover object from `pyliftover`
        library based on the provided `from_genome` and `to_genome` parameters.

        Parameters:
            from_genome (str): The source genome version, e.g., "hg19".
            to_genome (str): The target genome version, e.g., "hg38".

        Returns:
            list: A list containing the liftover object, source genome version,
                  and target genome version.

        Raises:
            ValueError: If the provided genome versions are not supported.
        """
        if hasattr(self, "variant"):
            variant = self.variant

        if self.from_genome == "hg19" and to_genome == "hg38":
            return ([pyliftover.LiftOver("hg19", "hg38"), from_genome, to_genome])
        elif from_genome == "hg38" and to_genome == "hg19":
            return ([pyliftover.LiftOver("hg38", "hg19"), from_genome, to_genome])
        else:
            msg = f"genome version {self.genome_version} not able to do liftover"
            raise ValueError(msg)
        
    def do_liftover(self, lo_obj):
        """
        Perform liftover for the variant using a pre-existing liftover object.

        This method takes a liftover object (`lo_obj`) from the self.get_lo_obj()
        library and performs liftover for the current variant.
        self.get_lo_obj() is not called inside this function as it's time consuming
        Parameters:
            lo_obj (list): A list containing the liftover object, source genome version,
                        and target genome version, obtained from `get_lo_obj` method.

        Returns:
            int or None: The lifted-over position in the target genome (1-based), or None
                        if the variant was not mapped.
        """
        # liftover already performed
        if self.pos and self.hg19_pos:
            return 0
        lo, from_genome, to_genome = lo_obj

        # vcf positions are 1-based and pyliftover needs 0-based positions 
        zero_based_pos = int(self.pos) - 1 

        results = lo.convert_coordinate(self.chromosome, zero_based_pos)

        
        if len(results) == 0:
            logger.warning(f"Variant not mapped: {self} from genome {from_genome} to genome {to_genome}")
            one_based_mapped_pos = None
        else:
            mapped_chrom, mapped_pos, _, _ = results[0]
            one_based_mapped_pos = int(mapped_pos) + 1

        if from_genome == "hg19":
            self.hg19_pos = self.pos
            self.pos = one_based_mapped_pos
            hgvsg = self.hgvsg
            if not hgvsg == ".":
                hg38_hgvsg = f"{self.chromosome}:g.{self.pos}{self.ref}>{self.alt}"
                self.hgvsg = hg38_hgvsg
                self.hg19_hgvsg = hgvsg



            logger.info(f"Setting Call \n\nliftovered coords: \n\thg19: {self.chromosome}:{self.hg19_pos}\n\thg38: {self.pos}\n\t id: {self.id}")
        elif from_genome == "hg38":
            self.hg19_pos = one_based_mapped_pos
            hgvsg = self.hgvsg
            if not hgvsg == ".":
                hg19_hgvsg = f"{self.chromosome}:g.{self.hg19_pos}{self.ref}>{self.alt}"
                self.hg19_hgvsg = hg19_hgvsg
            logger.info(f"Setting liftovered coords: \n\thg19: {self.hg19_pos}\n\thg38: {self.pos}")

        self.save()
        return (one_based_mapped_pos)

