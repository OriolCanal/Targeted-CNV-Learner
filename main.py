import pyliftover
import os
import sys
import pandas as pd
# from modules.create_mongo import MongoDBManager
from modules.mongo_classes import Variant, Call, MongoDBManager, Annotation
from modules.mongo_classes import Run as Mongo_Run
from modules.mongo_classes import Sample as Mongo_Sample
from modules.bed import Bed

from modules.CNV_detection_algorithms.gatk_gCNV import Gatk_gCNV, Case_Gatk_gCNV, Cohort_Gatk_gCNV
from modules.params import load_config
from modules.log import logger
from modules.picard import Picard, Metrics_Df
from modules.mosdepth import Mosdepth, Mosdepth_df
from modules.clustering import Clustering_Mosdepth
from modules.run_class import Analysis_Run, Analysis_Sample 
from modules.utils import get_timestamp

# # initializing mongodb
mongo_manager = MongoDBManager("master_db")

logger.info(f"loading configuration files.")
ann_conf, ref_conf, docker_conf, isoforms_conf = load_config()

mongo_db_created = False # if mongodb has been created in this script the variable is True
if not mongo_manager.check_database_exists():
    restore_file = ref_conf.mongo_restore_file
    mongo_manager.create_db_from_restore_file(restore_file)
    mongo_db_created = True

mongo_manager.init_mongodb()

Picard_Obj = Picard(docker_conf, ref_conf)
Picard_Obj.create_fasta_dict()

Mosdepth_Obj = Mosdepth(docker_conf)

Picard_Metrics_Df = Metrics_Df()
Mosdepth_Metrics_Df = Mosdepth_df(ref_conf)

Bed_obj = Bed(ref_conf.bed)
# -----------ITERATE OVER RUNS----------------
all_runs = Mongo_Run.objects
time_ordered_runs = sorted(all_runs, key=get_timestamp, reverse=True)

run_ids = list()
for run in time_ordered_runs:
    if run.panel == "SUDD_147":
        run_ids.append(run.run_id)
        logger.info(f"{run.run_id} will be analysed")
    if len(run_ids) == 7:
        break

logger.info(f"{len(run_ids)} runs will be downloaded")

analysis_samples = dict()

for run in run_ids:
    Run_obj = Analysis_Run(run)
    Run_obj.get_samples(Mongo_Run, Mongo_Sample)

    for sample_147 in Run_obj.samples_147:

        Bam = sample_147.get_bam_bai_from_compendi(ref_conf)

        # check if files have been downloaded correctly
        if sample_147.bam:
            analysis_samples[sample_147.sample_id] = sample_147
        else:
            continue
        # Picard
        Picard_Obj.run_collectHsMetrics(Bam, Bed_obj)
        info_line, value_line = Picard_Obj.get_picard_metrics(Bam)
        if not Picard_Metrics_Df.has_df_header():
            Picard_Metrics_Df.add_metrics_header(info_line)
        Picard_Metrics_Df.add_metrics_line(value_line)

        # Mosdepth
        Mosdepth_Obj.run_mosdepth(Bam.path,force=False)
        exon_coverage, sample_mean_coverage = Mosdepth_Obj.parse_mosdepth_regions_bed(Run_obj.run_id)
        Mosdepth_Metrics_Df.add_normalized_mean_coverage_dict(exon_coverage, sample_mean_coverage)


mosdepth_excel_path = os.path.join(ref_conf.main_dir, "mosdepth_excel.xlsx")
exons_mean_coverage_df = Mosdepth_Metrics_Df.get_df_from_normalized_exons_coverage()
logger.info(f"Creating mosdepth coverage dataframe in {mosdepth_excel_path}")
exons_mean_coverage_df.to_excel(mosdepth_excel_path, index=False)

heatmap_path = os.path.join(ref_conf.main_dir, "Mosdepth_heatmap.png")
Mosdepth_Metrics_Df.apply_pca(normalized=True)
# closes_points = Mosdepth_Metrics_Df.get_n_pca_closest_samples("RB36157", n=5)
Mosdepth_Metrics_Df.apply_umap(normalized=True)
Cluster_Pca = Clustering_Mosdepth(Mosdepth_Metrics_Df, "PCA", ref_conf)
# Cluster_Pca.apply_bayesian_gaussian_mixture()
# Cluster_Pca.do_hierarchical_clustering()
Cluster_Pca.apply_hdbscan(min_samples=15, min_cluster_size=5)

Cluster_Umap = Clustering_Mosdepth(Mosdepth_Metrics_Df, "UMAP", ref_conf)
# Cluster_Umap.apply_bayesian_gaussian_mixture()
# Cluster_Umap.do_hierarchical_clustering()
Cluster_Umap.apply_hdbscan(min_samples=15, min_cluster_size=3)

cluster_samples = Cluster_Pca.get_hdbscan_cluster_samples(analysis_samples)

# if force_run = True, if output file of gatk exists, the command will run and the file will be overwritten
force_run = False
# # Run GATK gCNV
Picard_Obj.run_bed_to_interval_list(Bed_obj)
Gatk_obj = Gatk_gCNV(docker_conf, ref_conf, Bed_obj, force_run)
Gatk_obj.run_preprocess_intervals()
Gatk_obj.run_index_feature_file()
Gatk_obj.run_annotate_intervals()

# creating hdf5 file for each sample
for sample in analysis_samples.values():

    Bam = sample.bam
    Gatk_obj.run_collect_read_counts(Bam)

# creating a model for each sample
for cluster in cluster_samples.keys():
    if cluster == -1:
        continue
    cohort_samples = cluster_samples[cluster]
    print(cohort_samples)
    bams = list()
    for sample in cohort_samples:
        Gatk_obj.run_filter_intervals(sample, cohort_samples)
        Gatk_obj.run_interval_list_tools()
        Gatk_obj.run_filter_intervals(sample, cohort_samples)
        Gatk_Cohort = Cohort_Gatk_gCNV(sample, docker_conf, ref_conf, Bed_obj, cohort_samples, force_run)
        Sample_Case_Gatk_gCNV = Case_Gatk_gCNV(sample, docker_conf, ref_conf, Bed_obj, force_run)
        Gatk_Cohort.run_determine_germline_contig_ploidy(sample)
        Sample_Case_Gatk_gCNV.run_determine_germline_contig_ploidy(Gatk_Cohort)
        Gatk_Cohort.run_germline_cnv_caller(sample)
        Sample_Case_Gatk_gCNV.run_germline_cnv_caller(Gatk_Cohort)
        Sample_Case_Gatk_gCNV.run_postprocess_germline_calls(Gatk_Cohort)
        bams.append(sample.bam.path)
    print(f"bams of cluster {cluster}:\n {bams}")

# should create a 
# Gatk_obj.run_determine_germline_contig_ploidy(Bed_obj)

# --------------FINISH ITERATING OVER RUNS---------


        
# for sample in samples:
#     bam_path, bai_path = get_bam_bai_from_compendi(sample, "RUN20240208-CGC64001", ref_conf)

# annotations = Annotation.objects
# variants_cnvs = list()
# cnvs_panel_147 = set()
# for i, ann in enumerate(annotations): # [87999:]
#     if i % 1000 == 0:
#         logger.info(f"{i} annotations analysed, {len(variants_cnvs)} found in panel 147")
#     if not hasattr(ann, "TYPE"):
#         continue
#     if ann.TYPE != "CNV":
#         continue
#     # print("CNV found")
#     cnv_variant = ann.variant
#     # print(cnv_variant.id)

#     variant_calls = cnv_variant.calls
    
#     variant_calls_ids = set()
#     cnv_info = dict()
#     for variant_call in variant_calls:
#         try:
#             call = Call.objects.get(id=variant_call.id)
#             # print(call.id)
#             if hasattr(call, "panel"):
#                 # print(call.panel)
#                 if call.panel == "147" or call.panel == "SUDD_147":
            
#                     cnv_info = {
#                         "variant_id": cnv_variant.id,
#                         "chr": cnv_variant.chr,
#                         "pos": cnv_variant.pos,
#                         "ref": cnv_variant.ref,
#                         "alt": cnv_variant.alt,
#                         "genome_version": cnv_variant.genome_version,
#                         "run_id": variant_call.run_id,
#                         "panel": variant_call.panel,
#                         "lab_id": variant_call.lab_id,
#                         "subpanel": variant_call.subpanel,
#                         "panel_version": variant_call.panel_version,
#                     }
#                     # print(cnv_info)
#                     # print("CNV 147 panel")
                    
#                     variants_cnvs.append(cnv_info)
#                     # print(f"variants cnv: {variants_cnvs}")
#         except:
#             pass
#             # logger.info(f"Call ID {variant_call.id} does not exists, it has been removed due to false positive interpretation.")

# # Initialize an empty set to store hashable representations of dictionaries
# seen = set()

# # Initialize an empty list to store unique dictionaries
# unique_variants_cnvs = []

# # Iterate over the list of dictionaries
# for variant_cnv in variants_cnvs:
#     # Convert the dictionary to a tuple of its items for hashing
#     variant_cnv_tuple = tuple(variant_cnv.items())
    
#     # Check if the tuple representation has been seen before
#     if variant_cnv_tuple not in seen:
#         # If not seen, add it to the set and the list of unique dictionaries
#         seen.add(variant_cnv_tuple)
#         unique_variants_cnvs.append(variant_cnv)
        
# df = pd.DataFrame(unique_variants_cnvs)
# df.to_excel("/home/ocanal/Desktop/master_thesis_CNV/cnv_excel.xlsx", index=False)
# print(len(unique_variants_cnvs))

# df = pd.read_excel("/home/ocanal/Desktop/master_thesis_CNV/cnv_excel.xlsx")

# # Function to filter rows within each group
# def filter_within_group(group):
#     # Sort the group by 'pos'
#     sorted_group = group.sort_values(by='pos')
    
#     # Compute the difference between consecutive 'pos' values
#     diff = sorted_group['pos'].diff().fillna(0)
    
#     # Keep rows where the difference is greater than 200 or it's the first row
#     return sorted_group[diff >= 200]
#         # print("CNV found with call")
#         # count = Call.objects(id=variant_call.id).count()
#         # if count > 0:
#         #     variants_cnvs.add(cnv_variant)



# # Apply the filter function to each group
# filtered_df = df.groupby(['lab_id', 'chr']).apply(filter_within_group).reset_index(drop=True)
# for cnv_variant in variants_cnvs:
#     cnv_info = dict()
#     cnv_info = {
#         "variant_id": cnv_variant.id,
#         "pos": cnv_variant.pos,
#         "ref": cnv_variant.ref,
#         "alt": cnv_variant.alt,
#         "genome_version": cnv_variant.genome_version
#     }







    
    
    # if not hasattr(ann, "CSQ"):
    #     continue
    # # print("has CSQ")
    # if not "type" in ann.CSQ:
    #     continue
    # print("has type")
    # variant_type = ann.CSQ["Type"]
    # if variant_type != "CNV":
    #     continue
    # print("type is CNV")
    # cnv_variant = ann.variant

    # variant_calls = cnv_variant.calls

    # if any(call for call in variant_calls if call.panel == "147"):
    #     cnvs_panel_147.add(cnv_variant)

# print(len(cnvs_panel_147))
        
















# remove_annotations_without_csq(Annotation=Annotation, Variant=Variant)

# # get gene synonyms
# Gene_Synonyms = GeneSynonyms(ann_conf)

# # parsejar-ho amb batches
# all_variants = Variant.objects
# all_anns = Annotation.objects
# all_calls = Call.objects


# logger.info(f"Getting Gene Synonyms from {Gene_Synonyms.gene_synonyms_path} and including them into annotations")
# for ann in all_anns:
#     ann.get_synonyms(Gene_Synonyms)

# hg19_lo_obj = [pyliftover.LiftOver("hg19", "hg38"), "hg19", "hg38"]
# hg38_lo_obj = [pyliftover.LiftOver("hg38", "hg19"), "hg38", "hg19"]

# # liftover over calls based on it's variant associated genome version
# for call in all_calls:
#     if hasattr(call, "variant"):
#         if call.variant is None:
#             call.do_liftover(hg19_lo_obj)
#             continue
#         variant = call.variant
#         if variant.genome_version == "hg19":
#             call.do_liftover(hg19_lo_obj)
#         elif variant.genome_version == "hg38":
#             call.do_liftover(hg38_lo_obj)

# vcf_creator = VcfFile()

# for variant in all_variants:

#     # performing liftover to all variants
#     variant.set_genome_v_grch37_to_hg19()
#     if variant.genome_version == "hg19":
#         if not hasattr(variant, "hg19_pos"):
#             variant.do_liftover(hg19_lo_obj)
#     elif variant.genome_version == "hg38":
#         if not hasattr(variant, "hg19_pos"):
#             variant.do_liftover(hg38_lo_obj)

#     # checking if variant format is correct
#     variant.check_variant(ref_conf.hg38.fasta_path)

#     # if variant contains invalid format, skip it
#     if variant.invalid_format:
#         logger.info(f"{variant} with invalid format, it won't be added to the vcf file")
#         continue

#     vcf_creator.add_variant(variant)

# vcf_creator.create_vcf_file()

# logger.info(f"correcting reference hg38 allele from fasta file")
# snp_vcf, indel_vcf, sv_vcf = vcf_creator.get_new_vcf_with_corrected_ref(ref_conf)


# vep_obj = Vep(ref_conf, ann_conf, docker_conf, "hg38")

# # run vep over all_snp.vcf
# vep_vcf = vep_obj.run_vep(snp_vcf)
# # run vep over sv
# sv_vep_vcf =  vep_obj.run_vep(sv_vcf, structural_variants=True)

# # vep_vcf = "/home/ocanal/Desktop/reannotate_annotations_db/vcf_files/vep_no_sv_all_variants.vcf"
# # sv_vep_vcf = "/home/ocanal/Desktop/reannotate_annotations_db/vcf_files/vep_sv_all_variants.vcf"
# snp_vep_obj = Vep_vcf(vep_vcf)
# sv_vep_obj = Vep_vcf(sv_vep_vcf)

# logger.info(f"Updating annotations obtained from VEP in MongoDB:")
# # update annotations in mongodb from the vep vcf into mongodb
# snp_vep_obj.update_ann_mongo(Gene_Synonyms, Variant=Variant, isoform_conf=isoforms_conf)
# sv_vep_obj.update_ann_mongo(Gene_Synonyms, Variant=Variant, isoform_conf=isoforms_conf, sv=True)
# # indel_vcf = "/home/ocanal/Desktop/reannotate_annotations_db/vcf_files/trial_spliceai_indels.vcf"

# # Run SpliceAI over indels vcf
# splice_vcf = SpliceaiVcf(indel_vcf, docker_conf, ann_conf, ref_conf)
# splice_vcf.run_spliceai()
# for splice_variant in splice_vcf.get_spliceai_variants():
#     # update spliceai fields into mongodb
#     splice_variant.change_spliceai_field_in_annotation(Gene_Synonyms)

# # compare vep_fields annotations for the same transcript and if there are any change it is added to the annotations_log collection
# log_annotations_changes(all_variants, vep_fields=["EXON", "INTRON", "HGVSc", "HGVSp", "HGVSg"])


        # Actually we don't have all the transcript versions as most of annotations are not referenced by a transcript version
# # Check if any of the results are not None
# if any(result is not None for result in [exon_result, intron_result, hgvsc_result, hgvsp_result, hgvsg_result]):
#     print("At least one result is not None")
# else:
#     if not prev_ann_v == actual_ann_v:
#         annotation_log = AnnotationsLog(
#             variant = variant.id,
#             transcript_id = transcript_id,
#             field = "Transcript_version",
#             previous_value = prev_ann_v,
#             actual_value = actual_ann_v,
#             original_transcript_version =  prev_ann_v,
#             new_transcript_version = actual_ann_v,
#             timestamp=datetime.now()
#         )
#         logger.info(f"changed transcript version in {transcript_id} from v.{prev_ann_v} to v.{actual_ann_v} \
#             but any change in exon, intron, HGVSc, HGVSp, HGVSg")
#         annotation_log.save()
#         variant.update(push__changes=annotation_log)

#         # Save the updated variant
#         variant.save()

# compare_anns_values(ann, prev_ann, "CANONICAL")
# compare_anns_values(ann, prev_ann, "HGVSg")
# compare_anns_values(ann, prev_ann, "GNOMAD_AC")
# compare_anns_values(ann, prev_ann, "ClinVar_CLNSIG")
# compare_aconns_values(ann, prev_ann, "Consequence")
# compare_anns_values(ann, prev_ann, "MANE_SELECT")
#   compare_anns_values(ann, prev_ann, "RefSeq")
# compare_anns_values(ann, prev_ann, "gnomADe_AF")
# compare_anns_values(ann, prev_ann, "MAX_AF")


# example of annotation id that should be removed from annotations and variant collections
# # { "annotations": ObjectId("6307da5bb1c4467a24133bf6") }
# { "timestamp": { "$gt": ISODate("2024-04-01T00:00:00.000Z") } }