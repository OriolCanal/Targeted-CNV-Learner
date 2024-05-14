# CNV Detection on Targeted Sequencing
Descripció:
The master's thesis will be conducted at the Unitat de Diagnòstic Molecular i Medicina Personalitzada (UDMMP) at Hospital Josep Trueta. Our team specializes in genetic diagnostics for various diseases with a significant genetic component, placing particular emphasis on cardiovascular diseases as we are recognized as a leading center for cardiovascular disease diagnosis in Catalunya.

 

The focus of the master's thesis will be on developing a model to enhance the identification of copy number variants (CNVs). CNVs are commonly defined as a subtype of Structural Variant (SV) larger than 50 bp that involves a net change in the copy number of specific DNA regions. CNVs include duplications (where an additional copy of the DNA segment is present), deletions (where a segment is missing), insertions (where a segment is inserted from another part of the genome), and other forms of complex rearrangements.

 

CNVs have been reported in various genes associated with inherited cardiomyopathies, playing a pathogenic role in approximately 2% of patients [1]. Therefore, accurate CNV detection is a crucial aspect of genetic analysis pipelines.

 

Targeted sequencing applications, including exome and gene panels, have been widely established for clinical genetic diagnosis due to its reduced costs of sequencing and data storage compared to whole genome sequencing. CNV detection on targeted sequencing data imposes additional difficulties compared with genome sequencing for two main reasons. First, the sparseness of DNA capturing approaches hampers the use of breakpoint signals. And second, a highly variable depth of coverage caused by sampling, GC (guanine cytosine)-content and bait capturing efficiency.

 

While many algorithms have been developed for CNV detection based on targeted sequencing data in genetic diagnostics [1], most excel in reliably detecting large CNVs (on the order of megabases). However, they often exhibit poor performance when dealing with small CNVs that affect only one or a few exons, which are typically implicated in various genetic diseases.

 

Recent studies have shown that the integration of multiple CNV detection algorithms can significantly improve the identification of high-confidence CNVs [2]. Therefore, my proposal is to develop an artificial intelligence model specifically trained on samples from our lab. This model will leverage genetic data commonly used for CNV detection, including predictions and scores from various software (tools such as DECON[3], GATK-gCNV [4]), as well as parameters that can influence the quality of the call as quality parameters of the sample and mappability and GC content of the CNV region.

 

The primary goal of this model is to outperform individual software predictions and provide more accurate CNV predictions. By doing so, we aim to integrate this model into clinical practice to reduce the need for manual inspection of CNVs and their subsequent orthogonal validations using other techniques like Multiplex Ligation-dependent Probe Amplification (MLPA), which can be costly both economically and in terms of personnel.

 

The scope of the project include the following task:

Determination of the State-of-the-Art CNV algorithms: research and identify the most advanced and effective algorithms for CNV detection.

Obtaining validated CNVs: Gathering CNVs that have been validated orthogonally by our group to establish a ground truth for evaluating the performance of the model.

Introduction of in-silico CNVs: Generating synthetic CNVs into samples allow to create a diverse training dataset that is essential for training a robust model.

Designing the machine learning pipeline: Training the model incrementally with in-silico dataset to build a robust model for CNV detection.

Validation of the discovered CNVs to evaluate the performance of the model

 

 

Mates, J., Mademont-Soler, I., Del Olmo, B., Ferrer-Costa, C., Coll, M., Pérez-Serra, A., Picó, F., Allegue, C., Fernandez-Falgueras, A., Álvarez, P., Yotti, R., Espinosa, M. A., Sarquella-Brugada, G., Cesar, S., Carro, E., Brugada, J., Arbelo, E., Garcia-Pavia, P., Borregan, M., Tizzano, E., … Brugada, R. (2018). Role of copy number variants in sudden cardiac death and related diseases: genetic analysis and translation into clinical practice. European journal of human genetics : EJHG, 26(7), 1014–1025. https://doi.org/10.1038/s41431-018-0119-1

Moreno-Cabrera, J.M., del Valle, J., Castellanos, E. et al. Evaluation of CNV detection tools for NGS panel data in genetic diagnostics. Eur J Hum Genet 28, 1645–1655 (2020). https://doi.org/10.1038/s41431-020-0675-z

Pounraja, V. K., Jayakar, G., Jensen, M., Kelkar, N., & Girirajan, S. (2019). A machine-learning approach for accurate detection of copy number variants from exome sequencing. Genome research, 29(7), 1134–1143. https://doi.org/10.1101/gr.245928.118

Molyneaux, B. J., Goff, L. A., Brettler, A. C., Chen, H. H., Hrvatin, S., Rinn, J. L., & Arlotta, P. (2015). DeCoN: genome-wide analysis of in vivo transcriptional dynamics during pyramidal neuron fate selection in neocortex. Neuron, 85(2), 275–288. https://doi.org/10.1016/j.neuron.2014.12.024

Babadi, M., Fu, J.M., Lee, S.K. et al. GATK-gCNV enables the discovery of rare copy number variants from exome sequencing data. Nat Genet 55, 1589–1597 (2023). https://doi.org/10.1038/s41588-023-01449-0