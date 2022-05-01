#
# Cancer Biomarkers database
#
# https://www.cancergenomeinterpreter.org/biomarkers
#

The Cancer Biomarkers database is curated and maintained by several clinical and scientific experts [https://www.cancergenomeinterpreter.org/curators] 
in the field of precision oncology supported by the European Union’s Horizon 2020 funded project [http://www.medbioinformatics.eu/]. This database
 is currently being integrated with databases developed by other institutions in a collaborative effort [http://cancervariants.org/] under the 
umbrella of the Global Alliance for Genomics and Health [https://genomicsandhealth.org/]. The contribution of the community to the resource 
--both, in annotation and curation-- is encouraged. Any user may suggest new entries, edition or comments about existing entries by contacting 
us here or by using the feedback icon located at the left of each entry of the table. This table provides a summary of the content of the 
database that can be interactively browsed. Additional information, including the genomic coordinates of the variants, can be accessed by 
downloading it. This database is licensed under a Creative Commons Public Domain Dedication (CC0 1.0 Universal) [https://creativecommons.org/publicdomain/zero/1.0/]. 
When using this database, please cite: Cancer Genome Interpreter Annotates The Biological And Clinical Relevance Of Tumor Alterations; doi: https://doi.org/10.1101/140475.


Description
--------------------------
The Cancer Biomarkers database is provided through two files. The first (cgi_biomarkers.tsv) is a simplified version in which variants 
that constitute biomarkers of the same type of response response to the same drug in a given cancer type(s) are grouped in the same row. 
In the second file (cgi_biomarkers_per_variant.tsv), each variant is contained in a separate row and --for point mutations-- the genomic 
coordinates are included. For the latter, note that the correspondence between genomic and protein coordinates is computed using the gene’s 
longest transcript (except for a subset of 109 genes in which the annotation transcript was manually selected) see CGI paper for the gene 
transcript mapping. 

cgi_biomarkers.tsv
------------------
Each row describes a biomarker of drug response in a given cancer type(s).

File columns:

Biomarker: Alteration that has been linked experimentally or through clinical observation to a given drug response
Gene: Name of gene bearing the alteration 
Alteration type: type of alteration that constitutes the biomarker
Alteration: Representation of the alteration in format gene:alteration. Note that  any oncogenic mutation affecting a given gene is represented as a dot. 
Targeting (deprecated): Whether the drug targets the gene directly or indirectly
Drug status (deprecated): Stage of drug development (e.g., approved, in clinical trials, etc)
Drug family (deprecated): Family of the drug (e.g. BRAF inhibitor)
Drug: drug name
Association: Effect of the biomarker on the effect of the drug on a tumor (e.g., resistant, responsive, etc.)
Evidence level: Degree of clinical support of the effect of the biomarker (e.g., guidelines, clinical trials, etc)
Assay type: Experimental system in which biomarkers pre-clinical support level have been observed.
Source: Source of the biomarker (e.g., FDA guidelines, PMID of publication, etc)
Curator: Name of the curator(s) of the biomarker in the database
Curation date: Latest date of curation 
Primary Tumor type (Cancer type): Primary tumor type(s) in which the biomarker has been observed
Metastatic Tumor Type (deprecated): Metastatic tumor type(s) in which the biomarker has been observed
TCGI included:
Drug full name: Full name of the drug
Primary Tumor type full name: Full name of the primary tumor type(s) in which the biomarker has been observed



cgi_biomarkers_per_variant.tsv
------------------------------
Describes all variants observed in tumors that fulfill the criteria of the biomarkers of the Cancer Biomarkers Database.

Columns:

Alteration: Representation of the alteration in format gene:alteration (i.e., any oncogenic mutation affecting a given gene is represented as a dot)
Alteration type: type of alteration that constitutes the biomarker
Assay type: Experimental system in which biomarkers pre-clinical support level have been observed
Association: Effect of the biomarker on the effect of the drug on a tumor (e.g., resistant, responsive, etc.)
Biomarker: Alteration that has been linked experimentally or through clinical observation to a given drug response
Curator: Name of the curator(s) of the biomarker in the database
Drug: drug name
Drug family: Family of the drug (e.g. BRAF inhibitor)
Drug full name: Full name of the drug
Drug status: Stage of drug development (e.g., approved, in clinical trials, etc)
Evidence level: Degree of clinical support of the effect of the biomarker (e.g., guidelines, clinical trials, etc)
Gene: Name of gene that suffers the alteration
Metastatic Tumor Type: Metastatic tumor type(s) in which the biomarker has been observed
TCGI included
Primary Tumor type: Primary tumor type(s) in which the biomarker has been observed
Primary Tumor acronym:  Name of primary tumor type(s) in which the biomarker has been observed
Source: Source of the biomarker (e.g., FDA guidelines, PMID of publication, etc)
Targeting: Whether the drug targets the gene directly or indirectly
Individual_mutation: Individual mutations observed in tumors that fulfill the criteria of the biomarker
Transcript: Transcript used to map the variant to protein coordinates
Gene: Gene name
Strand: Strand of the variant
Region: region of the gene where the variant is located
Info: Information on the consequence of the variant according to Ensembl 
cDNA: location of the variant in cDNA coordinates and nucleotide substitution
gDNA: location of the variant in genomic coordinates and nucleotide substitution

