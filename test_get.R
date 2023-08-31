#This script was added to demonstrate how the recently added get_ functions were tested with different parameter values, parameter combinations, etc.
#Considered removing this script after the PR from which it's included in is merged into master (or leave as a testing resource).

#ID_EASE
#using bundled metadata with these_samples_metadata together with these_sample_ids (1 actual sample, 3 fake samples)
id_ease_test1 = id_ease(these_samples_metadata = get_gambl_metadata(), 
                        these_sample_ids = c("r2d2","c3P0","Luke","Reddy_3812T"), 
                        verbose = TRUE)

dim(id_ease_test1) #check dimensions of returned metadata table
unique(id_ease_test1$sample_id) #how many unique sample IDs do we have?

#only providing the 3 fake samples and 1 actual sample with these_sample_ids
id_ease_test2 = id_ease(these_sample_ids = c("r2d2","c3P0","Luke","Reddy_3812T"),
                        verbose = TRUE)

dim(id_ease_test2) #check dimensions of returned metadata table
unique(id_ease_test2$sample_id) #how many unique sample IDs do we have?

#only one sample with these_sample_ids (non-existing sample ID, meant to fail)
id_ease_test3 = id_ease(these_sample_ids = c("r2d2"), 
                        verbose = TRUE)

#test if seq_type is honored
id_ease_test4 = id_ease(this_seq_type = "genome", 
                        verbose = TRUE)

unique(id_ease_test4$seq_type) #check what seq types we have in the returned metadata table

id_ease_test5 = id_ease(this_seq_type = "capture", 
                        verbose = TRUE)

unique(id_ease_test5$seq_type) #check what seq types we have in the returned metadata table

#don't give the function anything
id_ease_test6 = id_ease(verbose = FALSE)

dim(id_ease_test6) #check dimensions for the returned metadata table, with seq_type set to default (genome + capture)


#give the function a sample ID (capture) with genome metadata and the other way around (meant to fail)
id_ease_test7 = id_ease(these_sample_ids = "00-22011_tumorB", 
                        this_seq_type = "genome", 
                        verbose = TRUE)

id_ease_test8 = id_ease(these_sample_ids = "00-16220_tumorB", 
                        this_seq_type = "capture", 
                        verbose = TRUE)


##########################################################################################################################################################################################################


#GET_MANTA_SV
#non-default parameter combinations #1
manta_sv_test1 = get_manta_sv(these_sample_ids = "DOHH-2", 
                              projection = "hg38", 
                              pairing_status = "unmatched", 
                              min_vaf = 0.3)

unique(manta_sv_test1$tumour_sample_id) #what's the returned sample ID?
length(unique(manta_sv_test1$tumour_sample_id)) #how many samples do we have manta calls for?
head(manta_sv_test1$CHROM_A, 1) #check if chr are prefixed, should be, it's hg38

#non-default parameter combinations #2
manta_sv_test2 = get_manta_sv(pairing_status = "matched", 
                              min_vaf = 0.5, 
                              min_score = 100)

unique(manta_sv_test2$pair_status) #what's the pairing status of returned calls?
all(manta_sv_test2$VAF_tumour < 0.5) #do we have any VAF scores lower than 05?
all(manta_sv_test2$SCORE < 100) #do we have any variant scores less than 100?


#give it region with chr prefix when the manta calls are not prefixed in that projection and the other way around
manta_sv_test3 = get_manta_sv(projection = "grch37", 
                              region = "chr8:128723128-128774067")

head(manta_sv_test3$CHROM_A, 1) #are the chromosomes prefixed or not? Shouldn't be, we're requesting grch37

manta_sv_test4 = get_manta_sv(projection = "hg38", 
                              region = "8:127736231-127742951")

head(manta_sv_test4$CHROM_A, 1) #are the chromosomes prefixed or not? Should be, we're requesting hg38

#use the chromosome, qstart and qend parameters
manta_sv_test5 = get_manta_sv(these_sample_ids = "DOHH-2",
                              projection = "grch37", 
                              chromosome = "8", 
                              qstart = 128723128, 
                              qend = 128774067)

dplyr::select(manta_sv_test5, CHROM_A, CHROM_B, START_A, START_B, END_A, END_B) #are the returned calls within the specified region?

#use the chromosome, qstart and qend parameters, but with wrong chromosome prefix
manta_sv_test6 = get_manta_sv(these_sample_ids = "DOHH-2",
                              projection = "grch37", 
                              chromosome = "chr8", 
                              qstart = 128723128, 
                              qend = 128774067)

dplyr::select(manta_sv_test6, CHROM_A, CHROM_B) #check the if chromosomes are prefixed

#use the chromosome, qstart and qend parameters AND region (region should take precedence over the other parameters)
manta_sv_test7 = get_manta_sv(these_sample_ids = "DOHH-2",
                              projection = "grch37", 
                              chromosome = "5", 
                              qstart = 100, 
                              qend = 100000, 
                              region = "8:128723128-128774067")

dplyr::select(manta_sv_test7, CHROM_A, CHROM_B, START_A, START_B, END_A, END_B) #are the returned calls within the specified region?

#use these_samples_metadata and these_sample_ids parameters
manta_sv_test8 = get_manta_sv(these_sample_ids = "DOHH-2")

unique(manta_sv_test8$tumour_sample_id) #do we have the requested sample ID in the return?

manta_sv_test9 = get_manta_sv(these_samples_metadata = get_gambl_metadata())

length(unique(manta_sv_test9$tumour_sample_id)) #how many samples do we have Manta calls for?

#use both sample/metadata parameters (they're matching)
manta_sv_test10 = get_manta_sv(these_sample_ids = "DOHH-2", 
                              these_samples_metadata = get_gambl_metadata() %>% 
                                dplyr::filter(sample_id == "DOHH-2"))

unique(manta_sv_test10$tumour_sample_id) #do we have the requested sample ID in the return?

#use both sample/metadata parameters (they are NOT matching, should fail)
manta_sv_test11 = get_manta_sv(these_sample_ids = "DOHH-2", 
                              these_samples_metadata = get_gambl_metadata() %>% 
                                dplyr::filter(sample_id == "00-22011_tumorB"))

#provide sample ID and give it the complete metadata table and add verbose output
manta_sv_test12 = get_manta_sv(these_sample_ids = "DOHH-2", 
                               these_samples_metadata = get_gambl_metadata(), 
                               verbose = TRUE)

unique(manta_sv_test12$tumour_sample_id) #do we have the requested sample ID in the return?

#give the function nothing
manta_sv_test14 = get_manta_sv()

length(unique(manta_sv_test14$tumour_sample_id)) #how many samples do we return Manta calls for with default parameters?

#use all parameters in non-default mode
manta_sv_test15 = get_manta_sv(these_sample_ids = c("DOHH-2", "BLGSP-71-27-00422-01A-01E", "FL1008T2"),
                               these_samples_metadata = get_gambl_metadata(),
                               verbose = TRUE, 
                               projection = "hg38", 
                               chromosome = "8", 
                               qstart = 127736231, 
                               qend = 127742951, 
                               pairing_status = "unmatched", 
                               min_vaf = 0.3, 
                               min_score = 50, 
                               pass = FALSE)

#give the function non-sense parameters (meant to fail)
manta_sv_test16 = get_manta_sv(this_is_not_a_parameter = TRUE)

#give the function a parameter that is not intended for this version of the function (meant to fail)
manta_sv_test16 = get_manta_sv(write_to_file = TRUE)
                               
                               
##########################################################################################################################################################################################################


#GET_CN_SEGMENTS
#give the function nothing (should fail, no regions to return segments for)
get_cn_segments_test1 = get_cn_segments()

#provide a region
get_cn_segments_test2 = get_cn_segments(region = "8:128723128-128774067")

get_cn_segments_test2 #are the CN segments matching the requested regions?

#use chromosome, qstart and qend parameters for the same region
get_cn_segments_test3 = get_cn_segments(chromosome = "8", 
                                        qstart = 128723128, 
                                        qend = 128774067)

identical(get_cn_segments_test3, get_cn_segments_test2)  #are the two data frames identical? (they should be)

#provide a chr prefixed chromosome in region for a projection that does not have that prefix
get_cn_segments_test4 = get_cn_segments(chromosome = "chr8", 
                                        qstart = 128723128, 
                                        qend = 128774067)

head(get_cn_segments_test4$chrom, 5) #are the returned regions chr prefixed or not? They should not be, we're dealing with grch37

#try a different projection
get_cn_segments_test5 = get_cn_segments(region = "chr8:127736231-127742951", 
                                        projection = "hg38")

#try an unsupported seq type for this version of the function (meant to fail)
get_cn_segments_test6 = get_cn_segments(region = "8:128723128-128774067", 
                                        this_seq_type = "capture")

#see if we can get chr prefixes added for a projection that does not have it by default
get_cn_segments_test7 = get_cn_segments(region = "8:128723128-128774067", 
                                        with_chr_prefix = TRUE)

head(get_cn_segments_test7$chrom, 5) #are the returned regions chr prefixed or not?

#try out the streamlined option
get_cn_segments_test8 = get_cn_segments(region = "8:128723128-128774067", 
                                        streamlined = TRUE)

ncol(get_cn_segments_test8) #should be two columns...

get_cn_segments_test9 = get_cn_segments(region = "8:128723128-128774067", 
                                        streamlined = FALSE)

ncol(get_cn_segments_test9) #should be 7

#try different parameter combinations
get_cn_segments_test10 = get_cn_segments(region = "chr8:127736231-127742951", 
                                         with_chr_prefix = TRUE, 
                                         projection = "hg38")

get_cn_segments_test11 = get_cn_segments(region = "8:128723128-128774067", 
                                         streamlined = TRUE, 
                                         with_chr_prefix = TRUE)

#give the function non-existing parameters (meant to fail)
get_cn_segments_test12 = get_cn_segments(region = "8:128723128-128774067", 
                                         this_is_not_a_para = TRUE)

#give the function both a region and region in chromsome, qstart and qend format (defaults to the supplied region)
get_cn_segments_test13 = get_cn_segments(region = "8:127736231-127742951", 
                                         chromosome = "1", 
                                         qstart = 100, 
                                         qend = 100000, 
                                         projection = "hg38")

sum(str_detect(get_cn_segments_test13$chrom, '1')) > 0 #check if the returned SSMs are on chr1 (should be FALSE)
sum(str_detect(get_cn_segments_test13$chrom, '8')) > 0 #check if the returned SSMs are on chr8 (should be TRUE)

#call with all parameters
get_cn_segments_test14 = get_cn_segments(region = "chr8:127736231-127742951", 
                                         chromosome = "chr8", 
                                         qstart = 127736231, 
                                         qend = 127742951, 
                                         projection = "hg38", 
                                         with_chr_prefix = TRUE, 
                                         streamlined = TRUE, 
                                         this_seq_type = "genome")


##########################################################################################################################################################################################################


#GET_CODING_SSM
#give the function nothing
coding_ssm_test1 = get_coding_ssm()

dim(coding_ssm_test1) #check dimensions of the returned SSM calls

#test out another projection
coding_ssm_test2 = get_coding_ssm(projection = "hg38")

dim(coding_ssm_test2) #check dimensions of the returned SSM calls

#give it a metadata table with samples IDs of interest
coding_ssm_test3 = get_coding_ssm(these_samples_metadata = get_gambl_metadata() %>% 
                                    dplyr::filter(cohort == "DLBCL_cell_lines"))

unique(coding_ssm_test3$Tumor_Sample_Barcode) #what sample do we have?

#limit the cohort
coding_ssm_test4 = get_coding_ssm(limit_cohort = "DLBCL_cell_lines")

unique(coding_ssm_test4$Tumor_Sample_Barcode) #what sample do we have?
identical(coding_ssm_test3, coding_ssm_test4) #is it identical with the previous test, it should be

#exclude the same cohort
coding_ssm_test5 = get_coding_ssm(exclude_cohort = "DLBCL_cell_lines")

#limit pathology
coding_ssm_test6 = get_coding_ssm(limit_pathology = "DLBCL")

#request the DLBCL cohort, but limit to DOHH-2
coding_ssm_test7 = get_coding_ssm(limit_cohort = "DLBCL_cell_lines", 
                                  limit_samples = "DOHH-2")

unique(coding_ssm_test7$Tumor_Sample_Barcode) #do we ahve the sample we requested?

#return an object with specific columns
coding_ssm_test8 = get_coding_ssm(basic_columns = FALSE, 
                                  maf_cols = c("Hugo_Symbol", "Variant_Type", "Variant_Classification", "Tumor_Sample_Barcode"))

colnames(coding_ssm_test8) #do we have the requested columns returned?

#use the force unmatched samples parameter
coding_ssm_test9 = get_coding_ssm(force_unmatched_samples = "DOHH-2",
                                  these_samples_metadata = get_gambl_metadata() %>% 
                                    dplyr::filter(cohort == "DLBCL_cell_lines"))

#use less stringent `min read support` filter but remove silent mutations
coding_ssm_test10 = get_coding_ssm(min_read_support = 1, include_silent = FALSE)

all(coding_ssm_test10$t_alt_count < 100) #did our filters work?

#use different parameter combinations
coding_ssm_test11 = get_coding_ssm(limit_cohort = "DLBCL_cell_lines", 
                                   limit_pathology = "DLBCL", 
                                   limit_samples = c("OCI-Ly10", "DOHH-2"), 
                                   these_samples_metadata = GAMBLR.data::gambl_metadata, 
                                   force_unmatched_samples = "DOHH-2", 
                                   projection = "hg38", 
                                   basic_columns = FALSE, 
                                   maf_cols = c("Hugo_Symbol", "Variant_Type", "Variant_Classification", "Tumor_Sample_Barcode"), 
                                   min_read_support = 1, 
                                   include_silent = FALSE)

#try a non-existing parameter (meant to fail)
coding_ssm_test12 = get_coding_ssm(this_is_not_a_para = TRUE)


##########################################################################################################################################################################################################


#GET_SAMPLE_SEGMENTS
#give the function nothing
sample_segs_test1 = get_sample_cn_segments()

dim(sample_segs_test1) #what does the return look like?

#use these_samples_metadata and these_sample_ids parameters
sample_segs_test2 = get_sample_cn_segments(these_sample_ids = "DOHH-2")

unique(sample_segs_test2$ID) #do we have the requested sample ID?

sample_segs_test3 = get_sample_cn_segments(these_samples_metadata = get_gambl_metadata())

length(unique(sample_segs_test3$ID)) #how many samples are avaialble?

#use both sample/metadata parameters (they're matching)
sample_segs_test4 = get_sample_cn_segments(these_sample_ids = "DOHH-2", 
                                           these_samples_metadata = get_gambl_metadata() %>% 
                                              dplyr::filter(sample_id == "DOHH-2"))

unique(sample_segs_test4$ID) #do we have the requested sample ID?

#use both sample/metadata parameters (they are NOT matching, should fail)
sample_segs_test5 = get_sample_cn_segments(these_sample_ids = "DOHH-2", 
                                           these_samples_metadata = GAMBLR.data::gambl_metadata %>% 
                                              dplyr::filter(sample_id == "00-22011_tumorB"))

#try non-default projection
sample_segs_test6 = get_sample_cn_segments(projection = "hg38")

#try a non-supported seq type (meant to fail)
sample_segs_test7 = get_sample_cn_segments(this_seq_type = "capture")

#try streamlined
sample_segs_test8 = get_sample_cn_segments(streamlined = TRUE)

colnames(sample_segs_test8) #do we get thhe expected columns? (ID, CN)

#try with chr prefix, for a projection that does not have it
sample_segs_test9 = get_sample_cn_segments(with_chr_prefix = TRUE)

sample_segs_test9[1,2] #are they prefixed?

#try a non-existing parameter (meant to fail)
coding_ssm_test10 = get_sample_cn_segments(this_is_not_a_para = TRUE)

#pass multiple non-default parameters to the function
coding_ssm_test10 = get_sample_cn_segments(these_sample_ids = "DOHH-2",
                                           these_samples_metadata = get_gambl_metadata(), 
                                           projection = "hg38", 
                                           this_seq_type = "genome", 
                                           with_chr_prefix = TRUE)


##########################################################################################################################################################################################################


#GET_SSM_BY_SAMPLE/SAMPLES
#give it nothing!
sample_ssm_test1 = get_ssm_by_samples()

dim(sample_ssm_test1) #dimensions?

#test these_sample_ids parameter (notice the alias function used here...)
sample_ssm_test2 = get_ssm_by_sample(these_sample_ids = "DOHH-2")

unique(sample_ssm_test2$Tumor_Sample_Barcode) #do we get the sample we want?

#test these_sample_ids parameter with multiple sample IDs
sample_ssm_test3 = get_ssm_by_samples(these_sample_ids = c("DOHH-2", "OCI-Ly10"))

unique(sample_ssm_test3$Tumor_Sample_Barcode) #do we get the samples we want?

#test these_samples_metadata parameter
sample_ssm_test4 = get_ssm_by_samples(these_samples_metadata = get_gambl_metadata() %>% 
                                        dplyr::filter(sample_id == "DOHH-2"))

unique(sample_ssm_test4$Tumor_Sample_Barcode) #do we get the sample we want?

#try a non-valid seq type (meant to fail)
sample_ssm_test5 = get_ssm_by_samples(this_seq_type = "capture")

#try different gene formats
sample_ssm_test6 = get_ssm_by_samples(these_genes = "MYC")

unique(sample_ssm_test6$Hugo_Symbol) #do we get MYC and MYC only?

sample_ssm_test7 = get_ssm_by_samples(these_genes = c("MYC", "BCL2"))

unique(sample_ssm_test7$Hugo_Symbol) #do we get the genes we requested?

#request specific MAF columns
sample_ssm_test8 = get_ssm_by_samples(basic_columns = FALSE, 
                                      maf_cols = c("Hugo_Symbol", "Variant_Type", "Variant_Classification", "Tumor_Sample_Barcode"))

colnames(sample_ssm_test8) #do we get the requested columns?

#parameter combination tests
sample_ssm_test9 = get_ssm_by_samples(these_sample_ids = "DOHH-2", 
                                      these_samples_metadata = GAMBLR.data::gambl_metadata, 
                                      this_seq_type = "genome", 
                                      projection = "hg38", 
                                      these_genes = c("MYC", "BCL2"), 
                                      min_read_support = 1, 
                                      verbose = TRUE, 
                                      basic_columns = FALSE, 
                                      maf_cols = c("Hugo_Symbol", "Variant_Type", "Variant_Classification", "Tumor_Sample_Barcode"))

#try a non-existing parameter (meant to fail)
sample_ssm_test10 = get_ssm_by_samples(this_is_not_a_para = TRUE)


##########################################################################################################################################################################################################


#GET_SSM_BY_PATIENTS
#get ssm for one patient and compare to the same return but using get_ssm_by_sample
patient_ssm_test1 = get_ssm_by_patients(these_patient_ids = "DOHH-2")
patient_ssm_test2 = get_ssm_by_sample(these_sample_ids = "DOHH-2")

identical(patient_ssm_test1, patient_ssm_test2) #are they the same (they should be)?

#try a non-valid seq type (meant to fail)
patient_ssm_test3 = get_ssm_by_patients(these_patient_ids = "DOHH-2", projection = "hg38", this_seq_type = "capture")

#non-default projection
patient_ssm_test4 = get_ssm_by_patients(these_patient_ids = "DOHH-2", projection = "hg38")

identical(patient_ssm_test1, patient_ssm_test4) #are they the same, they should not be!

#use these_samples_metadata parameter
patient_ssm_test5 = get_ssm_by_patients(these_samples_metadata = get_gambl_metadata() %>% 
                                          dplyr::filter(cohort == "DLBCL_cell_lines"))

unique(patient_ssm_test5$Tumor_Sample_Barcode) #what samples are returned?

#test multiple samples as a vector with these_sample_ids
patient_ssm_test6 = get_ssm_by_patients(these_patient_ids = c("DOHH-2", "OCI-Ly10"))

unique(patient_ssm_test6$Tumor_Sample_Barcode) #what samples are returned?

#try some of the filtering options
patient_ssm_test7 = get_ssm_by_patients(these_patient_ids = c("DOHH-2", "OCI-Ly10"), 
                                        min_read_support = 50)

all(patient_ssm_test7$t_alt_count < 50) #do we have any variant with read count less than 50? (i.e did our filter work)

#return a set number of columns in the MAF
patient_ssm_test8 = get_ssm_by_patients(these_patient_ids = c("DOHH-2", "OCI-Ly10"), 
                                        basic_columns = FALSE, 
                                        min_read_support = 20, 
                                        maf_cols = c("Variant_Classification", "Tumor_Sample_Barcode"))

colnames(patient_ssm_test8) #do we get the requested columns?

#give the function sample IDs and a metadata table (matching)
cell_lines_patient = get_gambl_metadata() %>% 
  dplyr::filter(cohort == "DLBCL_cell_lines") %>%
  pull(patient_id)

cell_lines_meta = get_gambl_metadata() %>% 
  dplyr::filter(cohort == "DLBCL_cell_lines")

patient_ssm_test9 = get_ssm_by_patients(these_samples_metadata = cell_lines_meta, these_patient_ids = cell_lines_patient)

unique(patient_ssm_test9$Tumor_Sample_Barcode) #what samples are kept?

#give the function sample IDs and a metadata table (not matching) +verbose (should be 0 variants)
some_meta = get_gambl_metadata() %>% 
  dplyr::filter(sample_id == "00-16220_tumorB")

patient_ssm_test10 = get_ssm_by_patients(these_samples_metadata = some_meta, 
                                         these_patient_ids = cell_lines_patient, 
                                         verbose = TRUE)

#give the function non-existing parameters (meant to fail)
patient_ssm_test11 = get_ssm_by_patients(these_patient_idsss = "DOHH-2",
                                         this_is_not_a_parameter = NULL, 
                                         jimi_hendrix_is_god = TRUE)


##########################################################################################################################################################################################################


#GET_SSM_BY_REGION
#use region parameter
region_ssm_test1 = get_ssm_by_region(region = "8:128,723,128-128,774,067")
region_ssm_test2 = get_ssm_by_region(region = "chr8:127736231-127742951", projection = "hg38", verbose = FALSE)

#use chromosome, qstart and qend instead
region_ssm_test3 = get_ssm_by_region(chromosome = "chr8", qstart = 128723128, qend = 128774067)

identical(region_ssm_test1, region_ssm_test3) #are the return the same for the two calls using the different parameters? They should be

#test streamlined option
region_ssm_test4 = get_ssm_by_region(region = "8:128,723,128-128,774,067", streamlined = TRUE)

colnames(region_ssm_test4) #did streamline work? (i.e Start_Position and Tumor_Sample_Barcode as the only two columns)

#test min read filtering options
region_ssm_test5 = get_ssm_by_region(region = "8:128,723,128-128,774,067", min_read_support = 10)

all(region_ssm_test5$t_alt_count < 10) #any variants with a read support of less than 10? (should be FALSE if the filtering worked)

#give the function more than one region (meant to fail and suggest user to call get_ssms_by_regions instead)
region_ssm_test6 = get_ssm_by_region(region = c("8:128,723,128-128,774,067", "chr8:127736231-127742951"))


##########################################################################################################################################################################################################


#GET_SSM_BY_REGIONS
#get some regions in different ways
#get regions in a bed format and add a region name column
this_bed = dplyr::mutate(GAMBLR.data::grch37_ashm_regions, name = paste(gene, region, sep = "_"))

#use regions bed
regions_ssm_test1 = get_ssm_by_regions(regions_bed = this_bed)

#use a vector of regions
regions_ssm_test2 = get_ssm_by_regions(regions_vector = c("chr1:6661482-6662702", 
                                                          "chr2:232572640-232574297", 
                                                          "chr17:75443766-75451177"), streamlined = TRUE)

colnames(regions_ssm_test2) #do we get the streamlined output (start, sample_id and region_name)

#return all MAF columns
regions_ssm_test3 = get_ssm_by_regions(regions_vector = c("chr1:6661482-6662702", 
                                                          "chr2:232572640-232574297", 
                                                          "chr17:75443766-75451177"), 
                                       basic_columns = TRUE)

colnames(regions_ssm_test3) #do we get all the MAF columns returned?

#test non-existing parameter (meant to fail)
regions_ssm_test4 = get_ssm_by_regions(region_bed = this_bed, streamlined = FALSE, use_name_column = TRUE)
