
#https://stackoverflow.com/questions/48262384/increase-r-studios-connection-timeout-limit
options(timeout = max(4000000, getOption("timeout")))

#https://support.bioconductor.org/p/9139740/
#https://www.biostars.org/p/429062/
library(biomaRt)
packageVersion("biomaRt") #‘2.50.3’

biolist <- as.data.frame(listMarts()) 
biolist #biolist contain version information
#biomart                version
#1 ENSEMBL_MART_ENSEMBL      Ensembl Genes 109
#2   ENSEMBL_MART_MOUSE      Mouse strains 109
#3     ENSEMBL_MART_SNP  Ensembl Variation 109
#4 ENSEMBL_MART_FUNCGEN Ensembl Regulation 109

#Every analysis with biomaRt starts with selecting a BioMart database to use.
#https://bioconductor.org/packages/release/bioc/vignettes/biomaRt/inst/doc/accessing_ensembl.html

##########################################################
## retrieve human-mouse orthologlous gene dn ds values
if(F){
    # 1. Selecting an Ensembl BioMart database and dataset
    # step1, identifying the database you need
    listEnsembl()
    ensembl <- useEnsembl(biomart = "genes")
    ensembl
    # step2, Choosing a dataset
    datasets <- listDatasets(ensembl)
    head(datasets)
    datasets[grep('musculus',datasets$dataset),]
    #dataset                    description     version
    #18  bmusculus_gene_ensembl Blue whale genes (mBalMus1.v2) mBalMus1.v2
    #107 mmusculus_gene_ensembl           Mouse genes (GRCm39)      GRCm39
    datasets[grep('sapien',datasets$dataset),]
    #dataset              description    version
    #80 hsapiens_gene_ensembl Human genes (GRCh38.p14) GRCh38.p14
    
    #Using archived versions of Ensembl
    listEnsemblArchives()
}

# as a used case show here: https://support.bioconductor.org/p/9141035/
## For the human data 
ensemblhsapiens = useEnsembl(version = 99, #this archived still contain dn ds values
                             biomart = 'ENSEMBL_MART_ENSEMBL', 
                             dataset = 'hsapiens_gene_ensembl')
ensemblhsapiens
#Object of class 'Mart':
#Using the ENSEMBL_MART_ENSEMBL BioMart database
#Using the hsapiens_gene_ensembl dataset

#https://rdrr.io/bioc/biomaRt/man/listAttributes.html
listAttributes(ensemblhsapiens) 
searchAttributes(ensemblhsapiens, 'homolog_dn') #yes, there is xxx_dn

## Input human gene list
hsapiens_GFList <- c("ENSG00000211979", "ENSG00000224650", "ENSG00000211973")

## human and mouse
hsapiens_mouse <- getBM(attributes = c('ensembl_gene_id', 
                                     'mmusculus_homolog_ensembl_gene',  
                                     'mmusculus_homolog_dn', 
                                     'mmusculus_homolog_ds',
                                     'mmusculus_homolog_orthology_type',
                                     'mmusculus_homolog_orthology_confidence'), 
                      filters = 'ensembl_gene_id', 
                      values = hsapiens_GFList, 
                      mart = ensemblhsapiens)
hsapiens_mouse


#quote"However I think it's telling that support for the dN/dS analysis has been dropped for over a year now, 
#from Ensembl 100 onwards (https://www.ensembl.info/2020/04/29/ensembl-100-has-been-released/)."


# retrieve dn/ds between M.musculus and the rat Rattus norvegicus using their ensembly genomes
#https://www.biostars.org/p/147351/
mouse= useEnsembl(version = 99, #this archived still contain dn ds values
                             biomart = 'ENSEMBL_MART_ENSEMBL', 
                             dataset = 'mmusculus_gene_ensembl')

listAttributes(mouse) 
searchAttributes(mouse,'hsapiens_homolog')
searchAttributes(mouse, 'hsapiens_homolog_dn') #yes, there is xxx_dn
searchAttributes(mouse, 'hsapiens_homolog_ds') #yes, there is xxx_dn


##################################################
## extract gene ID
## mouse
ensembl=useMart("ensembl") 
esemblist <- as.data.frame(listDatasets(ensembl))
esemblist[grep('musculus',esemblist$dataset),]
#dataset                              description  version
#mmusculus_gene_ensembl           Mouse genes (GRCm39)      GRCm39

mouse= useEnsembl(#version = 99, #this archived still contain dn ds values
                  version = 109,
                  biomart = 'ENSEMBL_MART_ENSEMBL', 
                  dataset = 'mmusculus_gene_ensembl')
filters = listFilters(mouse)
attributes = listAttributes(mouse)

annotLookup <- getBM(
  mart = mouse,
  attributes = c(
    'ensembl_gene_id',
    'gene_biotype',
    'external_gene_name',#'external_synonym','source_name',
    'uniprot_gn_symbol',
    'uniprot_gn_id',
    'mgi_symbol','description',
    'chromosome_name', 'start_position','end_position'),
  uniqueRows=TRUE)
head(annotLookup)

#head(subset(annotLookup, uniprot_gn_id != ''), 20)[,-4]
annotLookup2=annotLookup[annotLookup$uniprot_gn_id!='',]
dim(annotLookup) #86873
length(unique(annotLookup2$uniprot_gn_id)) #51102

saveRDS(annotLookup,'mmusculus_id_v109.rds')

######################
## rat
ensembl=useMart("ensembl") 
esemblist <- as.data.frame(listDatasets(ensembl))
esemblist[grep('norvegicus',esemblist$dataset),]
#dataset                              description  version
#rnorvegicus_gene_ensembl Rat genes (mRatBN7.2) mRatBN7.2

rat= useEnsembl(#version = 99, #this archived still contain dn ds values
  version = 109,
  biomart = 'ENSEMBL_MART_ENSEMBL', 
  dataset = 'rnorvegicus_gene_ensembl')
filters = listFilters(rat)
attributes = listAttributes(rat)
attributes[attributes$page=='feature_page',]

annotLookup <- getBM(
  mart = rat,
  attributes = c(
    'ensembl_gene_id',
    'gene_biotype',
    'external_gene_name','external_synonym','source_name',
    'uniprot_gn_symbol',
    'uniprot_gn_id',
    'mgi_symbol','description',
    'chromosome_name', 'start_position','end_position'),
  uniqueRows=TRUE)
head(annotLookup)

#head(subset(annotLookup, uniprot_gn_id != ''), 20)[,-4]
annotLookup2=annotLookup[annotLookup$uniprot_gn_id!='',]
dim(annotLookup) #50429
length(unique(annotLookup2$uniprot_gn_id)) #14983

saveRDS(annotLookup,'rnorvegicus_id_v109.rds')


######################
## human
ensembl=useMart("ensembl") 
esemblist <- as.data.frame(listDatasets(ensembl))
esemblist[grep('sapiens',esemblist$dataset),]
#dataset                              description  version
#hsapiens_gene_ensembl Human genes (GRCh38.p13) GRCh38.p13

human= useEnsembl(#version = 99, #this archived still contain dn ds values
  version = 109,
  biomart = 'ENSEMBL_MART_ENSEMBL', 
  dataset = 'hsapiens_gene_ensembl')
filters = listFilters(human)
attributes = listAttributes(human)
attributes[attributes$page=='feature_page',]

annotLookup <- getBM(
  mart = human,
  attributes = c(
    'ensembl_gene_id',
    'gene_biotype',
    'external_gene_name','external_synonym','source_name',
    'uniprot_gn_symbol',
    'uniprot_gn_id',
    #'mgi_symbol',
    'description',
    'chromosome_name', 'start_position','end_position'),
  uniqueRows=TRUE)
head(annotLookup)

#head(subset(annotLookup, uniprot_gn_id != ''), 20)[,-4]
annotLookup2=annotLookup[annotLookup$uniprot_gn_id!='',]
dim(annotLookup) #433733
length(unique(annotLookup2$uniprot_gn_id)) #71809

saveRDS(annotLookup,'hsapiens_id_v109.rds')

######################
## mouse
ensembl=useMart("ensembl") 
esemblist <- as.data.frame(listDatasets(ensembl))
esemblist[grep('musculus',esemblist$dataset),]
#dataset                              description  version
#mmusculus_gene_ensembl           Mouse genes (GRCm39)      GRCm39

mouse= useEnsembl(#version = 99, #this archived still contain dn ds values
                  version = 109,
                  biomart = 'ENSEMBL_MART_ENSEMBL', 
                  dataset = 'mmusculus_gene_ensembl')
filters = listFilters(mouse)
attributes = listAttributes(mouse)

annotLookup <- getBM(
  mart = mouse,
  attributes = c(
    'ensembl_gene_id',
    'gene_biotype',
    'external_gene_name',#'external_synonym','source_name',
    'uniprot_gn_symbol',
    'uniprot_gn_id',
    'mgi_symbol','description',
    'chromosome_name', 'start_position','end_position'),
  uniqueRows=TRUE)
head(annotLookup)

#head(subset(annotLookup, uniprot_gn_id != ''), 20)[,-4]
annotLookup2=annotLookup[annotLookup$uniprot_gn_id!='',]
dim(annotLookup) #86873
length(unique(annotLookup2$uniprot_gn_id)) #51102

saveRDS(annotLookup,'mmusculus_id_v109.rds')

######################
## rat
ensembl=useMart("ensembl") 
esemblist <- as.data.frame(listDatasets(ensembl))
esemblist[grep('norvegicus',esemblist$dataset),]
#dataset                              description  version
#rnorvegicus_gene_ensembl Rat genes (mRatBN7.2) mRatBN7.2

rat= useEnsembl(#version = 99, #this archived still contain dn ds values
  version = 109,
  biomart = 'ENSEMBL_MART_ENSEMBL', 
  dataset = 'rnorvegicus_gene_ensembl')
filters = listFilters(rat)
attributes = listAttributes(rat)
attributes[attributes$page=='feature_page',]

annotLookup <- getBM(
  mart = rat,
  attributes = c(
    'ensembl_gene_id',
    'gene_biotype',
    'external_gene_name','external_synonym','source_name',
    'uniprot_gn_symbol',
    'uniprot_gn_id',
    'mgi_symbol','description',
    'chromosome_name', 'start_position','end_position'),
  uniqueRows=TRUE)
head(annotLookup)

#head(subset(annotLookup, uniprot_gn_id != ''), 20)[,-4]
annotLookup2=annotLookup[annotLookup$uniprot_gn_id!='',]
dim(annotLookup) #50429
length(unique(annotLookup2$uniprot_gn_id)) #14983

saveRDS(annotLookup,'rnorvegicus_id_v109.rds')


######################
## human
ensembl=useMart("ensembl") 
esemblist <- as.data.frame(listDatasets(ensembl))
esemblist[grep('sapiens',esemblist$dataset),]
#dataset                              description  version
#hsapiens_gene_ensembl Human genes (GRCh38.p13) GRCh38.p13

human= useEnsembl(#version = 99, #this archived still contain dn ds values
  version = 109,
  biomart = 'ENSEMBL_MART_ENSEMBL', 
  dataset = 'hsapiens_gene_ensembl')
filters = listFilters(human)
attributes = listAttributes(human)
attributes[attributes$page=='feature_page',]

annotLookup <- getBM(
  mart = human,
  attributes = c(
    'ensembl_gene_id',
    'gene_biotype',
    'external_gene_name','external_synonym','source_name',
    'uniprot_gn_symbol',
    'uniprot_gn_id',
    #'mgi_symbol',
    'description',
    'chromosome_name', 'start_position','end_position'),
  uniqueRows=TRUE)
head(annotLookup)

#head(subset(annotLookup, uniprot_gn_id != ''), 20)[,-4]
annotLookup2=annotLookup[annotLookup$uniprot_gn_id!='',]
dim(annotLookup) #433733
length(unique(annotLookup2$uniprot_gn_id)) #71809

saveRDS(annotLookup,'hsapiens_id_v109.rds')

############################################################
## use biomaRt to get gene chr info 
# https://github.com/mingwhy/bioinfo_homemade_tools/tree/main/used.cases/biomaRt_usage
# add gene_biotype info: https://support.bioconductor.org/p/62441/
if(F){ #run once
#https://stackoverflow.com/questions/13012210/find-transcription-start-sites-with-biomart
library(biomaRt)
packageVersion("biomaRt") #‘2.50.3’

biolist <- as.data.frame(listMarts())
ensembl=useMart("ensembl")
esemblist <- as.data.frame(listDatasets(ensembl))
esemblist[grep('melanogaster',esemblist$description),]
#dataset                              description  version
#56 dmelanogaster_gene_ensembl Drosophila melanogaster genes (BDGP6.32) BDGP6.32

ensembl = useDataset("dmelanogaster_gene_ensembl",mart=ensembl)
filters = listFilters(ensembl)
attributes = listAttributes(ensembl)

attributes[grep('dmel',attributes$name),]
grep('transcript',attributes$name,value=TRUE,ignore.case = TRUE)

t2g<-getBM(attributes=c('ensembl_gene_id',"gene_biotype",
                        "ensembl_transcript_id","transcript_start","transcript_end",
                        #"ensembl_exon_id","exon_chrom_start","exon_chrom_end",
                        "chromosome_name","strand","transcript_biotype", 
                        'start_position','end_position',
                        'flybase_transcript_id',"transcript_length","transcript_count"), mart = ensembl)
dim(t2g) #41209    13
head(t2g)
saveRDS(t2g,'t2g_chr.coord.rds')
}
t2g=readRDS('t2g_chr.coord.rds')
# one gene has multiple rows as muliple transcripts
df.gene.length=t2g[!duplicated(t2g$ensembl_gene_id),]
table(df.gene.length$gene_biotype)
df.gene.length=df.gene.length[df.gene.length$gene_biotype == 'protein_coding',]
dim(df.gene.length) #13968    13


#####################################################################
## id conversion between ensemble and uniprot 
#https://www.biostars.org/p/429062/
library(biomaRt)
packageVersion("biomaRt") #‘2.50.3’

biolist <- as.data.frame(listMarts())
ensembl=useMart("ensembl")
esemblist <- as.data.frame(listDatasets(ensembl))
esemblist[grep('musculus',esemblist$dataset),]
#dataset                              description  version
#108 mmusculus_gene_ensembl           Mouse genes (GRCm39)      GRCm39

ensembl = useDataset("mmusculus_gene_ensembl",mart=ensembl)
filters = listFilters(ensembl)
attributes = listAttributes(ensembl)

annotLookup <- getBM(
  mart = ensembl,
  attributes = c(
    'ensembl_gene_id',
    'gene_biotype',
    'external_gene_name',
    'uniprot_gn_symbol',
    'uniprot_gn_id'),
  uniqueRows=TRUE)
head(annotLookup)
head(subset(annotLookup, uniprot_gn_id != ''), 20)[,-4]

annotLookup2=annotLookup[annotLookup$uniprot_gn_id!='',]
dim(annotLookup2) #51706
length(unique(annotLookup2$uniprot_gn_id)) #51102

saveRDS(annotLookup2,'mmus_id_ensembl2uniprot.rds')

###############################################################
## use biomaRt to get protein sequence
# https://github.com/mingwhy/bioinfo_homemade_tools/tree/main/used.cases/biomaRt_usage
# add gene_biotype info: https://support.bioconductor.org/p/62441/
if(F){ #run once
  #https://stackoverflow.com/questions/13012210/find-transcription-start-sites-with-biomart
  library(biomaRt)
  packageVersion("biomaRt") #‘2.50.3’
  
  biolist <- as.data.frame(listMarts())
  ensembl=useMart("ensembl")
  esemblist <- as.data.frame(listDatasets(ensembl))
  esemblist[grep('musculus',esemblist$dataset),]
  #dataset                              description  version
  #108 mmusculus_gene_ensembl           Mouse genes (GRCm39)      GRCm39
  
  ensembl = useDataset("mmusculus_gene_ensembl",mart=ensembl)
  filters = listFilters(ensembl)
  attributes = listAttributes(ensembl)
  
  attributes[grep('mmusculus',attributes$name),]
  attributes[grep('symbol',attributes$name),]
  grep('transcript',attributes$name,value=TRUE,ignore.case = TRUE) #feature_page
  attributes[grep('ensembl_transcript_id',attributes$name),]
  grep('peptide',attributes$name,value=TRUE,ignore.case = TRUE) #feature_page
  attributes[grep('cds_start',attributes$name),] #not in feature_page
  attributes[grep('exon_chrom_start',attributes$name),] #not in feature_page
  
  #https://www.biostars.org/p/9544009/#9544241
  ensemble_data<-getBM(attributes=c('mgi_symbol','ensembl_gene_id',
                          "ensembl_transcript_id","ensembl_peptide_id",
                          "strand","gene_biotype","ensembl_exon_id"),
                          #"cds_start","cds_end","exon_chrom_start","exon_chrom_end"),
                          mart = ensembl)
  
  dim(ensemble_data) #868930     7
  head(ensemble_data)
  
  table(ensemble_data$gene_biotype)
  protein_coding_genes<-ensemble_data[ensemble_data$gene_biotype=='protein_coding',]
  dim(protein_coding_genes) #759189
  length(unique(protein_coding_genes$mgi_symbol)) #21739
  length(unique(protein_coding_genes$ensembl_peptide_id)) #66310
  
  tmp=protein_coding_genes[protein_coding_genes$ensembl_peptide_id %in% 
    protein_coding_genes[duplicated(protein_coding_genes$ensembl_peptide_id),]$ensembl_peptide_id,]
  head(tmp)
  
  #https://rdrr.io/bioc/biomaRt/man/getSequence.html
  seq = getSequence(id = protein_coding_genes$mgi_symbol[1:2], 
                    type = "mgi_symbol", 
                    seqType = "peptide", 
                    mart = ensembl)
  show(seq)
  length(unique(protein_coding_genes$mgi_symbol)) #21739
  peptide_seqs <- getSequence(id = unique(protein_coding_genes$mgi_symbol),
                              type = "mgi_symbol", 
                              seqType = "peptide", 
                              mart = ensembl)
  dim(peptide_seqs) #69799
  
  saveRDS(peptide_seqs, 'mmusculus_peptide_seqs.rds')
}

peptide_seqs=readRDS('mmusculus_peptide_seqs.rds')


