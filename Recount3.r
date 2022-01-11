


set.seed("1234")
suppressMessages(library("recount3"))
suppressMessages(library("argparse"))

parser <- ArgumentParser()

parser$add_argument("--Project",
    help="Provide project names; options are'srid': Shoud start with SRP or tissue of tcga or GTEX")

parser$add_argument("--AnnotationDb",
    help="Provide Annotation database names; options are gencode_v29, gencode_v26")

parser$add_argument("--Feature",
    help="Provide featre yu want to analyse; options are 'exon','jxn','gene'")

parser$add_argument("--DataSource",
    help="Provide data source, if sr id provided tis option dont provde; options are 'gtex','tcga'")

parser$add_argument("--SampleName",
    help="A file with sample names to be included, column name shoud be: SampleName; No other colum is required")





args <- parser$parse_args()


Project=args$Project
Gannotation=args$AnnotationDb
Type=args$Feature
DataSource=args$DataSource
SampleNames=args$SampleName


#Project="SRP009615" 
#Gannotation="gencode_v29"
#Type="gene"
#DataSource="gtex"
#SampleNames="Samle_NamesToinclude.txt"



## Find all available human projects
human_projects <- available_projects()

write.csv(human_projects, 'Human_projects_Recount.csv')


TCGA<-c("ACC" ,"BLCA","BRCA","CESC","CHOL","COAD","DLBC","ESCA","GBM","HNSC",
        "KICH", "KIRC", "KIRP", "LAML", "LGG","LIHC","LUAD","LUSC","MESO","OV",
        "PAAD","PCPG","PRAD","READ","SARC","SKCM","STAD","TGCT","THCA","THYM")

GTEX<-c("ADIPOSE_TISSUE","MUSCLE","BLOOD_VESSEL","HEART","OVARY","UTERUS","VAGINA","BREAST",
     "SKIN","SALIVARY_GLAND","BRAIN","ADRENAL_GLAND","THYROID","LUNG","SPLEEN","PANCREAS",
     "ESOPHAGUS","STOMACH","COLON","SMALL_INTESTINE","PROSTATE","TESTIS","NERVE","PITUITARY",
     "BLOOD","LIVER","KIDNEY","CERVIX_UTERI","FALLOPIAN_TUBE","BLADDER","STUDY_NA","BONE_MARROW")  




## Find the project you are interested in

if   (grepl("SRP", Project)) {

#SRID
    proj_info <- subset(
                        human_projects,
                        project == Project & project_type == "data_sources")

#TCGA  
          }else if (Project %in% TCGA) {

                if (grepl("tcga", DataSource)) {
                    proj_info<-subset(human_projects, file_source == "tcga" & project_type == "data_sources" & project==Project)
                    }else {
                    print("Please provide tcga tissue name proerly; you requested samples from tcga project")
                    quit()
                        }
#GTEX
         }else if (Project %in% GTEX) {

          if (grepl("gtex", DataSource)) {
              proj_info<-subset(human_projects, file_source == "gtex" & project_type == "data_sources" & project==Project)

          } else {
                    print("Please provide gex tissue name proerly; you requested samples from gtex project")
                    quit()
                 }
         }


## Create a RangedSummarizedExperiment (RSE) object at the gene level
rse <- create_rse(proj_info,
                   annotation = Gannotation,
                    type = Type
                    )

metadata(rse)

#PROJECT SMAMRY
ProjectDetails<-as.data.frame(colData(rse))
write.csv(ProjectDetails,paste0(Project,"_ProjectSummary",Gannotation,".csv"))

#Anntation
Annotation=rowRanges(rse)
write.csv(Annotation,paste0(Project,"_",Type,"_Annotation",Gannotation,".csv"))


#Raw read cont
ReadCount=assay(rse)



if ( length(SampleNames)==0 ) { 
    write("All samples are reporting..\n", stderr())
    write.csv(ReadCount,paste0(Project,"_",Type,"_RawReadCount",Gannotation,".csv"))
}else {

Snames<-read.table(SampleNames,header=TRUE,sep = ",")
Snames$SampleName<-trimws(Snames$SampleName, which = c("both"))
SelectedCOlumns=as.vector(Snames$SampleName)

ReadCount<-subset(ReadCount, select=SelectedCOlumns)

write.csv(ReadCount,paste0(Project,"_",Type,"_RawReadCount",Gannotation,".csv"))

}



