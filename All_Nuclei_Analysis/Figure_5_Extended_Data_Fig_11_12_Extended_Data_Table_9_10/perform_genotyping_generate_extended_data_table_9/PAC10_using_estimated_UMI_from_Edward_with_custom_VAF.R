#!/bin/env Rscript

library(yyxMosaicHunter)
library(pryr)
library(stats)
library(data.table)

args<-commandArgs(TRUE)
verbose=TRUE
read_file=args[1]
write_file=args[2]
ignore_single_read_umi=args[3]
vaf=as.numeric(args[4])

MaxBaseQual=93
MaxDepth=1000

set.seed(1)

input=as.data.frame(fread(file=read_file,header=TRUE)) #changed this
input=input[nchar(input$estimated.umi)==12&is.na(input$misc)&input$call_description %in% c("ref","alt"),] #changed this 
input=input[order(input[,1],input[,3]),]

read=data.frame(input$best.bc,input$umi,input$estimated.umi,input$edit.distance,input$call_description,input$call_bp_quality,stringsAsFactors=FALSE)
colnames(read)=c("best.bc","raw.umi","estimated.umi","edit.distance","call.description","call.bp.quality")

umi=numeric(0)
best.bc=""
estimated.umi=""
ref.bases=""
alt.bases=""
for(i in 1:nrow(read))
{
	if(best.bc==read$best.bc[i]&estimated.umi==read$estimated.umi[i])
	{
		if(read$call.description[i]=="ref")
		{
			ref.bases=paste0(ref.bases,rawToChar(as.raw(read$call.bp.quality[i]+33)))
		}
		if(read$call.description[i]=="alt")
		{
			alt.bases=paste0(alt.bases,rawToChar(as.raw(read$call.bp.quality[i]+33)))
		}
	} else
	{
		if(i>1)
		{
			if(nchar(ref.bases)+nchar(alt.bases)>MaxDepth)
			{
				ref.num=round(nchar(ref.bases)/(nchar(ref.bases)+nchar(alt.bases))*MaxDepth)
				alt.num=round(nchar(alt.bases)/(nchar(ref.bases)+nchar(alt.bases))*MaxDepth)
				ref.bases=paste(sample(unlist(strsplit(ref.bases,"")),ref.num),collapse="")
				alt.bases=paste(sample(unlist(strsplit(alt.bases,"")),alt.num),collapse="")
			}
			print(paste0("ref: ",ref.bases))
			print(paste0("alt: ",alt.bases))
			mh_result=yyx_wrapped_mosaic_hunter_for_one_site(ref.bases,alt.bases,output_log10=TRUE,ref_het_alt_mosaic_prior=c(0.5,0,0.5,0))
			if(mh_result$ref_het_alt_mosaic_posterior[1]>mh_result$ref_het_alt_mosaic_posterior[3])
			{
				genotype="ref"
				genoqual=round((mh_result$ref_het_alt_mosaic_posterior[1]-mh_result$ref_het_alt_mosaic_posterior[3])*10)
			} else
			{
				genotype="alt"
				genoqual=round((mh_result$ref_het_alt_mosaic_posterior[3]-mh_result$ref_het_alt_mosaic_posterior[1])*10)
			}
			umi=rbind(umi,c(best.bc,estimated.umi,genotype,genoqual,ref.bases,alt.bases))
		}
		
		best.bc=read$best.bc[i]
		estimated.umi=read$estimated.umi[i]
		ref.bases=""
		alt.bases=""
		
		if(read$call.description[i]=="ref")
		{
			ref.bases=paste0(ref.bases,rawToChar(as.raw(read$call.bp.quality[i]+33)))
		}
		if(read$call.description[i]=="alt")
		{
			alt.bases=paste0(alt.bases,rawToChar(as.raw(read$call.bp.quality[i]+33)))
		}
	}
	print(paste0("read: ",i))
}

if(nchar(ref.bases)+nchar(alt.bases)>MaxDepth)
{
	ref.num=round(nchar(ref.bases)/(nchar(ref.bases)+nchar(alt.bases))*MaxDepth)
	alt.num=round(nchar(alt.bases)/(nchar(ref.bases)+nchar(alt.bases))*MaxDepth)
	ref.bases=paste(sample(unlist(strsplit(ref.bases,"")),ref.num),collapse="")
	alt.bases=paste(sample(unlist(strsplit(alt.bases,"")),alt.num),collapse="")
}
mh_result=yyx_wrapped_mosaic_hunter_for_one_site(ref.bases,alt.bases,output_log10=TRUE,ref_het_alt_mosaic_prior=c(0.5,0,0.5,0))
if(mh_result$ref_het_alt_mosaic_posterior[1]>mh_result$ref_het_alt_mosaic_posterior[3])
{
	genotype="ref"
	genoqual=round((mh_result$ref_het_alt_mosaic_posterior[1]-mh_result$ref_het_alt_mosaic_posterior[3])*10)
} else
{
	genotype="alt"
	genoqual=round((mh_result$ref_het_alt_mosaic_posterior[3]-mh_result$ref_het_alt_mosaic_posterior[1])*10)
}
umi=rbind(umi,c(best.bc,estimated.umi,genotype,genoqual,ref.bases,alt.bases))

umi=as.data.frame(umi,stringsAsFactors=F)
colnames(umi)=c("best.bc","estimated.umi","genotype","genoqual","ref.bases","alt.bases")
umi$genoqual=as.numeric(umi$genoqual)
umi$genoqual[umi$genoqual>MaxBaseQual]=MaxBaseQual
umi=umi[umi$genoqual>0,]
if(ignore_single_read_umi==TRUE)
{
	umi=umi[nchar(umi$ref.bases)+nchar(umi$alt.bases)>1,]
}

bc=numeric(0)
best.bc=""
ref.bases=""
alt.bases=""
for(i in 1:nrow(umi))
{
	if(best.bc==umi$best.bc[i])
	{
		if(umi$genotype[i]=="ref")
		{
			ref.bases=paste0(ref.bases,rawToChar(as.raw(umi$genoqual[i]+33)))
		}
		if(umi$genotype[i]=="alt")
		{
			alt.bases=paste0(alt.bases,rawToChar(as.raw(umi$genoqual[i]+33)))
		}
	} else
	{
		if(i>1)
		{
			if(nchar(ref.bases)+nchar(alt.bases)>MaxDepth)
			{
				ref.num=round(nchar(ref.bases)/(nchar(ref.bases)+nchar(alt.bases))*MaxDepth)
				alt.num=round(nchar(alt.bases)/(nchar(ref.bases)+nchar(alt.bases))*MaxDepth)
				ref.bases=paste(sample(unlist(strsplit(ref.bases,"")),ref.num),collapse="")
				alt.bases=paste(sample(unlist(strsplit(alt.bases,"")),alt.num),collapse="")
			}
			print(paste0("ref: ",ref.bases))
			print(paste0("alt: ",alt.bases))
			mh_result=yyx_wrapped_mosaic_hunter_for_one_site(ref.bases,alt.bases,output_log10=TRUE,ref_het_alt_mosaic_prior=c(1-vaf,vaf,0,0))
			if(mh_result$ref_het_alt_mosaic_posterior[1]>mh_result$ref_het_alt_mosaic_posterior[2])
			{
				genotype="ref-hom"
				genoqual=round((mh_result$ref_het_alt_mosaic_posterior[1]-mh_result$ref_het_alt_mosaic_posterior[2])*10)
			} else
			{
				genotype="hetero"
				genoqual=round((mh_result$ref_het_alt_mosaic_posterior[2]-mh_result$ref_het_alt_mosaic_posterior[1])*10)
			}
			bc=rbind(bc,c(best.bc,genotype,genoqual,ref.bases,alt.bases))
		}
		
		best.bc=umi$best.bc[i]
		estimated.umi=umi$estimated.umi[i]
		ref.bases=""
		alt.bases=""
		
		if(umi$genotype[i]=="ref")
		{
			ref.bases=paste0(ref.bases,rawToChar(as.raw(umi$genoqual[i]+33)))
		}
		if(umi$genotype[i]=="alt")
		{
			alt.bases=paste0(alt.bases,rawToChar(as.raw(umi$genoqual[i]+33)))
		}
	}
	print(paste0("umi: ",i))
}

if(nchar(ref.bases)+nchar(alt.bases)>MaxDepth)
{
	ref.num=round(nchar(ref.bases)/(nchar(ref.bases)+nchar(alt.bases))*MaxDepth)
	alt.num=round(nchar(alt.bases)/(nchar(ref.bases)+nchar(alt.bases))*MaxDepth)
	ref.bases=paste(sample(unlist(strsplit(ref.bases,"")),ref.num),collapse="")
	alt.bases=paste(sample(unlist(strsplit(alt.bases,"")),alt.num),collapse="")
}
mh_result=yyx_wrapped_mosaic_hunter_for_one_site(ref.bases,alt.bases,output_log10=TRUE,ref_het_alt_mosaic_prior=c(1-vaf,vaf,0,0))
if(mh_result$ref_het_alt_mosaic_posterior[1]>mh_result$ref_het_alt_mosaic_posterior[2])
{
	genotype="ref-hom"
	genoqual=round((mh_result$ref_het_alt_mosaic_posterior[1]-mh_result$ref_het_alt_mosaic_posterior[2])*10)
} else
{
	genotype="hetero"
	genoqual=round((mh_result$ref_het_alt_mosaic_posterior[2]-mh_result$ref_het_alt_mosaic_posterior[1])*10)
}
bc=rbind(bc,c(best.bc,genotype,genoqual,ref.bases,alt.bases))

bc=as.data.frame(bc,stringsAsFactors=F)
colnames(bc)=c("best.bc","genotype","genoqual","ref.bases","alt.bases")
bc$genoqual=as.numeric(bc$genoqual)
bc$genoqual[bc$genoqual>MaxBaseQual]=MaxBaseQual

write.table(bc,file=write_file,sep="\t",row.names=FALSE,quote=FALSE,col.names=TRUE)
