## Creator: KK (Kunaphas Kongkititmanon)
## script version: 1
##
## runtime: GCP+docker+GCS
##


## Input


## WORKFLOW DEFINITION
# This WDL pipeline implements Variant Filtration with hard filtering approach according to filtering recommendations  on GATK V 3.5.x.
# Reference: https://gatk.broadinstitute.org/hc/en-us/articles/360037499012
# SNP (whole genome):  'QD < 2.0 || MQ < 30.0 || FS > 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || DP < 8.0'
# SNP (exomes) :  'QD < 2.0 || MQ < 30.0 || FS > 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 '
# INDEL : 'QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0'
# INDEL (10 or more samples): 'QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0 || InbreedingCoeff < -0.8'
# For GATK 3.6.x - 4.x
#  --filterExpression "QD < 2.0 || FS > 30.0 || SOR > 3.0 || MQ < 40.0 || MQRankSum < -3.0 || ReadPosRankSum < -3.0"
# However, We will use our old traditional way.

## Note:
# This pipeline uses traditional method from GATK 3.5.X,CombineGVCFs => GenotypeGVCFs, and no batch_size scatter applied yet.
# For CombineVariants ,we use GATK 3.7 
#

workflow VariantFiltration_HardFilter {
	Int disk_size = 100
	Int preemptible = 3
	String gatk_docker
	String gatk_path

	File RefFastaref_fasta
	File ref_fasta_index
	File ref_dict

	File dbsnp_vcf
	File dbsnp_vcf_index

	Array[File] input_gvcfs
	Array[File] input_gvcfs_indices

	# whole genome by default
	Boolean exomes = false
	Boolean gt_10samples = false

	# HardFilter Condition
	String SNPfilterExpression =  if (exomes) then   "QD < 2.0 || MQ < 30.0 || FS > 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" else  "QD < 2.0 || MQ < 30.0 || FS > 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || DP < 8.0"
	String INDELfilterExpression = if (gt_10samples) then " QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0 || InbreedingCoeff < -0.8 " else  " QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0 "

	call CombineGVCFs  {
		input:
			input_gvcfs = input_gvcfs,
			input_gvcfs_indices = input_gvcfs_indices,
			RefFasta = RefFasta,
			RefIndex = RefIndex,
			RefDict = RefDict,
			disk_size = disk_size,
			docker = gatk_docker,
			gatk_path = gatk_path,
			preemptible = preemptible

	}

	call GenotypeGVCFs {
		input:
			vcf_file =CombineGVCFs.output_cohort_g_vcf,
			RefFasta = RefFasta,
			RefIndex = RefIndex,
			RefDict = RefDict,
			disk_size = disk_size,
			docker = gatk_docker,
			gatk_path = gatk_path,
			preemptible = preemptible
	}

	call SelectVariants as selectSNPs {
		input:
			type="SNP",
			rawVCF=GenotypeGVCFs.output_RAW_vcf,
			RefFasta = RefFasta,
			RefIndex = RefIndex,
			RefDict = RefDict,
			disk_size = disk_size,
			docker = gatk_docker,
			gatk_path = gatk_path,
			preemptible = preemptible
	}

	call SelectVariants as selectIndels {
		input:
			type="INDEL",
			rawVCF=GenotypeGVCFs.output_RAW_vcf,
			RefFasta = RefFasta,
			RefIndex = RefIndex,
			RefDict = RefDict,
			disk_size = disk_size,
			docker = gatk_docker,
			gatk_path = gatk_path,
			preemptible = preemptible
	}

	call HardFilterSNP {
		input:
			rawSNPs=selectSNPs.rawSubset,
			filterExpression =	SNPfilterExpression ,
			RefFasta = RefFasta,
			RefIndex = RefIndex,
			RefDict = RefDict,
			disk_size = disk_size,
			docker = gatk_docker,
			gatk_path = gatk_path,
			preemptible = preemptible
	}

	call HardFilterIndel {
		input:
			rawIndels=selectIndels.rawSubset,
			filterExpression =	INDELfilterExpression ,
			RefFasta = RefFasta,
			RefIndex = RefIndex,
			RefDict = RefDict,
			disk_size = disk_size,
			docker = gatk_docker,
			gatk_path = gatk_path,
			preemptible = preemptible
	}

	call Combine {
		input:
			filteredSNPs=HardFilterSNP.filteredSNPs,
			filteredIndels=HardFilterIndel.filteredIndels,
			RefFasta = RefFasta,
			RefIndex = RefIndex,
			RefDict = RefDict,
			disk_size = disk_size,
			preemptible = preemptible
	}

	output {
		# outputs from the small callset path through the wdl
		Combine.filteredVCF
		GenotypeGVCFs.output_RAW_vcf
	}
}

task CombineGVCFs {
	File RefFasta
	File RefIndex
	File RefDict

	Array[File] input_gvcfs
	Array[File] input_gvcfs_indices

	String gatk_path
	String docker
	Int disk_size
	Int preemptible

	command {
		${gatk_path} \
			CombineGVCFs  \
			-R ${RefFasta} \
			-V ${sep=' -V ' input_gvcfs} \
			-O cohort.g.vcf
	}

	runtime {
		docker: docker
		memory: "7 GB"
		cpu: "2"
		disks: "local-disk " + disk_size + " HDD"
		preemptible: preemptible
	}

	output {
		File output_cohort_g_vcf= "cohort.g.vcf"
	}

}


task GenotypeGVCFs {

  File vcf_file

  File RefFasta
  File RefIndex
  File RefDict

  String gatk_path
  String docker
  Int disk_size
  Int preemptible

  command {
    ${gatk_path} --java-options "-Xmx5g -Xms5g" \
      GenotypeGVCFs \
     -R ${RefFasta} \
		 -V ${vcf_file} \
     -O cohort.RAW.vcf \
     -G StandardAnnotation \
  }
  runtime {

    docker: docker
    memory: "7 GB"
    cpu: "2"
    disks: "local-disk " + disk_size + " HDD"
    preemptible: preemptible

  }
  output {
    File output_RAW_vcf = "cohort.RAW.vcf"
  }
}

task SelectVariants {

	File RefFasta
	File RefIndex
	File RefDict

	String type
	File rawVCF

		String gatk_path
		String docker
		Int disk_size
		Int preemptible
	command {
	 ${gatk_path} \
			 SelectVariants \
			-R ${RefFasta} \
			-V ${rawVCF} \
			 --select-type ${type} \
			-O vcf_raw.${type}.vcf
	}
	output {
		File rawSubset = "vcf_raw.${type}.vcf"
	}
    runtime {
    	docker: docker
        memory: "2 GB"
        cpu: "1"
		disks: "local-disk " + disk_size + " HDD"
		preemptible: preemptible
    }
}

task HardFilterSNP {
	File RefFasta
	File RefIndex
	File RefDict

	File rawSNPs
  String filterExpression
	String gatk_path
	String docker
	Int disk_size
	Int preemptible
	command {
		 ${gatk_path} \
			VariantFiltration \
			-R ${RefFasta} \
			-V ${rawSNPs} \
			--filter-expression " ${filterExpression} " \
			--filter-name "snp_filter" \
			-O vcf_filtered.snps.vcf
	}
	output {
		File filteredSNPs = "vcf_filtered.snps.vcf"
	}
    runtime {
    	docker: docker
        memory: "2 GB"
        cpu: "1"
		disks: "local-disk " + disk_size + " HDD"
		preemptible: preemptible
    }
}

task HardFilterIndel {
	File RefFasta
	File RefIndex
	File RefDict

	File rawIndels
	String filterExpression
	String gatk_path
	String docker
	Int disk_size
	Int preemptible
	command {
		 ${gatk_path} \
			 VariantFiltration \
			-R ${RefFasta} \
			-V ${rawIndels} \
			--filter-expression " ${filterExpression} " \
			--filter-name "indel_filter" \
			-O vcf_filtered.indels.vcf
	}
	output {
		File filteredIndels = "vcf_filtered.indels.vcf"
	}
    runtime {
    	 docker: docker
		 memory: "2 GB"
		 cpu: "1"
	     disks: "local-disk " + disk_size + " HDD"
		 preemptible: preemptible
	 }
}

task Combine {

	File RefFasta
	File RefIndex
	File RefDict

	File filteredSNPs
	File filteredIndels


	Int disk_size
	Int preemptible

	command {
		java -jar /usr/GenomeAnalysisTK.jar \
			-T CombineVariants \
			-R ${RefFasta} \
			-V ${filteredSNPs} \
			-V ${filteredIndels} \
			--genotypemergeoption UNSORTED \
			-o vcf_filtered.snps.indels.vcf
	}
	output {
		File filteredVCF = "vcf_filtered.snps.indels.vcf"
	}
    runtime {
    	docker: "broadinstitute/gatk3:3.7-0"
	    memory: "2 GB"
			disks: "local-disk " + disk_size + " HDD"
	 		preemptible: preemptible
    }
}