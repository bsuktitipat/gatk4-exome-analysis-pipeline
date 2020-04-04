version 1.0

# Author: Bhoom Suktitipat, Kunaphas Kongkittimanon
# script version: 1.1

# Input

# WORKFLOW DEFINITION
This WDL pipeline implements Variant Filtration with hard filtering approach according to filtering recommendations  on GATK V 3.5.x.
Ref 1) VariantFiltration Criteria: https://gatk.broadinstitute.org/hc/en-us/articles/360037499012 (Update on 2020-03-25)
Ref 2) Imporvoing the filter: https://gatk.broadinstitute.org/hc/en-us/articles/360035890471 (Update on 2020-03-31)
Ref 3) VariantFiltration Doc: https://gatk.broadinstitute.org/hc/en-us/articles/360037226192
QualByDepth (QD)
FisherStrand (FS)
StrandOddsRatio (SOR)
RMSMappingQuality (MQ)
MappingQualityRankSumTest (MQRankSum)
ReadPosRankSumTest (ReadPosRankSum)
SNP :  'QD < 2.0 || MQ < 40.0 || FS > 60.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0' 
INDEL : 'QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0 || SOR > 10.0'
INDEL (10 or more samples): 'QD < 2.0 || FS > 200.0 || SOR > 10.0 || ReadPosRankSum < -20.0 || InbreedingCoeff < -0.8'
For GATK 3.6.x - 4.x
 --filterExpression "QD < 2.0 || FS > 30.0 || SOR > 3.0 || MQ < 40.0 || MQRankSum < -3.0 || ReadPosRankSum < -3.0"
However, We will use our old traditional way.

# Note:
This pipeline uses traditional method from GATK 3.5.X,CombineGVCFs => GenotypeGVCFs, and no batch_size scatter applied yet.
For CombineVariants ,we use GATK 3.7 


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
	Boolean gt_10samples = false

	# HardFilter Condition
	String SNPfilterExpression =  "QD < 2.0 || MQ < 40.0 || FS > 60.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" 
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
			vcf_file = CombineGVCFs.output_cohort_g_vcf,
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
			type = "SNP",
			rawVCF = GenotypeGVCFs.output_RAW_vcf,
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
			rawIndels = selectIndels.rawSubset,
			filterExpression = INDELfilterExpression ,
			RefFasta = RefFasta,
			RefIndex = RefIndex,
			RefDict = RefDict,
			disk_size = disk_size,
			docker = gatk_docker,
			gatk_path = gatk_path,
			preemptible = preemptible
	}

	call CombineSnvIndel {
		input:
			filteredSNPs = HardFilterSNP.filteredSNPs,
			filteredIndels = HardFilterIndel.filteredIndels,
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

    docker: "broadinstitute/gatk:4.0.10.1"
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
    	docker: "broadinstitute/gatk:4.0.10.1"
        memory: "4 GB"
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

task CombineSnvIndel {

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
    	docker: "broadinstitute/gatk:4.0.10.1"
	    memory: "8 GB"
			disks: "local-disk " + disk_size + " HDD"
	 		preemptible: preemptible
    }
}