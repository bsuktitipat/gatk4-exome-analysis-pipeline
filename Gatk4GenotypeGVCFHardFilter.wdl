version 1.0

# **Author:** Bhoom Suktitipat, Kunaphas Kongkittimanon
# **script version:** 1.1
# **Github**: https://github.com/hypotheses/gatk4-exome-analysis-pipeline/blob/master/Gatk4GenotypeGVCFHardFilter.wdl
# # Input

# # WORKFLOW DEFINITION
# This WDL pipeline implements Variant Filtration with hard filtering approach according to filtering recommendations on GATK V 4

# - QualByDepth (QD)
# - FisherStrand (FS)
# - StrandOddsRatio (SOR)
# - RMSMappingQuality (MQ)
# - MappingQualityRankSumTest (MQRankSum)
# - ReadPosRankSumTest (ReadPosRankSum)

# SNP :  'QD < 2.0 || MQ < 40.0 || FS > 60.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0' 
# INDEL : 'QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0 || SOR > 10.0'
# INDEL (10 or more samples): 'QD < 2.0 || FS > 200.0 || SOR > 10.0 || ReadPosRankSum < -20.0 || InbreedingCoeff < -0.8'

# For GATK 3.6.x - 4.x
#  `--filterExpression "QD < 2.0 || FS > 30.0 || SOR > 3.0 || MQ < 40.0 || MQRankSum < -3.0 || ReadPosRankSum < -3.0"`
# However, We will use our old traditional way.

# # Reference
# 1. VariantFiltration Criteria: https://gatk.broadinstitute.org/hc/en-us/articles/360037499012 (Update on 2020-03-25)
# 2. Imporving the filter: https://gatk.broadinstitute.org/hc/en-us/articles/360035890471 (Update on 2020-03-31)
# 3. VariantFiltration Doc: https://gatk.broadinstitute.org/hc/en-us/articles/360037226192

# # Note:
# This pipeline uses traditional method from GATK 3.5.X,CombineGVCFs => GenotypeGVCFs, and no batch_size scatter applied yet.
# For CombineVariants ,we use GATK 3.7 


workflow VariantFiltration_HardFilter {
	input {
		Int disk_size = 100
		Int preemptible = 3
		String gatk_docker = select_first([gatk_docker,"us.gcr.io/broad-gatk/gatk:4.0.10.1"])
		String gatk_path = select_first([gatk_path,"gatk"])
        String output_prefix = select_first([output_prefix, ""])

		File RefFasta
		File RefIndex
		File RefDict

		Array[File] input_gvcfs
		Array[File] input_gvcfs_indices

		# whole genome by default
		Boolean gt_10samples = false

		# HardFilter Condition
		String SNPfilterExpression =  "QD < 2.0 || MQ < 40.0 || FS > 60.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" 
		String INDELfilterExpression = if (gt_10samples) then " QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0 || InbreedingCoeff < -0.8 " else  " QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0 "
	}
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
			prefix = output_prefix,
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
			rawSNPs = selectSNPs.rawSubset,
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
			prefix = output_prefix,
            filteredSNPs = HardFilterSNP.filteredSNPs,
			filteredIndels = HardFilterIndel.filteredIndels,
			RefFasta = RefFasta,
			RefIndex = RefIndex,
			RefDict = RefDict,
			disk_size = disk_size,
			docker = gatk_docker,
			preemptible = preemptible
	}

	output {
		# outputs from the small callset path through the wdl
		File Combine_filteredVCF = CombineSnvIndel.filteredVCF
		File GenotypeGVCFs_output_RAW_vcf = GenotypeGVCFs.output_RAW_vcf
	}
}

task CombineGVCFs {
	
	input {
		File RefFasta
		File RefIndex
		File RefDict

		Array[File] input_gvcfs
		Array[File] input_gvcfs_indices

		String gatk_path
		String docker
		Int disk_size
		Int preemptible
	}
	
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
		File output_cohort_g_vcf = "cohort.g.vcf"
	}

}


task GenotypeGVCFs {

	input {
		File vcf_file
		File RefFasta
		File RefIndex
		File RefDict
        String prefix
		String gatk_path
		String docker
		Int disk_size
		Int preemptible
	}

	command {
			${gatk_path} --java-options "-Xmx5g -Xms5g" \
			GenotypeGVCFs \
			-R ${RefFasta} \
			-V ${vcf_file} \
			-O ${prefix}.RAW.vcf \
			-G StandardAnnotation
	}
	runtime {
		docker: docker
		memory: "7 GB"
		cpu: "2"
		disks: "local-disk " + disk_size + " HDD"
		preemptible: preemptible
	}
	
	output {
		File output_RAW_vcf = "${prefix}.RAW.vcf"
	}
}

task SelectVariants {
	input {
		File RefFasta
		File RefIndex
		File RefDict

		String type
		File rawVCF

		String gatk_path
		String docker
		Int disk_size
		Int preemptible
	}

	command {
			${gatk_path} SelectVariants \
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
	input {
		File RefFasta
		File RefIndex
		File RefDict
		File rawSNPs
		String filterExpression
		String gatk_path
		String docker
		Int disk_size
		Int preemptible
	}

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
        memory: "4 GB"
        cpu: "1"
		disks: "local-disk " + disk_size + " HDD"
		preemptible: preemptible
    }
}

task HardFilterIndel {
	input {
		File RefFasta
		File RefIndex
		File RefDict

		File rawIndels
		String filterExpression
		String gatk_path
		String docker
		Int disk_size
		Int preemptible
	}

	command {
			${gatk_path} \
			VariantFiltration \
			-R ${RefFasta} \
			-V ${rawIndels} \
			--filter-expression "${filterExpression}" \
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
	input {
		File RefFasta
		File RefIndex
		File RefDict
		File filteredSNPs
		File filteredIndels
		String prefix
        String docker

		Int disk_size
		Int preemptible
	}
	
	command {
			java -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xms8000m \
			-jar /usr/gitc/GATK35.jar \
			-T CombineVariants \
			-R ${RefFasta} \
			-V ${filteredSNPs} \
			-V ${filteredIndels} \
			--genotypemergeoption UNSORTED \
			-o ${prefix}.filtered.snps.indels.vcf
	}
	output {
		File filteredVCF = "${prefix}.filtered.snps.indels.vcf"
	}
    runtime {
    	docker: "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.4.3-1564508330"
    	memory: "8 GB"
    	disks: "local-disk " + disk_size + " HDD"
    	preemptible: preemptible
    }
}