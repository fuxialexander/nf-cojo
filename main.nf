#!/usr/bin/env nextflow
nextflow.preview.dsl=2


/*
 * Prepare input files for cojo
 */

// Target format
// SNP A1 A2 freq b se p N 
// rs1001 A G 0.8493 0.0024 0.0055 0.6653 129850 
// rs1002 C G 0.0306 0.0034 0.0115 0.7659 129799 
// rs1003 A C 0.5128 0.0045 0.0038 0.2319 129830

process prepare_cojo_input {
    tag "${name}"
    publishDir "${params.tempdir}"

    input:
        tuple file(range), val(name)
    output:
        tuple val(name), file("${name}.ma")
    script:
    """
    plink --extract range ${range} \
    --bfile ${params.raw_dir}/${name} --make-just-bim --freq --out output
    tr -s ' ' < output.frq \
    | cut -d" " -f1,2 --complement  \
    | tail -n +2 \
    | awk -F" " '{print \$1,\$2,\$3,\$4}' \
    | tr ' ' '\t' \
    | sort -k1,1 > output.part1
    zcat ${params.sum_stat_dir}/${name}.P1.assoc.logistic.gz \
    | tr -s ' ' \
    | cut -d" " -f1,2,4,6,10,11,12 --complement \
    | tr ' ' '\t' \
    | tail -n +2 \
    | sort -k1,1 > output.part2
    join -1 1 -2 1 output.part1 output.part2 \
    | awk -F " " '{OFS="\t"; print \$1,\$2,\$3,\$4,\$7,\$8,\$9,\$6}' > ${name}.ma
    """
}

process get_sig_list {
    tag "${name}"
    publishDir "${params.tempdir}"
    output:
    file("*")
    script:
    """
    cat ${params.tempdir}/test_gwas_*.ma | awk '\$7<5e-8' | cut -f1 | sort -u > significant_snps.snplist
    """
}

process get_sig_names {
    tag "${name}"
    publishDir "${params.tempdir}"
    output:
    file("*")
    script:
    """
    cat ${params.tempdir}/test_gwas_*.ma | awk '\$7<5e-8' | cut -f1 -d: | sort -u | awk '{print "${params.plinkpre}"\$0}' > sig_names
    """
}

process run_cojo {
    tag "${name}"
    cpus 6
    maxForks 1
    publishDir "${params.outdir}"
    input:
        tuple val(name), file(cojo), file(snplist)
    output:
        file("*.cojo")
    script:
    """
    gcta64 --bfile ${params.raw_dir}/${name} --chr ${name.split("_")[2]} --maf 0.01 --cojo-file ${cojo} --cojo-slct --out ${name} --thread-num 6
    """
}

workflow {
    main:
    names = Channel.fromPath( "${params.raw_dir}/${params.plinkpat}.bim" )
                   .map{ file -> file.baseName }
    range = Channel.fromPath( "${params.data_dir}/plink_range.txt" )
                   .combine(names)
    cojo = prepare_cojo_input(range)
    sig_list = get_sig_list()
    get_sig_names().splitCsv().join(cojo).combine(sig_list) | run_cojo
}