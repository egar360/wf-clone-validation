import groovy.json.JsonBuilder

// processes required for assembly using canu

process assembleCore_canu {
    errorStrategy = { task.attempt <= 4 ? 'retry' : 'ignore' }
    maxRetries 4
    label "canu"
    cpus params.threads
    memory "7GB"

    input:
        tuple val(meta), path(fastq)

    output:
        tuple val(meta), path("${meta.alias}.reconciled.fasta"), optional: true, emit: assembly
        tuple val(meta), path("${meta.alias}.downsampled.fastq"), optional: true, emit: downsampled
        tuple val(meta.alias), env(STATUS), emit: status

    script:
        String cluster_dir = "trycycler/cluster_001"
        int coverage_target = params.assm_coverage * 3
        int min_dep = (params.assm_coverage / 3) * 2
        int min_len = 100
        int max_len = (meta.approx_size as Integer) * 1.2
        def fast = params.canu_fast == true ? '-fast' : ''
        // WSL does not support named pipes used by Canu, setting these parameters avoids their use
        def windows_params = System.properties['os.version'].toLowerCase().contains("wsl") ? """\
        -mhapPipe=false \
        -purgeOverlaps=false \
        -saveOverlaps=true """ : ""
        def seqkit_threads = params.threads >= 6 ? 2 : 1

    """
    ############################################################
    # Trimming
    ############################################################
    STATUS="Failed to trim reads"
    (
        if [[ $params.trim_length -gt 0 ]]; then
            seqkit subseq -j $seqkit_threads -r $params.trim_length:-$params.trim_length $fastq
        else
            cat $fastq
        fi \\
        | seqkit subseq -j $seqkit_threads -r 1:$max_len \\
        | seqkit seq -j $seqkit_threads -m $min_len -Q $params.min_quality -g > "${meta.alias}.trimmed.fastq"
    ) &&

    ############################################################
    # Downsampling
    ############################################################
    STATUS="Failed to downsample reads" &&
    (rasusa \\
        --coverage $coverage_target \\
        --genome-size ${meta.approx_size} \\
        --input "${meta.alias}.trimmed.fastq" > "${meta.alias}.downsampled.fastq") &&


    ############################################################
    # Subsampling + Canu + Trycycler with internal retries
    ############################################################
    attempt=1
    max_attempts=4
    approx_size=${meta.approx_size}
    final_ok=0

    mkdir -p fallback_candidates
    echo -e "attempt\\tstage\\tmessage" > fallback_candidates/attempt_log.tsv

    while [ \$attempt -le \$max_attempts ]; do
        echo "[attempt \$attempt/\$max_attempts] starting subsample -> canu -> trycycler" 1>&2

        # Clean per-attempt dirs from prior attempt (but keep fallbacks)
        rm -rf sets assm_sample_* trycycler *.deconcat.fasta *.trimmed.fasta 2>/dev/null || true

        ############################################################
        # Subsetting (nondeterministic; this is why retries help)
        ############################################################
        STATUS="Failed to Subset reads (attempt \$attempt)"
        if ! trycycler subsample \\
            --count 3 \\
            --min_read_depth $min_dep \\
            --reads "${meta.alias}.downsampled.fastq" \\
            --out_dir sets \\
            --genome_size ${meta.approx_size}; then
            echo -e "\$attempt\\tsubsample\\tcommand_failed" >> fallback_candidates/attempt_log.tsv
            attempt=\$((attempt+1))
            continue
        fi

        if ! ls sets/sample_*.fastq >/dev/null 2>&1; then
            echo "[attempt \$attempt] subsample produced no read sets" 1>&2
            echo -e "\$attempt\\tsubsample\\tno_sets_created" >> fallback_candidates/attempt_log.tsv
            attempt=\$((attempt+1))
            continue
        fi

        ############################################################
        # Assembly (Canu) on the subsets
        ############################################################
        STATUS="Failed to assemble using Canu (attempt \$attempt)"
        canu_ok=1
        for SUBSET in \$(ls sets/sample_*.fastq); do
            SUBSET_NAME=\$(basename -s .fastq \$SUBSET)
            canu \\
                -p \$SUBSET_NAME \\
                -d assm_\${SUBSET_NAME} \\
                -maxThreads=$task.cpus \\
                genomeSize=${meta.approx_size} \\
                $fast \\
                -nanopore \$SUBSET \\
                $windows_params || canu_ok=0
        done

        if [ \$canu_ok -eq 0 ] && ! ls assm_sample_0*/*.contigs.fasta >/dev/null 2>&1; then
            echo "[attempt \$attempt] canu failed (no contigs)" 1>&2
            echo -e "\$attempt\\tcanu\\tno_contigs" >> fallback_candidates/attempt_log.tsv
            attempt=\$((attempt+1))
            continue
        fi

        ############################################################
        # Trim assemblies + deconcatenate
        ############################################################
        STATUS="Failed to trim/deconcatenate assemblies (attempt \$attempt)"
        for assembly in \$(ls assm_sample_0*/*.contigs.fasta 2>/dev/null); do
            assembly_name=\$(basename -s .fasta \$assembly)
            trim.py "\$assembly" -o "\${assembly_name}.trimmed.fasta" || true
            deconcatenate.py \\
                "\${assembly_name}.trimmed.fasta" \\
                -o "\${assembly_name}.deconcat.fasta" \\
                --approx_size ${meta.approx_size} || true
        done

        if ! ls *.deconcat.fasta >/dev/null 2>&1; then
            echo "[attempt \$attempt] no deconcat contigs produced" 1>&2
            echo -e "\$attempt\\tdeconcat\\tno_deconcat_fastas" >> fallback_candidates/attempt_log.tsv
            attempt=\$((attempt+1))
            continue
        fi

        ############################################################
        # Trycycler cluster
        ############################################################
        STATUS="Failed to cluster contigs (attempt \$attempt)"
        if ! trycycler cluster \\
            --assemblies *.deconcat.fasta \\
            --reads "${meta.alias}.downsampled.fastq" \\
            --out_dir trycycler; then
            echo -e "\$attempt\\tcluster\\tcommand_failed" >> fallback_candidates/attempt_log.tsv
            attempt=\$((attempt+1))
            continue
        fi

        if [ ! -d "trycycler/cluster_001/1_contigs" ]; then
            echo "[attempt \$attempt] trycycler cluster did not produce cluster_001" 1>&2
            echo -e "\$attempt\\tcluster\\tmissing_cluster_001" >> fallback_candidates/attempt_log.tsv
            attempt=\$((attempt+1))
            continue
        fi

        ############################################################
        # Save 1_contigs candidates for this attempt (always)
        ############################################################
        for f in $cluster_dir/1_contigs/*.fasta; do
            [ -e "\$f" ] || continue
            cp -f "\$f" "fallback_candidates/attempt\${attempt}_\$(basename "\$f")" 2>/dev/null || true
        done

        ############################################################
        # Filter contigs BEFORE reconcile (required for stability)
        ############################################################
        STATUS="Failed to filter contigs (attempt \$attempt)"
        if ! python3 ${workflow.projectDir}/bin/filter_trycycler_contigs.py \\
            --approx_size ${meta.approx_size} \\
            --cluster_dir $cluster_dir \\
            --phylip trycycler/contigs.phylip \\
            --max_keep 26; then
            echo "[attempt \$attempt] filter failed" 1>&2
            echo -e "\$attempt\\tfilter\\tcommand_failed" >> fallback_candidates/attempt_log.tsv
            attempt=\$((attempt+1))
            continue
        fi

        ############################################################
        # Trycycler reconcile + msa + partition + consensus
        ############################################################
        STATUS="Failed to reconcile assemblies (attempt \$attempt)"
        if ! trycycler reconcile \\
            --reads "${meta.alias}.downsampled.fastq" \\
            --cluster_dir $cluster_dir \\
            --max_trim_seq_percent 5 \\
            --max_add_seq_percent 10; then
            echo -e "\$attempt\\treconcile\\tcommand_failed" >> fallback_candidates/attempt_log.tsv
            attempt=\$((attempt+1))
            continue
        fi

        if [ ! -s "$cluster_dir/2_all_seqs.fasta" ]; then
            echo "[attempt \$attempt] reconcile did not produce 2_all_seqs.fasta" 1>&2
            echo -e "\$attempt\\treconcile\\tmissing_2_all_seqs" >> fallback_candidates/attempt_log.tsv
            attempt=\$((attempt+1))
            continue
        fi

        # Keep the reconciled set too (useful fallback if later steps fail)
        cp -f "$cluster_dir/2_all_seqs.fasta" "fallback_candidates/attempt\${attempt}_2_all_seqs.fasta" 2>/dev/null || true

        STATUS="Failed to run Trycycler MSA/partition/consensus (attempt \$attempt)"
        if ! trycycler msa --cluster_dir $cluster_dir; then
            echo -e "\$attempt\\tmsa\\tcommand_failed" >> fallback_candidates/attempt_log.tsv
            attempt=\$((attempt+1))
            continue
        fi

        if ! trycycler partition --reads "${meta.alias}.downsampled.fastq" --cluster_dirs $cluster_dir; then
            echo -e "\$attempt\\tpartition\\tcommand_failed" >> fallback_candidates/attempt_log.tsv
            attempt=\$((attempt+1))
            continue
        fi

        if ! trycycler consensus --cluster_dir $cluster_dir; then
            echo -e "\$attempt\\tconsensus\\tcommand_failed" >> fallback_candidates/attempt_log.tsv
            attempt=\$((attempt+1))
            continue
        fi

        if [ -s "$cluster_dir/7_final_consensus.fasta" ]; then
            mv "$cluster_dir/7_final_consensus.fasta" "${meta.alias}.reconciled.fasta"
            STATUS="Completed successfully (attempt \$attempt)"
            final_ok=1
            break
        else
            echo "[attempt \$attempt] consensus not produced" 1>&2
            echo -e "\$attempt\\tconsensus\\tmissing_7_final_consensus" >> fallback_candidates/attempt_log.tsv
            attempt=\$((attempt+1))
            continue
        fi
    done

    ############################################################
    # Final fallback selection (if no attempt succeeded)
    ############################################################
    if [ \$final_ok -ne 1 ]; then
        echo "[final] no consensus after \$max_attempts attempts, selecting best fallback" 1>&2
        STATUS="FALLBACK after \${max_attempts} attempts"

        best=""
        best_diff=999999999

        for f in fallback_candidates/attempt*_*.fasta; do
            [ -e "\$f" ] || continue
            len=\$(seqkit fx2tab -n -l "\$f" 2>/dev/null | awk '{print \$2}' | head -n 1)
            [ -n "\$len" ] || continue
            diff=\$(( len > approx_size ? len-approx_size : approx_size-len ))
            if [ \$diff -lt \$best_diff ]; then
                best_diff=\$diff
                best="\$f"
            fi
        done

        if [ -n "\$best" ]; then
            cp -f "\$best" "${meta.alias}.reconciled.fasta"
            echo "[final] selected \$best (diff=\$best_diff)" 1>&2
        else
            echo "[final] no fallback candidates found" 1>&2
        fi
    fi
    """
}
