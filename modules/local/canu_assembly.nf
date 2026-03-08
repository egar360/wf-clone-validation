import groovy.json.JsonBuilder

// processes required for assembly using canu


process assembleCore_canu {
    errorStrategy = {task.attempt <= 4 ? 'retry' : 'ignore'}
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
        int exit_number = task.attempt <= 4 ? 1 : 0
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
        fi \
        | seqkit subseq -j $seqkit_threads -r 1:$max_len \
        | seqkit seq -j $seqkit_threads -m $min_len -Q $params.min_quality -g > "${meta.alias}.trimmed.fastq"
    ) &&

    ############################################################
    # Downsampling
    ############################################################

    STATUS="Failed to downsample reads" &&
    (rasusa \
        --coverage $coverage_target \
        --genome-size ${meta.approx_size} \
        --input "${meta.alias}.trimmed.fastq" > "${meta.alias}.downsampled.fastq") &&
        

    ############################################################
    # Subsetting
    ############################################################
    STATUS="Failed to Subset reads" &&
    (trycycler subsample \
        --count 3 \
        --min_read_depth $min_dep \
        --reads "${meta.alias}.downsampled.fastq" \
        --out_dir sets \
        --genome_size ${meta.approx_size}) && 

    ############################################################
    # Assembly
    ############################################################
    STATUS="Failed to assemble using Canu" &&
    (for SUBSET in \$(ls sets/sample_*.fastq)
    do
        SUBSET_NAME=\$(basename -s .fastq \$SUBSET)
        canu \
            -p \$SUBSET_NAME \
            -d assm_\${SUBSET_NAME} \
            -maxThreads=$task.cpus \
            genomeSize=${meta.approx_size} \
            $fast \
            -nanopore \$SUBSET \
            $windows_params
    done) && 

    ############################################################
    # Trim assemblies
    ############################################################
    STATUS="Failed to trim Assembly" &&
    (for assembly in \$(ls assm_sample_0*/*.contigs.fasta)
    do  
        echo \$assembly
        assembly_name=\$(basename -s .fasta \$assembly)
        trim.py \
            \$assembly \
            -o \${assembly_name}.trimmed.fasta
        deconcatenate.py \
            \${assembly_name}.trimmed.fasta \
            -o \${assembly_name}.deconcat.fasta \
            --approx_size ${meta.approx_size}
    done
    ls *.deconcat.fasta > /dev/null 2>&1) &&

    ############################################################
    # Reconciliation
    ############################################################
    STATUS="Failed to reconcile assemblies" &&
    (trycycler cluster \
        --assemblies *.deconcat.fasta \
        --reads "${meta.alias}.downsampled.fastq" \
        --out_dir trycycler) &&
    
    ############################################################
    # Filter Trycycler contigs (fix >26 contigs + length outliers)
    ############################################################
    STATUS="Failed to filter contigs for Trycycler reconcile" &&
    cluster_dir="trycycler/cluster_001" &&
    approx_size="${meta.approx_size}" &&
    
    python3 - <<PY
    import re, statistics
    from pathlib import Path
    
    approx_size = int("${meta.approx_size}")
    cluster_dir = Path("${cluster_dir}")
    contig_dir = cluster_dir / "1_contigs"
    phylip = Path("trycycler/contigs.phylip")
    max_keep = 26
    
    reject_dir = cluster_dir / "1_contigs_rejected"
    reject_dir.mkdir(parents=True, exist_ok=True)
    
    meta = {}
    rx = re.compile(r'_(?P<tig>[A-Z]_tig\\d+?)_(?P<bp>\\d+)_bp_(?P<depth>[0-9.]+)x$')
    
    lines = phylip.read_text().splitlines()
    for line in lines[1:]:
        parts = line.split()
        if not parts:
            continue
        nm = parts[0]
        m = rx.search(nm)
        if not m:
            continue
        tig = m.group("tig")
        bp = int(m.group("bp"))
        depth = float(m.group("depth"))
        if tig not in meta or depth > meta[tig][1]:
            meta[tig] = (bp, depth)
    
    fasta_files = sorted(contig_dir.glob("*.fasta"))
    if not fasta_files:
        raise SystemExit("No contigs found in trycycler/cluster_001/1_contigs")
    
    items = []
    missing = []
    for fp in fasta_files:
        tig = fp.stem
        if tig in meta:
            bp, depth = meta[tig]
            items.append((tig, fp, bp, depth))
        else:
            missing.append(tig)
            items.append((tig, fp, 0, 0.0))
    
    bps_nonzero = [bp for _, _, bp, _ in items if bp > 0]
    if not bps_nonzero:
        center = approx_size
    else:
        med = statistics.median(bps_nonzero)
        center = approx_size if abs(med - approx_size) / max(approx_size, 1) > 0.10 else int(med)
    
    def filter_by_window(pct):
        lo = int(center * (1 - pct))
        hi = int(center * (1 + pct))
        kept = [it for it in items if lo <= it[2] <= hi]
        return kept, lo, hi
    
    pct_list = [0.08, 0.12, 0.20]
    kept = items
    used_pct = None
    lo = hi = None
    for pct in pct_list:
        kept2, lo, hi = filter_by_window(pct)
        if len(kept2) >= 3:
            kept = kept2
            used_pct = pct
            break
    if used_pct is None:
        kept = items
        used_pct = -1
        lo = hi = None
    
    def rank_key(it):
        tig, fp, bp, depth = it
        return (abs(bp - approx_size) if bp > 0 else 10**12, -depth, -bp)
    
    kept_sorted = sorted(kept, key=rank_key)
    kept_final = kept_sorted[:max_keep]
    keep_names = set(tig for tig, *_ in kept_final)
    
    report = []
    for tig, fp, bp, depth in sorted(items, key=lambda x: x[0]):
        keep = tig in keep_names
        if not keep:
            reason = f"REJECT:length_outlier (window={lo}-{hi})" if (used_pct != -1 and (bp < lo or bp > hi)) else "REJECT:ranked_out (cap_26)"
            fp.rename(reject_dir / fp.name)
        else:
            reason = "KEEP"
        report.append((tig, bp, depth, "KEEP" if keep else "REJECT", reason))
    
    out = cluster_dir / "filter_report.tsv"
    with out.open("w") as f:
        f.write("tig\\tbp\\tdepth_x\\tdecision\\treason\\n")
        for tig, bp, depth, decision, reason in report:
            f.write(f"{tig}\\t{bp}\\t{depth:.3f}\\t{decision}\\t{reason}\\n")
    
    print(f"[filter] approx_size={approx_size} center={center} pct={used_pct} kept={len(keep_names)} total={len(items)}")
    if missing:
        print("[filter] WARNING missing in contigs.phylip:", ", ".join(missing[:10]) + (" ..." if len(missing)>10 else ""))
    PY
    &&
    
    ############################################################
    # Reconcile, MSA, partition, consensus
    ############################################################
    STATUS="Failed to reconcile assemblies" &&
    (trycycler reconcile \
        --reads "${meta.alias}.downsampled.fastq" \
        --cluster_dir "${cluster_dir}" \
        --max_trim_seq_percent 5 \
        --max_add_seq_percent 10) &&
    test -s "${cluster_dir}/2_all_seqs.fasta" &&
    (trycycler msa --cluster_dir "${cluster_dir}") &&
    (trycycler partition --reads "${meta.alias}.downsampled.fastq" --cluster_dirs "${cluster_dir}") &&
    (trycycler consensus --cluster_dir "${cluster_dir}")
    ############################################################
    # Exit handling
    ############################################################

    if [ ! -f "${cluster_dir}/7_final_consensus.fasta" ]; then
        if ls ${cluster_dir}/1_contigs/*.fasta 1> /dev/null 2>&1; then
            STATUS="Completed but failed to reconcile"
            (seqkit sort ${cluster_dir}/1_contigs/*.fasta --by-length \
                | seqkit head -n 1 > "${meta.alias}.reconciled.fasta") \
                && echo "Trycycler failed, outputting un-reconciled assembly"
        elif [ "$exit_number" == "1" ]; then
            echo \$STATUS
            echo "Assembly failed, retrying process"
            exit 1
        elif [ "$exit_number" == "0" ]; then
            echo \$STATUS
            echo "Failed final attempt"
        fi
    else
        mv "${cluster_dir}/7_final_consensus.fasta" "${meta.alias}.reconciled.fasta"
        STATUS="Completed successfully"
    fi
    """
}
