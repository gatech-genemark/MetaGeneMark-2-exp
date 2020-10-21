function get_5p_correct() {
    pf_reference="$1"
    pf_prediction="$2"
    pf_output="$3"

    compp -a $pf_reference -b $pf_prediction -q -n -i | grep -v "#" > $pf_output
}

function get_5p_incorrect() {
    pf_reference="$1"
    pf_prediction="$2"
    pf_output="$3"

    compp -a $pf_reference -b $pf_prediction -q -n -s -l | grep -v "#" > $pf_output
}

function create_3p_key() {
    local accession="$1"
    local left="$2"
    local right="$3"
    local strand="$4"

    if [[ $strand == "+" ]]; then
        echo "$accession;$strand;$right"
    else
        echo "$accession;$strand;$left"
    fi
}

function extract_motif_info_from_mgm2_file() {
    awk '{if (NF) print}' "$1" | while read -r line; do
        accession=$(echo "$line" | awk -F "\t" '{print $1}');
        left=$(echo "$line" | awk -F "\t" '{print $4}');
        right=$(echo "$line" | awk -F "\t" '{print $5}');
        strand=$(echo "$line" | awk -F "\t" '{print $7}');

        motif=$(echo "$line" | sed -E 's/^.+start_stop_score[^;]+;\s+([ACGT]+).+$/\1/g')
        pos=$(echo "$line" | sed -E 's/^.+start_stop_score[^;]+;\s+[ACGT]+\s+([0-9]+)\s.+$/\1/g')
        type=$(echo "$line" | sed -E 's/^.+start_stop_score[^;]+;\s+[ACGT]+\s+[0-9]++\s+([0-9]++)\s.+$/\1/g')

        key=$(create_3p_key $accession $left $right $strand)

        echo "$key,$motif,$pos,$type" 
    done | sort -t, -k1
}


function extract_motif_info_from_mprodigal_file() {
    awk '{if (NF) print}' "$1" | while read -r line; do
        accession=$(echo "$line" | awk -F "\t" '{print $1}');
        left=$(echo "$line" | awk -F "\t" '{print $4}');
        right=$(echo "$line" | awk -F "\t" '{print $5}');
        strand=$(echo "$line" | awk -F "\t" '{print $7}');

        motif=$(echo "$line" | sed -E 's/^.+rbs_motif=([^;]+).+$/\1/g')
        pos=$(echo "$line" | sed -E 's/^.+spacer=([^;]+).+$/\1/g')
        type="0"

        key=$(create_3p_key $accession $left $right $strand)

        echo "$key,$motif,$pos,$type" 
    done | sort -t, -k1
}

function extract_motif_info_for_same_gene() {
    local pf_mgm2_view="$1"
    local pf_mprodigal_view="$2"
    local gcfid="$3"
    local gc="$4"
    local tag="$5"

    out_mgm2=$(extract_motif_info_from_mgm2_file $pf_mgm2_view)
    out_mprodigal=$(extract_motif_info_from_mprodigal_file $pf_mprodigal_view)

    paste <(echo "$out_mgm2") <(echo "$out_mprodigal") --delimiters ',' | while read -r line; do
        echo "$gcfid,$gc,$tag,$line"
    done



}

function analyze_mgm2_vs_mprodigal() {
    
    gcfid="$1"

    pf_verified=$data/$gcfid/verified.gff
    pf_mgm2=$runs/$gcfid/mgms_rerun/prediction.gff
    pf_mprodigal=$runs/$gcfid/mprodigal/prediction.gff

    # get correct false predictions by each 
    pf_mgm2_correct="mgm2_correct.gff"
    pf_mgm2_incorrect="mgm2_incorrect.gff"

    pf_mprodigal_correct="mprodigal_correct.gff"
    pf_mprodigal_incorrect="mprodigal_incorrect.gff"

    get_5p_correct $pf_verified $pf_mgm2 $pf_mgm2_correct
    get_5p_incorrect $pf_verified $pf_mgm2 $pf_mgm2_incorrect

    get_5p_correct $pf_verified $pf_mprodigal $pf_mprodigal_correct
    get_5p_incorrect $pf_verified $pf_mprodigal $pf_mprodigal_incorrect


    pf_mgm2_correct_mprodigal_incorrect_view_mgm2="mgm2_correct_mprodigal_incorrect_view_mgm2.gff"
    pf_mgm2_correct_mprodigal_incorrect_view_mprodigal="mgm2_correct_mprodigal_incorrect_view_mprodigal.gff"

    compp -a $pf_mgm2_correct -b $pf_mprodigal_incorrect -q -n -S -L | grep -v "#" > $pf_mgm2_correct_mprodigal_incorrect_view_mgm2
    compp -a $pf_mgm2_correct -b $pf_mprodigal_incorrect -q -n -l -s | grep -v "#" > $pf_mgm2_correct_mprodigal_incorrect_view_mprodigal


    pf_mprodigal_correct_mgm2_incorrect_view_mprodigal="mprodigal_correct_mgm2_incorrect_view_mprodigal.gff"
    pf_mprodigal_correct_mgm2_incorrect_view_mgm2="mprodigal_correct_mgm2_incorrect_view_mgm2.gff"

    compp -a $pf_mprodigal_correct -b $pf_mgm2_incorrect -q -n -S -L | grep -v "#" > $pf_mprodigal_correct_mgm2_incorrect_view_mprodigal
    compp -a $pf_mprodigal_correct -b $pf_mgm2_incorrect -q -n -l -s | grep -v "#" > $pf_mprodigal_correct_mgm2_incorrect_view_mgm2

    gc=$(probuild --stat --gc --seq $data/$gcfid/sequence.fasta | awk '{print $3}')

    echo "Genome,GC,Tag,AKey,AMotif,ASpacer,AMotifType,BKey,BMotif,BSpacer,BMotifType"
    extract_motif_info_for_same_gene $pf_mgm2_correct_mprodigal_incorrect_view_mgm2 $pf_mgm2_correct_mprodigal_incorrect_view_mprodigal $gcfid $gc "MGMS"
    extract_motif_info_for_same_gene $pf_mprodigal_correct_mgm2_incorrect_view_mgm2 $pf_mprodigal_correct_mgm2_incorrect_view_mprodigal $gcfid $gc "MetaProdigal"
}



 analyze_mgm2_vs_mprodigal Mycobacterium_tuberculosis_H37Rv_uid57777 > results.txt
 analyze_mgm2_vs_mprodigal Deinococcus_deserti_VCD115_uid58615 | grep -v "BKey" >> results.txt
 analyze_mgm2_vs_mprodigal Halobacterium_salinarum_R1_uid61571 | grep -v "BKey" >> results.txt
