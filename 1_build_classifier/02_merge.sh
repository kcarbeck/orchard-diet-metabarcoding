# import all COI database sources into CRABS format, then merge them into a single database
# katherine carbeck
# 09 sep 2025

# The taxonomic lineage CRABS creates is based on the NCBI taxonomy and CRABS requires the three files downloaded using the --download-taxonomy function, i.e., --names, --nodes, and --acc2tax. From version v 1.0.0, CRABS is capable of resolving synonym and unaccepted names to incorporate a larger number of sequences and diversity in the local reference database. The taxonomic ranks to be included in the taxonomic lineage can be specified using the --ranks parameters. While any taxonomic rank can be included, we recommend using the following input to include all necessary information for most taxonomic classifiers --ranks 'superkingdom;phylum;class;order;family;genus;species'

##############################################################################
#*                 STEP 1 - import sequences
##############################################################################
export PATH="/lustre2/home/lc736_0001/diet/songbird_coi_database/reference_database_creator:$PATH"

# import MIDORI sequences
crabs --import \
  --import-format midori \
  --input raw_downloads/cytb_total_267.fasta \
  --names taxonomy/names.dmp \
  --nodes taxonomy/nodes.dmp \
  --acc2tax taxonomy/nucl_gb.accession2taxid \
  --output processed/midori_total_imported.txt \
  --ranks 'superkingdom;phylum;class;order;family;genus;species'
#  Imported 3027227 out of 3044836 sequences into CRABS format (99.42%)

# import NCBI sequences
crabs --import \
  --import-format ncbi \
  --input raw_downloads/ncbi_arthropoda_coi.fasta \
  --names taxonomy/names.dmp \
  --nodes taxonomy/nodes.dmp \
  --acc2tax taxonomy/nucl_gb.accession2taxid \
  --output processed/ncbi_imported.txt \
  --ranks 'superkingdom;phylum;class;order;family;genus;species'
#  Imported 3029239 out of 3029239 sequences into CRABS format (100.0%)

##############################################################################
#*                 STEP 1A.import bold in parallel
##############################################################################
# function to process one order
process_order() {
    local order=$1
    local output_file="processed/bold_${order}_imported.txt"
    
    #check if output file already exists AND that its not empty
    if [[ -s "$output_file" ]]; then
        echo "Skipping $order - output file already exists: $output_file"
        return 0
    fi
    
    echo "Processing $order..."
    
    #create output dir if it doesnt exist
    mkdir -p "processed"
    
    #run the crabs import command
    if crabs --import \
        --import-format bold \
        --input "raw_downloads/bold_${order}.fasta" \
        --names taxonomy/names.dmp \
        --nodes taxonomy/nodes.dmp \
        --acc2tax taxonomy/nucl_gb.accession2taxid \
        --output "$output_file" \
        --ranks 'superkingdom;phylum;class;order;family;genus;species'; then
        echo "Successfully completed $order"
    else
        echo "ERROR: Failed to process $order"
        # remove the output file if it exists but the command failed
        [[ -f "$output_file" ]] && rm "$output_file"
        return 1
    fi
}

# export the function to parallel 
export -f process_order

#orders array
orders=("Entomobryomorpha" "Neelipleona" "Poduromorpha" "Symphypleona" "Rhabdura" "Dicellurata" "Acerentomata" "Eosentomata" "Sinentomata"  "Zygentoma" "Archaeognatha" "Lepidoptera" "Trichoptera" "Coleoptera" "Diptera" "Hymenoptera" "Mecoptera" "Megaloptera" "Neuroptera" "Raphidioptera" "Siphonaptera" "Strepsiptera" "Hemiptera" "Psocodea" "Thysanoptera" "Dermaptera" "Dictyoptera" "Embioptera" "Grylloblattidae" "Mantophasmatodea" "Austrophasmatidae" "Mantophasmatidae" "Orthoptera" "Phasmatodea" "Plecoptera" "Zoraptera" "Ephemeroptera" "Odonata" "Anomopoda" "Ctenopoda" "Haplopoda" "Onychopoda" "Cyclestherida" "Laevicaudata" "Spinicaudata" "Notostraca" "Anostraca" "Brachypoda" "Copepoda" "Tantulocarida" "Deoterthridae" "Microdajidae" "Amphionidacea" "Decapoda" "Euphausiacea" "Stomatopoda" "Amphipoda" "Cumacea" "Ingolfiellida" "Isopoda" "Mictacea" "Mysidacea" "Spelaeogriphacea" "Tanaidacea" "Thermosbaenacea" "Anaspidacea" "Bathynellacea" "Leptostraca" "Dendrogastrida" "Laurida" "Cryptophialida" "Lithoglyptida" "Chthamalophilidae" "Clistosaccidae" "Lernaeodiscidae" "Mycetomorphidae" "Peltogastridae" "Polyascidae" "Polysaccidae" "Sacculinidae" "Sylonidae" "Thompsoniidae" "Iblomorpha" "Balanomorpha" "Calanticomorpha" "Pollicipedomorpha" "Scalpellomorpha" "Verrucomorpha" "Facetotecta" "Hansenocaris" "Branchiura" "Cephalobaenida" "Porocephalida" "Raillietiellida" "Reighardiida" "Mystacocaridida" "Halocyprida" "Myodocopida" "Platycopida" "Podocopida" "Nectiopoda" "Sarcoptiformes" "Trombidiformes" "Opilioacariformes" "Opilioacaridae" "Holothyrida" "Ixodida" "Mesostigmata" "Araneae" "Opiliones" "Amblypygi" "Palpigradi" "Pseudoscorpiones" "Ricinulei" "Schizomida" "Scorpiones" "Solifugae" "Uropygi" "Xiphosura" "Limulidae" "Pantopoda" "Chilopoda" "Diplopoda" "Pauropoda" "Symphyla" "Apioceridae" "Apsilocephalidae" "Asilidae" "Bombyliidae" "Evocoidae" "Hilarimorphidae" "Mydidae" "Mythicomyiidae" "Ocoidae" "Scenopinidae" "Therevidae" "Ironomyiidae" "Lonchopteridae" "Phoridae" "Platypezidae" "Sciadoceridae" "Pipunculidae" "Syrphidae" "Australimyzidae" "Braulidae" "Canacidae" "Carnidae" "Chloropidae" "Cryptochetidae" "Inbiomyiidae" "Milichiidae" "Tethinidae" "Conopidae" "Diopsidae" "Megamerinidae" "Nothybidae" "Psilidae" "Somatiidae" "Strongylophthalmyiidae" "Syringogastridae" "Tanypezidae" "Camillidae" "Campichoetidae" "Curtonotidae" "Diastatidae" "Drosophilidae" "Ephydridae" "Celyphidae" "Chamaemyiidae" "Lauxaniidae" "Cypselosomatidae" "Micropezidae" "Neriidae" "Pseudopomyzidae" "Acartophthalmidae" "Agromyzidae" "Anthomyzidae" "Asteiidae" "Aulacigastridae" "Clusiidae" "Fergusoninidae" "Marginidae" "Nannodastiidae" "Neminidae" "Neurochaetidae" "Odiniidae" "Opomyzidae" "Paraleucopidae" "Periscelididae" "Teratomyzidae" "Xenasteiidae" "Coelopidae" "Dryomyzidae" "Helcomyzidae" "Helosciomyzidae" "Heterocheilidae" "Natalimyzidae" "Sciomyzidae" "Sepsidae" "Chyromyidae" "Heleomyzidae" "Mormotomyiidae" "Sphaeroceridae" "Ctenostylidae" "Lonchaeidae" "Pallopteridae" "Piophilidae" "Platystomatidae" "Pyrgotidae" "Richardiidae" "Tephritidae" "Ulidiidae" "Glossinidae" "Hippoboscidae" "Nycteribiidae" "Streblidae" "Anthomyiidae" "Fanniidae" "Muscidae" "Scathophagidae" "Calliphoridae" "Mystacinobiidae" "Oestridae" "Polleniidae" "Rhiniidae" "Rhinophoridae" "Sarcophagidae" "Tachinidae" "Ulurumyiidae" "Atelestidae" "Brachystomatidae" "Dolichopodidae" "Empididae" "Hybotidae" "Acroceridae" "Nemestrinidae" "Stratiomyidae" "Xylomyidae" "Athericidae" "Austroleptidae" "Oreoleptidae" "Pelecorhynchidae" "Rhagionidae" "Tabanidae" "Vermileonidae" "Pantophthalmidae" "Xylophagaidae" "Axymyiidae" "Bibionidae" "Hesperinidae" "Pleciidae" "Pachyneuridae" "Bolitophilidae" "Cecidomyiidae" "Diadocidiidae" "Ditomyiidae" "Keroplatidae" "Mycetophilidae" "Rangomaramidae" "Sciaridae" "Blephariceridae" "Deuterophlebiidae" "Nymphomyiidae" "Ceratopogonidae" "Chironomidae" "Simuliidae" "Thaumaleidae" "Chaoboridae" "Culicidae" "Dixidae" "Psychodidae" "Canthyloscelidae" "Scatopsidae" "Anisopodidae" "Perissommatidae" "Synneuridae" "Trichoceridae" "Ptychopteridae" "Tanyderidae" "Cylindrotomidae" "Limoniidae" "Pediciidae" "Tipulidae" "Annelida" "Mollusca" "Nematoda" "Tardigrada")


# count total and already completed orders
total_orders=${#orders[@]}
completed_count=0
for order in "${orders[@]}"; do
    if [[ -s "processed/bold_${order}_imported.txt" ]]; then
        ((completed_count++))
    fi
done

echo "progress: $completed_count/$total_orders orders already completed"
echo "timestamp: $(date)"

# Run in parallel, adjust -j for number of cores
# Removed --halt option since we now handle empty files gracefully
printf '%s\n' "${orders[@]}" | parallel --progress -j 20 process_order

echo "completed at: $(date)"

# Final summary
final_completed=0
failed_orders=()
for order in "${orders[@]}"; do
    if [[ -s "processed/bold_${order}_imported.txt" ]]; then
        ((final_completed++))
    else
        failed_orders+=("$order")
    fi
done

echo "successfully completed: $final_completed/$total_orders"
if [[ ${#failed_orders[@]} -gt 0 ]]; then
    echo "  failed orders: ${failed_orders[*]}"
fi


##############################################################################
#*                 STEP 2 - merge into a single database
##############################################################################
# now combine all BOLD imports
bold_files=$(ls processed/bold_*_imported.txt | tr '\n' ';' | sed 's/;$//')
crabs --merge \
  --input "${bold_files}" \
  --uniq \
  --output processed/bold_merged.txt
# Results | Written 12601911 sequences to processed/bold_merged.txt by merging 280 files containing 12677411 sequences (99.4%)

# then, merge all data sources with duplicate removal
# the--uniq parameter retains only a single version of each accession number
crabs --merge \
  --input 'processed/bold_merged.txt;processed/ncbi_imported.txt;processed/midori_total_imported.txt' \
  --uniq \
  --output processed/full_merged.txt
# Results | Written 13851732 sequences to processed/full_merged.txt by merging 3 files containing 18658377 sequences (74.24%)

