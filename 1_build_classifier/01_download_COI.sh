# download COI database using CRABS
# katherine carbeck
# 07 sep 2025

##############################################################################
#*          STEP 1 - download NCBI taxonomy
##############################################################################
#export path
export PATH="/lustre2/home/lc736_0001/diet/songbird_coi_database/reference_database_creator:$PATH"

# direct path to crabs
crabs --help

cd songbird_coi_database

# first, download NCBI taxonomy files - these are essential for taxonomic validation
crabs --download-taxonomy --output taxonomy/

##############################################################################
#*                 STEP 2 - download COI seqs
##############################################################################
# BOLD Systems download (primary COI source)
orders=("Entomobryomorpha" "Neelipleona" "Poduromorpha" "Symphypleona" "Rhabdura" "Dicellurata" "Acerentomata" "Eosentomata" "Sinentomata"  "Zygentoma" "Archaeognatha" "Lepidoptera" "Trichoptera" "Coleoptera" "Diptera" "Hymenoptera" "Mecoptera" "Megaloptera" "Neuroptera" "Raphidioptera" "Siphonaptera" "Strepsiptera" "Hemiptera" "Psocodea" "Thysanoptera" "Dermaptera" "Dictyoptera" "Embioptera" "Grylloblattidae" "Mantophasmatodea" "Austrophasmatidae" "Mantophasmatidae" "Orthoptera" "Phasmatodea" "Plecoptera" "Zoraptera" "Ephemeroptera" "Odonata" "Anomopoda" "Ctenopoda" "Haplopoda" "Onychopoda" "Cyclestherida" "Laevicaudata" "Spinicaudata" "Notostraca" "Anostraca" "Brachypoda" "Copepoda" "Tantulocarida" "Deoterthridae" "Microdajidae" "Amphionidacea" "Decapoda" "Euphausiacea" "Stomatopoda" "Amphipoda" "Cumacea" "Ingolfiellida" "Isopoda" "Mictacea" "Mysidacea" "Spelaeogriphacea" "Tanaidacea" "Thermosbaenacea" "Anaspidacea" "Bathynellacea" "Leptostraca" "Dendrogastrida" "Laurida" "Cryptophialida" "Lithoglyptida" "Chthamalophilidae" "Clistosaccidae" "Lernaeodiscidae" "Mycetomorphidae" "Peltogastridae" "Polyascidae" "Polysaccidae" "Sacculinidae" "Sylonidae" "Thompsoniidae" "Iblomorpha" "Balanomorpha" "Calanticomorpha" "Pollicipedomorpha" "Scalpellomorpha" "Verrucomorpha" "Facetotecta" "Hansenocaris" "Branchiura" "Cephalobaenida" "Porocephalida" "Raillietiellida" "Reighardiida" "Mystacocaridida" "Halocyprida" "Myodocopida" "Platycopida" "Podocopida" "Nectiopoda" "Sarcoptiformes" "Trombidiformes" "Opilioacariformes" "Opilioacaridae" "Holothyrida" "Ixodida" "Mesostigmata" "Araneae" "Opiliones" "Amblypygi" "Palpigradi" "Pseudoscorpiones" "Ricinulei" "Schizomida" "Scorpiones" "Solifugae" "Uropygi" "Xiphosura" "Limulidae" "Pantopoda" "Chilopoda" "Diplopoda" "Pauropoda" "Symphyla" "Apioceridae" "Apsilocephalidae" "Asilidae" "Bombyliidae" "Evocoidae" "Hilarimorphidae" "Mydidae" "Mythicomyiidae" "Ocoidae" "Scenopinidae" "Therevidae" "Ironomyiidae" "Lonchopteridae" "Phoridae" "Platypezidae" "Sciadoceridae" "Pipunculidae" "Syrphidae" "Australimyzidae" "Braulidae" "Canacidae" "Carnidae" "Chloropidae" "Cryptochetidae" "Inbiomyiidae" "Milichiidae" "Tethinidae" "Conopidae" "Diopsidae" "Megamerinidae" "Nothybidae" "Psilidae" "Somatiidae" "Strongylophthalmyiidae" "Syringogastridae" "Tanypezidae" "Camillidae" "Campichoetidae" "Curtonotidae" "Diastatidae" "Drosophilidae" "Ephydridae" "Celyphidae" "Chamaemyiidae" "Lauxaniidae" "Cypselosomatidae" "Micropezidae" "Neriidae" "Pseudopomyzidae" "Acartophthalmidae" "Agromyzidae" "Anthomyzidae" "Asteiidae" "Aulacigastridae" "Clusiidae" "Fergusoninidae" "Marginidae" "Nannodastiidae" "Neminidae" "Neurochaetidae" "Odiniidae" "Opomyzidae" "Paraleucopidae" "Periscelididae" "Teratomyzidae" "Xenasteiidae" "Coelopidae" "Dryomyzidae" "Helcomyzidae" "Helosciomyzidae" "Heterocheilidae" "Natalimyzidae" "Sciomyzidae" "Sepsidae" "Chyromyidae" "Heleomyzidae" "Mormotomyiidae" "Sphaeroceridae" "Ctenostylidae" "Lonchaeidae" "Pallopteridae" "Piophilidae" "Platystomatidae" "Pyrgotidae" "Richardiidae" "Tephritidae" "Ulidiidae" "Glossinidae" "Hippoboscidae" "Nycteribiidae" "Streblidae" "Anthomyiidae" "Fanniidae" "Muscidae" "Scathophagidae" "Calliphoridae" "Mystacinobiidae" "Oestridae" "Polleniidae" "Rhiniidae" "Rhinophoridae" "Sarcophagidae" "Tachinidae" "Ulurumyiidae" "Atelestidae" "Brachystomatidae" "Dolichopodidae" "Empididae" "Hybotidae" "Acroceridae" "Nemestrinidae" "Stratiomyidae" "Xylomyidae" "Athericidae" "Austroleptidae" "Oreoleptidae" "Pelecorhynchidae" "Rhagionidae" "Tabanidae" "Vermileonidae" "Pantophthalmidae" "Xylophagaidae" "Axymyiidae" "Bibionidae" "Hesperinidae" "Pleciidae" "Pachyneuridae" "Bolitophilidae" "Cecidomyiidae" "Diadocidiidae" "Ditomyiidae" "Keroplatidae" "Mycetophilidae" "Rangomaramidae" "Sciaridae" "Blephariceridae" "Deuterophlebiidae" "Nymphomyiidae" "Ceratopogonidae" "Chironomidae" "Simuliidae" "Thaumaleidae" "Chaoboridae" "Culicidae" "Dixidae" "Psychodidae" "Canthyloscelidae" "Scatopsidae" "Anisopodidae" "Perissommatidae" "Synneuridae" "Trichoceridae" "Ptychopteridae" "Tanyderidae" "Cylindrotomidae" "Limoniidae" "Pediciidae" "Tipulidae" "Annelida" "Mollusca" "Nematoda" "Tardigrada")


# bash script to loop through each taxonomic group and download sequences using CRABS:
{
for i in "${orders[@]}"; do
  echo "Downloading: ${i}"
  crabs --download-bold \
    --taxon "${i}" \
    --marker 'COI-5P' \
    --output "raw_downloads/bold_${i}.fasta" || echo "Failed to download ${i}"
done 
} &
# Bold_raw_seqs.fasta

# NCBI GenBank comprehensive arthropod COI
{
crabs --download-ncbi \
  --query '("Arthropoda"[Organism] OR Arthropoda[All Fields]) AND (mitochondrion[filter] AND ("100"[SLEN] : "25000"[SLEN]) AND (COI[Gene Name] OR "cytochrome oxidase subunit I"[Gene Name] OR "cytochrome c oxidase subunit 1"[Gene Name]))' \
  --output raw_downloads/ncbi_arthropoda_coi.fasta \
  --email kmc464@cornell.edu \
  --database nucleotide \
  --batchsize 5000
} &

# MIDORI 
{
crabs --download-midori \
  --output raw_downloads/cytb_total_267.fasta \
  --gb-number 267_2025-06-19 \
  --gene CO1 \
  --gb-type total
} &








