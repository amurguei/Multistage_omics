#Andromeda 
ssh -l amalia.murgueitiocal ssh3.hac.uri.edu
cd /data/putnamlab/Amalia_Murgueitio/orthogroups/refs

#Hive: 
ssh amalia@hive02.haifa.ac.il
cd /lustre1/home/mass/amalia/orthofinder_refs

/lustre1/home/mass/amalia/.conda/envs/orthofinder #orthofinder location


interactive
curl -O http://cyanophora.rutgers.edu/montipora/Montipora_capitata_HIv3.genes.pep.faa.gz && gunzip Montipora_capitata_HIv3.genes.pep.faa.gz
curl -O http://cyanophora.rutgers.edu/Pocillopora_acuta/Pocillopora_acuta_HIv2.genes.pep.faa.gz && gunzip Pocillopora_acuta_HIv2.genes.pep.faa.gz
curl -O http://spis.reefgenomics.org/download/Spis.genome.annotation.pep.longest.fa.gz && gunzip Spis.genome.annotation.pep.longest.fa.gz
#Note: I was trying to download from the URL, it´s much better to do so from the ftp at NCBI (I found out late)
wget -O - https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/949/145/GCF_001949145.1_OKI-Apl_1.0/GCF_001949145.1_OKI-Apl_1.0_protein.faa.gz | gunzip > Acanthaster_planci_protein.faa
wget -O - https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/009/805/615/GCA_009805615.1_SNU_Aamp_1/GCA_009805615.1_SNU_Aamp_1_protein.faa.gz | gunzip > Amphibalanus_amphitrite_protein.faa
wget -O - https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/090/795/GCF_000090795.2_v1.1/GCF_000090795.2_v1.1_protein.faa.gz | gunzip > Amphimedon_queenslandica_protein.faa
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/222/465/GCF_000222465.1_Adig_1.1/GCF_000222465.1_Adig_1.1_protein.faa.gz -O - | gunzip > GCF_000222465.1_Acropora_digitifera_1.1_protein.faa
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/602/425/GCF_009602425.1_ASM960242v1/GCF_009602425.1_ASM960242v1_protein.faa.gz -O - | gunzip > GCF_009602425.1_Actinia_tenebrosa_ASM960242v1_protein.faa
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/417/965/GCF_001417965.1_Aiptasia_genome_1.1/GCF_001417965.1_Aiptasia_genome_1.1_protein.faa.gz -O - | gunzip > GCF_001417965.1_Exaiptasia_diaphana_protein.faa
wget "https://aniseed.fr/aniseed/download/?file=data%2Fboleac%2FBoleac_proteins_v4_fasta.zip&module=aniseed&action=download:index" -O - | funzip > Botrylloides_leachii_proteins_v4.fasta
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/013/753/865/GCF_013753865.1_Amil_v2.1/GCF_013753865.1_Amil_v2.1_protein.faa.gz -O - | gunzip > GCF_013753865.1_Acropora_millepora_v2.1_protein.faa
wget -O - http://corallimorpharia.reefgenomics.org/download/afen.prot.fa.gz | gunzip > Amplexidiscus_fenestrafer_protein.faa
wget -O - http://corallimorpharia.reefgenomics.org/download/dspp.prot.fa.gz | gunzip > Discosoma_sp_protein.faa
wget -O - http://aten.reefgenomics.org/download/aten_0.11.maker_post_001.proteins.fasta.gz | gunzip > Acropora_tenuifolia_0.11.maker_post_001.proteins.fasta
wget -O - http://ffun.reefgenomics.org/download/ffun_1.0.proteins.fasta.gz | gunzip > Fungia_spp_1.0_proteins.fasta
wget -O - http://gfas.reefgenomics.org/download/gfas_1.0.proteins.fasta.gz | gunzip > Galaxea_fascicularis_1.0_proteins.fasta
wget -O - http://gasp.reefgenomics.org/download/gasp_1.0.proteins.fasta.gz | gunzip > Goniastrea_aspera_1.0_proteins.fasta
wget -O - http://plut.reefgenomics.org/download/plut2v1.1.proteins.fasta.gz | gunzip > Porites_lutea_2v1.1_proteins.fasta
wget -O - http://rmue.reefgenomics.org/download/renilla_predicted_proteins.fa.gz | gunzip > Renilla_muelleri_predicted_proteins.fa
wget -O - https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/009/887/335/GCA_009887335.1_TAU_Hsal_1/GCA_009887335.1_TAU_Hsal_1_protein.faa.gz | gunzip > Henneguya_salminicola_protein.faa
wget -O - https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/004/324/835/GCF_004324835.1_DenGig_1.0/GCF_004324835.1_DenGig_1.0_protein.faa.gz | gunzip > Dendronephthya_gigantea_protein.faa
wget -O - https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/038/396/675/GCF_038396675.1_HydraT2T_AEP/GCF_038396675.1_HydraT2T_AEP_protein.faa.gz | gunzip > Hydra_vulgaris_protein.faa
wget -O - https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/571/385/GCF_002571385.2_Stylophora_pistillata_v1.1/GCF_002571385.2_Stylophora_pistillata_v1.1_protein.faa.gz | gunzip > GCF_002571385.2_Stylophora_pistillata_v1.1_protein.faa
wget -O - https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/010/108/815/GCA_010108815.2_TAU_Msqu_1.1/GCA_010108815.2_TAU_Msqu_1.1_protein.faa.gz | gunzip > GCA_010108815.2_Myxobolus_squamalis_protein.faa
wget -O - https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/932/526/225/GCF_932526225.1_jaNemVect1.1/GCF_932526225.1_jaNemVect1.1_protein.faa.gz | gunzip > GCF_932526225.1_Nematostella_vectensis_protein.faa
wget -O - https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/042/975/GCF_002042975.1_ofav_dov_v1/GCF_002042975.1_ofav_dov_v1_protein.faa.gz | gunzip > GCF_002042975.1_Orbicella_faveolata_protein.faa
wget -O - https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/704/095/GCF_003704095.1_ASM370409v1/GCF_003704095.1_ASM370409v1_protein.faa.gz | gunzip > GCF_003704095.1_Pocillopora_damicornis_protein.faa
wget -O - https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/827/895/GCA_000827895.1_ASM82789v1/GCA_000827895.1_ASM82789v1_protein.faa.gz | gunzip > GCA_000827895.1_Thelohanellus_kitauei_protein.faa
wget -O - https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/039/355/GCF_001039355.2_LinAna2.0/GCF_001039355.2_LinAna2.0_protein.faa.gz | gunzip > GCF_001039355.2_Lingula_anatina_protein.faa
wget -O - https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/865/GCF_000002865.3_V1.0/GCF_000002865.3_V1.0_protein.faa.gz | gunzip > GCF_000002865.3_Monosiga_brevicollis_protein.faa
wget -O - https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/188/695/GCF_000188695.1_Proterospongia_sp_ATCC50818/GCF_000188695.1_Proterospongia_sp_ATCC50818_protein.faa.gz | gunzip > GCF_000188695.1_Salpingoeca_rosetta_protein.faa
wget -O - https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/235/GCF_000002235.5_Spur_5.0/GCF_000002235.5_Spur_5.0_protein.faa.gz | gunzip > GCF_000002235.5_Strongylocentrotus_purpuratus_protein.faa
wget -O - https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/004/195/GCF_000004195.4_UCB_Xtro_10.0/GCF_000004195.4_UCB_Xtro_10.0_protein.faa.gz | gunzip > GCF_000004195.4_Xenopus_tropicalis_protein.faa
wget -O "Xenia sp..fa" https://cmo.carnegiescience.edu/endosymbiosis/genome/xenSp1.proteins.fa
wget -O Danio_rerio_protein.fa.gz ftp://ftp.ensembl.org/pub/release-113/fasta/danio_rerio/pep/Danio_rerio.GRCz11.pep.all.fa.gz
wget -O - ftp://ftp.ensemblgenomes.ebi.ac.uk/pub/metazoa/release-60/fasta/drosophila_melanogaster/pep/Drosophila_melanogaster.BDGP6.46.pep.all.fa.gz | gunzip > Drosophila_melanogaster_protein.faa
wget -O - https://ftp.ensemblgenomes.ebi.ac.uk/pub/fungi/release-60/fasta/fungi_ascomycota2_collection/eremothecium_gossypii_fdag1_gca_000968835/pep/Eremothecium_gossypii_fdag1_gca_000968835.ASM96883v1.pep.all.fa.gz | gunzip > Eremothecium_gossypii_protein.faa
wget -O - https://ftp.ensemblgenomes.ebi.ac.uk/pub/fungi/release-60/fasta/fungi_basidiomycota1_collection/rhodotorula_taiwanensis_gca_002922495/pep/Rhodotorula_taiwanensis_gca_002922495.ASM292249v1.pep.all.fa.gz | gunzip > Rhodotorula_taiwanensis_ASM292249v1_protein.faa
wget -O - https://ftp.ensemblgenomes.ebi.ac.uk/pub/metazoa/release-60/fasta/apis_mellifera/pep/Apis_mellifera.Amel_HAv3.1.pep.all.fa.gz | gunzip > Apis_mellifera_protein.faa
wget -O - https://ftp.ensemblgenomes.ebi.ac.uk/pub/metazoa/release-60/fasta/caenorhabditis_elegans/pep/Caenorhabditis_elegans.WBcel235.pep.all.fa.gz | gunzip > Caenorhabditis_elegans_protein.faa
wget -O - https://ftp.ensemblgenomes.ebi.ac.uk/pub/fungi/release-60/fasta/candida_glabrata/pep/Candida_glabrata.GCA000002545v2.pep.all.fa.gz | gunzip > Candida_glabrata_GCA000002545v2_protein.faa
wget -O - https://ftp.ensemblgenomes.ebi.ac.uk/pub/fungi/release-60/fasta/saccharomyces_cerevisiae/pep/Saccharomyces_cerevisiae.R64-1-1.pep.all.fa.gz | gunzip > Saccharomyces_cerevisiae_R64-1-1_protein.faa
wget -O - https://ftp.ensemblgenomes.ebi.ac.uk/pub/fungi/release-60/fasta/fungi_ascomycota1_collection/scedosporium_apiospermum_gca_000732125/pep/Scedosporium_apiospermum_gca_000732125.ScApio1.0.pep.all.fa.gz | gunzip > Scedosporium_apiospermum_GCA000732125_Proten.faa
wget -O - https://ftp.ensemblgenomes.ebi.ac.uk/pub/fungi/release-60/fasta/schizosaccharomyces_pombe/pep/Schizosaccharomyces_pombe.ASM294v2.pep.all.fa.gz | gunzip > Schizosaccharomyces_pombe_ASM294v2_protein.faa
wget -O - https://ftp.ensemblgenomes.ebi.ac.uk/pub/fungi/release-60/fasta/yarrowia_lipolytica/pep/Yarrowia_lipolytica.GCA_000002525.1.pep.all.fa.gz | gunzip > Yarrowia_lipolytica_GCA000002525_Protein.faa
wget -O - https://ftp.ensembl.org/pub/release-113/fasta/homo_sapiens/pep/Homo_sapiens.GRCh38.pep.all.fa.gz | gunzip > Homo_sapiens_GRCh38_protein.faa
wget -O - https://ftp.ensemblgenomes.ebi.ac.uk/pub/metazoa/release-60/fasta/lottia_gigantea/pep/Lottia_gigantea.Lotgi1.pep.all.fa.gz | gunzip > Lottia_gigantea_Lotgi1_protein.faa
wget -O - https://ftp.ensemblgenomes.ebi.ac.uk/pub/metazoa/release-60/fasta/mnemiopsis_leidyi/pep/Mnemiopsis_leidyi.MneLei_Aug2011.pep.all.fa.gz | gunzip > Mnemiopsis_leidyi_MneLei_Aug2011_protein.faa
wget -O - https://ftp.ensembl.org/pub/release-113/fasta/mus_musculus/pep/Mus_musculus.GRCm39.pep.all.fa.gz | gunzip > Mus_musculus_GRCm39_protein.faa
wget -O - https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/516/915/GCA_000516915.1_OphHan1.0/GCA_000516915.1_OphHan1.0_protein.faa.gz | gunzip > GCA_000516915.1_Ophiophagus_hannah_protein.faa
wget -O - https://ftp.ensemblgenomes.ebi.ac.uk/pub/metazoa/release-60/fasta/strigamia_maritima/pep/Strigamia_maritima.Smar1.pep.all.fa.gz | gunzip > Strigamia_maritima_Smar1_protein.faa

#My additions: 
#Corals
wget -O - https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/029/204/205/GCA_029204205.1_Loph_1.0/GCA_029204205.1_Loph_1.0_protein.faa.gz | gunzip > Desmophyllum_pertusum_Loph_1.0_protein.faa
wget -O - https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/942/486/035/GCA_942486035.1_PLOB_v1/GCA_942486035.1_PLOB_v1_protein.faa.gz | gunzip > Porites_lobata_PLOB_v1_protein.faa
wget -O - https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/942/486/025/GCA_942486025.1_PEVE_v1/GCA_942486025.1_PEVE_v1_protein.faa.gz | gunzip > Porites_evermanni_PEVE_v1_protein.faa
wget -O - https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/036/669/915/GCF_036669915.1_ASM3666991v2/GCF_036669915.1_ASM3666991v2_protein.faa.gz | gunzip > Pocillopora_verrucosa_ASM3666991v2_protein.faa
wget -O - https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/704/095/GCF_003704095.1_ASM370409v1/GCF_003704095.1_ASM370409v1_protein.faa.gz | gunzip > Pocillopora_damicornis_ASM370409v1_protein.faa
wget -O - https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/942/486/045/GCA_942486045.1_PMEA_v1/GCA_942486045.1_PMEA_v1_protein.faa.gz | gunzip > Pocillopora_meandrina_GCA942486045_protein.faa
wget -O - https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/036/669/905/GCF_036669905.1_ASM3666990v1/GCF_036669905.1_ASM3666990v1_protein.faa.gz | gunzip > Acropora_muricata_GCF036669905_protein.faa
wget -O - https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/042/850/425/GCA_042850425.1_Fanc_1.0/GCA_042850425.1_Fanc_1.0_protein.faa.gz | gunzip > Fimbriaphyllia_ancora_GCA042850425_protein.faa
wget -O - http://alor.reefgenomics.org/download/Acropora_loripes_predicted_proteins_v1.pep.fasta.gz | gunzip > Acropora_loripes_predicted_proteins_v1_protein.faa
wget -O - https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/032/359/415/GCA_032359415.1_NEU_Acer_K2/GCA_032359415.1_NEU_Acer_K2_protein.faa.gz | gunzip > Acropora_cervicornis_GCA032359415_protein.faa
wget -O - http://cyanophora.rutgers.edu/porites_compressa/Porites_compressa_HIv1.genes.pep.faa.gz | gunzip > Porites_compressa_HIv1_protein.faa

#Blue corals
wget -O Heliopora_caerulea_protein.faa https://figshare.com/ndownloader/files/39252302

#Anemones
wget -O - https://ftp.ensemblgenomes.ebi.ac.uk/pub/metazoa/release-60/fasta/actinia_equina_gca011057435/pep/Actinia_equina_gca011057435.equina_smartden.arrow4.noredun.pep.all.fa.gz | gunzip > Actinia_equina_protein.faa

#Hydroids
wget -O - https://ftp.ensemblgenomes.ebi.ac.uk/pub/metazoa/release-60/fasta/hydractinia_symbiolongicarpus_gca029227915v2rs/pep/Hydractinia_symbiolongicarpus_gca029227915v2rs.HSymV2.1.pep.all.fa.gz | gunzip > Hydractinia_symbiolongicarpus_GCA029227915_protein.faa
wget -O - https://ftp.ensemblgenomes.ebi.ac.uk/pub/metazoa/release-60/fasta/clytia_hemisphaerica_gca902728285/pep/Clytia_hemisphaerica_gca902728285.Clytia_hemisphaerica_genome_assembly.pep.all.fa.gz | gunzip > Clytia_hemisphaerica_gca902728285_protein.faa


nano orthofinder.sh
sbatch orthofinder.sh


#!/bin/bash
#SBATCH --job-name="Orthofinder"
#SBATCH --partition=hive7d
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=amalia.murgueitio@gmail.com
#SBATCH -D /lustre1/home/mass/amalia/.conda/envs/orthofinder

# Load miniconda3 module if required (some clusters need this to access conda)
module load miniconda3

# Activate the OrthoFinder Conda environment
source activate orthofinder

# Run OrthoFinder with the correct input directory
orthofinder -t 12 -a 12 -f /lustre1/home/mass/amalia/orthofinder_refs/

# Deactivate the environment after job completes
conda deactivate
