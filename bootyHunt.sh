# Find SNP region in mammal genomes

module load BLAST+/2.6.0-foss-2016uofa-Python-2.7.11

cd /data/rc003/atma/2016/LASTZ_extraction/Genomes/mammalia

for i in Ailuropoda.melanoleuca Balaenoptera.acutorostrata Bison.bison Bos.indicus Bos.mutus Bos.taurus Bubalus.bubalis Callithrix.jacchus Camelus.dromedarius Camelus.ferus Canis.lupus Capra.hircus Cavia.aperea Cavia.porcellus Ceratotherium.simum Chinchilla.lanigera Chlorocebus.sabaeus Choloepus.hoffmanni Chrysochloris.asiatica Condylura.cristata Cricetulus.griseus Dasypus.novemcinctus Dipodomys.ordii Echinops.telfairi Eidolon.helvum Elephantulus.edwardii Eptesicus.fuscus Equus.caballus.Mongolian Equus.caballus.Thoroughbred Equus.przewalskii Erinaceus.europaeus Felis.catus Fukomys.damarensis Galeopterus.variegatus Gorilla.gorilla Heterocephalus.glaber Homo.sapiens Ictidomys.tridecemlineatus Jaculus.jaculus Leptonychotes.weddellii Lipotes.vexillifer Loxodonta.africana Macaca.fascicularis Macaca.mulatta Macropus.eugenii Manis.pentadactyla Megaderma.lyra Mesocricetus.auratus Microcebus.murinus Microtus.ochrogaster Monodelphis.domestica Mus.musculus Mustela.putorius Myotis.brandtii Myotis.davidii Myotis.lucifugus Nannospalax.galili Nasalis.larvatus Nomascus.leucogenys Ochotona.princeps Octodon.degus Odobenus.rosmarus Orcinus.orca Ornithorhynchus.anatinus Orycteropus.afer Oryctolagus.cuniculus Otolemur.garnettii Ovis.aries.musimon Ovis.aries.Texel Pan.paniscus Panthera.tigris Pantholops.hodgsonii Pan.troglodytes Papio.anubis Peromyscus.maniculatus Physeter.catodon Pongo.abelii Procavia.capensis Pteronotus.parnellii Pteropus.alecto Pteropus.vampyrus Rattus.norvegicus Rhinolophus.ferrumequinum Rhinopithecus.roxellana Saimiri.boliviensis Sarcophilus.harrisii Sorex.araneus Sus.scrofa.Duroc Sus.scrofa.Minipig Sus.scrofa.Tibetan Tachyglossus.aculeatus Tarsius.syrichta Trichechus.manatus Tupaia.belangeri Tupaia.chinensis Tursiops.truncatus Ursus.maritimus Vicugna.pacos;
do 
blastn -query MEG3_DLK1_genes.fasta -db "$i"/*.fa -evalue 1e-5 -outfmt 6 -out "$i"_genes.out;
done

# Take the top hit for each gene and extract the whole region

module load BEDTools/2.25.0-foss-2015b

for i in Ailuropoda.melanoleuca Balaenoptera.acutorostrata Bison.bison Bos.indicus Bos.mutus Bos.taurus Bubalus.bubalis Callithrix.jacchus Camelus.dromedarius Camelus.ferus Canis.lupus Capra.hircus Cavia.aperea Cavia.porcellus Ceratotherium.simum Chinchilla.lanigera Chlorocebus.sabaeus Choloepus.hoffmanni Chrysochloris.asiatica Condylura.cristata Cricetulus.griseus Dasypus.novemcinctus Dipodomys.ordii Echinops.telfairi Eidolon.helvum Elephantulus.edwardii Eptesicus.fuscus Equus.caballus.Mongolian Equus.caballus.Thoroughbred Equus.przewalskii Erinaceus.europaeus Felis.catus Fukomys.damarensis Galeopterus.variegatus Gorilla.gorilla Heterocephalus.glaber Homo.sapiens Ictidomys.tridecemlineatus Jaculus.jaculus Leptonychotes.weddellii Lipotes.vexillifer Loxodonta.africana Macaca.fascicularis Macaca.mulatta Macropus.eugenii Manis.pentadactyla Megaderma.lyra Mesocricetus.auratus Microcebus.murinus Microtus.ochrogaster Monodelphis.domestica Mus.musculus Mustela.putorius Myotis.brandtii Myotis.davidii Myotis.lucifugus Nannospalax.galili Nasalis.larvatus Nomascus.leucogenys Ochotona.princeps Octodon.degus Odobenus.rosmarus Orcinus.orca Ornithorhynchus.anatinus Orycteropus.afer Oryctolagus.cuniculus Otolemur.garnettii Ovis.aries.musimon Ovis.aries.Texel Pan.paniscus Panthera.tigris Pantholops.hodgsonii Pan.troglodytes Papio.anubis Peromyscus.maniculatus Physeter.catodon Pongo.abelii Procavia.capensis Pteronotus.parnellii Pteropus.alecto Pteropus.vampyrus Rattus.norvegicus Rhinolophus.ferrumequinum Rhinopithecus.roxellana Saimiri.boliviensis Sarcophilus.harrisii Sorex.araneus Sus.scrofa.Duroc Sus.scrofa.Minipig Sus.scrofa.Tibetan Tachyglossus.aculeatus Tarsius.syrichta Trichechus.manatus Tupaia.belangeri Tupaia.chinensis Tursiops.truncatus Ursus.maritimus Vicugna.pacos;
do 
bedtools getfasta -s -fi "$i"/*.fa -bed "$i"_genes.bed -fo "$i"_genes.fasta;
done

# Align to 143bp SNP region from Callipyge sheep paper (Freking et al, 2002)
# And check if A or G
