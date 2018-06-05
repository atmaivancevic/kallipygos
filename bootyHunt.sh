# Find SNP region in mammal genomes

module load BLAST+/2.6.0-foss-2016uofa-Python-2.7.11

cd /data/rc003/atma/2016/LASTZ_extraction/Genomes/mammalia

for i in Ailuropoda.melanoleuca Balaenoptera.acutorostrata Bison.bison Bos.indicus Bos.mutus Bos.taurus Bubalus.bubalis Callithrix.jacchus Camelus.dromedarius Camelus.ferus Canis.lupus Capra.hircus Cavia.aperea Cavia.porcellus Ceratotherium.simum Chinchilla.lanigera Chlorocebus.sabaeus Choloepus.hoffmanni Chrysochloris.asiatica Condylura.cristata Cricetulus.griseus Dasypus.novemcinctus Dipodomys.ordii Echinops.telfairi Eidolon.helvum Elephantulus.edwardii Eptesicus.fuscus Equus.caballus.Mongolian Equus.caballus.Thoroughbred Equus.przewalskii Erinaceus.europaeus Felis.catus Fukomys.damarensis Galeopterus.variegatus Gorilla.gorilla Heterocephalus.glaber Homo.sapiens Ictidomys.tridecemlineatus Jaculus.jaculus Leptonychotes.weddellii Lipotes.vexillifer Loxodonta.africana Macaca.fascicularis Macaca.mulatta Macropus.eugenii Manis.pentadactyla Megaderma.lyra Mesocricetus.auratus Microcebus.murinus Microtus.ochrogaster Monodelphis.domestica Mus.musculus Mustela.putorius Myotis.brandtii Myotis.davidii Myotis.lucifugus Nannospalax.galili Nasalis.larvatus Nomascus.leucogenys Ochotona.princeps Octodon.degus Odobenus.rosmarus Orcinus.orca Ornithorhynchus.anatinus Orycteropus.afer Oryctolagus.cuniculus Otolemur.garnettii Ovis.aries.musimon Ovis.aries.Texel Pan.paniscus Panthera.tigris Pantholops.hodgsonii Pan.troglodytes Papio.anubis Peromyscus.maniculatus Physeter.catodon Pongo.abelii Procavia.capensis Pteronotus.parnellii Pteropus.alecto Pteropus.vampyrus Rattus.norvegicus Rhinolophus.ferrumequinum Rhinopithecus.roxellana Saimiri.boliviensis Sarcophilus.harrisii Sorex.araneus Sus.scrofa.Duroc Sus.scrofa.Minipig Sus.scrofa.Tibetan Tachyglossus.aculeatus Tarsius.syrichta Trichechus.manatus Tupaia.belangeri Tupaia.chinensis Tursiops.truncatus Ursus.maritimus Vicugna.pacos;
do 
blastn -query MEG3_DLK1_genes.fasta -db "$i"/*.fa -evalue 1e-5 -outfmt 6 -out "$i"_genes.out;
done

# Convert to bed format, merge MEG3 and DLK1 genes, extract the whole region
# Warning: includes some super dodgy text parsing

module load BEDTools/2.25.0-foss-2015b

for i in Ailuropoda.melanoleuca Balaenoptera.acutorostrata Bison.bison Bos.indicus Bos.mutus Bos.taurus Bubalus.bubalis Callithrix.jacchus Camelus.dromedarius Camelus.ferus Canis.lupus Capra.hircus Cavia.aperea Cavia.porcellus Ceratotherium.simum Chinchilla.lanigera Chlorocebus.sabaeus Choloepus.hoffmanni Chrysochloris.asiatica Condylura.cristata Cricetulus.griseus Dasypus.novemcinctus Dipodomys.ordii Echinops.telfairi Eidolon.helvum Elephantulus.edwardii Eptesicus.fuscus Equus.caballus.Mongolian Equus.caballus.Thoroughbred Equus.przewalskii Erinaceus.europaeus Felis.catus Fukomys.damarensis Galeopterus.variegatus Gorilla.gorilla Heterocephalus.glaber Homo.sapiens Ictidomys.tridecemlineatus Jaculus.jaculus Leptonychotes.weddellii Lipotes.vexillifer Loxodonta.africana Macaca.fascicularis Macaca.mulatta Macropus.eugenii Manis.pentadactyla Megaderma.lyra Mesocricetus.auratus Microcebus.murinus Microtus.ochrogaster Monodelphis.domestica Mus.musculus Mustela.putorius Myotis.brandtii Myotis.davidii Myotis.lucifugus Nannospalax.galili Nasalis.larvatus Nomascus.leucogenys Ochotona.princeps Octodon.degus Odobenus.rosmarus Orcinus.orca Ornithorhynchus.anatinus Orycteropus.afer Oryctolagus.cuniculus Otolemur.garnettii Ovis.aries.musimon Ovis.aries.Texel Pan.paniscus Panthera.tigris Pantholops.hodgsonii Pan.troglodytes Papio.anubis Peromyscus.maniculatus Physeter.catodon Pongo.abelii Procavia.capensis Pteronotus.parnellii Pteropus.alecto Pteropus.vampyrus Rattus.norvegicus Rhinolophus.ferrumequinum Rhinopithecus.roxellana Saimiri.boliviensis Sarcophilus.harrisii Sorex.araneus Sus.scrofa.Duroc Sus.scrofa.Minipig Sus.scrofa.Tibetan Tachyglossus.aculeatus Tarsius.syrichta Trichechus.manatus Tupaia.belangeri Tupaia.chinensis Tursiops.truncatus Ursus.maritimus Vicugna.pacos;
do 

cat "$i"_genes.out | awk '{print $2 " " $9 " " $10 " " $1}' | awk '{if ($2 < $3) print $0 " " "+"}' > "$i".plus.tmp;

cat "$i"_genes.out | awk '{print $2 " " $9 " " $10 " " $1}' | awk '{if ($2 > $3) print $1 " " $3 " " $2 " " $4 " " "-"}' > "$i".minus.tmp;

cat "$i".plus.tmp "$i".minus.tmp | awk '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" "1" "\t" $5}' > "$i".both.tmp;

bedtools sort -i "$i".both.tmp | bedtools merge -s -i test.both.sorted.tmp -c 4 -o collapse -d 1000000 | awk '{print $1 "\t" $2 "\t" $3 "\t" $5 "\t" "1" "\t" $4}' > "$i"_genes.bed;

bedtools getfasta -s -fi "$i"/*.fa -bed "$i"_genes.bed -fo "$i"_genes.fasta;

done

# remove all those temporary files
rm *.tmp

# Download to local comp
# Align to 143bp SNP region from Callipyge sheep paper (Freking et al, 2002)
# Check if A or G in SNP position
