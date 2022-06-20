# Illumina-paired-end-simulator

Simulator je kreiran u okviru projekta na kursu Genomske informatike 2022 - Elektrotehnički fakultet, Univerzitet u Beogradu

Simulator sekvencera kreira veštačke Illumina paired-end read-ove uzimajući nukleotide iz referentnog genoma koji je zadat na ulazu.
Kao ulazne parametre simulator prima prosečan kvalitet nukleotida (`avg_quality`), pokrivenost (`coverage`), dužinu read-ova (`read_length`), dužinu inserta (`insert_size`), kao i error rate - verovatnoća pojavljivanja pogrešnih nukleotida u okviru samih readova i to posebno za snipove (`error_rate_snip`), za insercije (`error_rate_ins`) i delecije (`error_rate_del`).
Za simulaciju umnožavanja klastera, uveden je i parametar `max_diff`, koji određuje koliko najviše nukleotida neki fragment može da porani ili zakasni pri očitavanju.

Izlaz simulatora su dva FASTQ fajla u kojima su upisani kreirani read-ovi i sam fajl koji sadrži poravnate sve read-ove iz FASTQ fajlova zajedno sa pozicijama iz referentnog genoma sa kojih su ekstrahovani nukleotidi.

Jedan referentni genom za testiranje je preuzet iz NIH baze (NCBI Reference Sequence: `NZ_OW968172.1`) i još dva genoma su preuzeta sa sbgenomics platforme (`PhiX_genome.fasta` i `example_human_reference.fasta`).

Koristeći fajlove koje je napravio simulator, izvršeno je testiranje kvaliteta alajmenta BWA-MEM i Bowtie alata i rezultati testiranja su prikazani grafički koristeći klasifikacione metrike.

## Sadržaj repozitorijuma:
`illumina_simulator.py` - opisani simulator

`compare.py` - vrši poređenje bwa-mem i bowtie fajlova sa referentnim sam fajlom i metrike upisuje u txt fajl

`graphs.py` - učitava metrike iz txt fajla i grafički ih prikazuje

## Prezentacija projekta:

[Simulator Illumina paired-end sekvencera](https://youtu.be/EmR55vqCvcc)

### Autori:

Sara Živković 2021/3285

Božidar Obradović 2021/3081
