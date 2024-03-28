# BacSeqer

# Author: Adela Fialova
# Python 3.10.
# Last update: 28. 03. 2024


import random
import re


def read_fasta(file_path):
    """
    Načte FASTA soubor a vrátí sekvence.

    Args:
        file_path (str): Cesta k FASTA souboru.

    Returns:
        seqs (dict): Slovník, kde klíč je název sekvence a hodnota je samotná sekvence.
    """
    with open(file_path, 'r') as file:                  # otevreni souboru
        seqs = {}                                       # prazdny slovnik pro ukladani sekvenci
        current_seq = ""                                # prazdny retezec pro sestaveni ctene sekvence
        for line in file:                               # prochazeni souboru po radcich
            if line.startswith('>'):                    # hledani hlavicky
                if current_seq:                         # ulozeni nactene sekvence do slovniku seqs s klicem seq_name
                    seqs[seq_name] = current_seq
                    current_seq = ""                    # resetovani aktualni sekvence na prazdny retezec
                seq_name = line.strip()                 # nazev nove sekvence extrahovan z aktualniho radku
            else:
                current_seq += line.strip()             # pokud radek neni hlavicka, je pripojen ke current_seq
        if current_seq:
            seqs[seq_name] = current_seq                # ulozeni posledni sekvence
    return seqs                                         # vraceni slovniku se vsemi sekvencemi


def calculate_gc_content(seq):
    """
    Spočítá GC content pro načtené sekvence.

    Args:
        seq (str): Sekvence.

    Returns:
        gc_content (float): %GC content načtené sekvence.
    """
    return (seq.count('G') + seq.count('C')) / len(seq)


def extract_operons_from_gff(gff_path):
    """
    Extrahuje operony a jejich umístění ze souboru GFF.

    Args:
        gff_path (str): Cesta k souboru GFF.

    Returns:
        operons (dict): Slovník s názvy operonů jako klíči a jejich umístěními jako hodnotami.
    """
    operons = {}
    with open(gff_path, 'r') as file:
        for line in file:
            if line.startswith('#') or line.strip() == "":
                continue  # preskoceni hlavicek a prazdnych radku
            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue  # preskoceni neuplnych radku

            attributes = parts[8]

            # kontrola pritomnosti "Operon" v atributu "Note"
            note_match = re.search(r'Note=Operon\s+(\d+)', attributes)
            if note_match:
                operon_number = note_match.group(1)
                operon_name = f"Operon {operon_number}"
                location = (parts[0], int(parts[3]), int(parts[4]))  # (seqname, začátek, konec)

                # pridat umisteni do slovniku, pripojit, pokud operon uz existuje
                if operon_name in operons:
                    operons[operon_name].append(location)
                else:
                    operons[operon_name] = [location]
    return operons


def extract_rrna_from_gff(gff_path):
    """
    Zpracovává soubor GFF za účelem extrakce anotací rRNA.

    Args:
        gff_path (str): Cesta k souboru GFF.

    Returns:
        rrna_annotations (list): Seznam slovníků, kde každý reprezentuje anotaci rRNA s klíči jako 'start', 'end',
                                 'type' a 'attributes'.
    """
    rrna_annotations = []

    with open(gff_path, 'r') as file:
        for line in file:
            if line.startswith('#') or line.strip() == '':
                continue  # preskocit komentare a prazdne radky

            fields = line.strip().split('\t')
            if len(fields) == 9:
                # rozdeleni radku na jednotlive casti
                seq_id, source, feature_type, start, end, score, strand, phase, attributes = fields

                if feature_type == 'rRNA':
                    attr_dict = {}
                    for attr in attributes.split(';'):
                        # rozdeleni atributu na klice a hodnoty
                        key, value = attr.split('=')
                        attr_dict[key] = value

                    # pridani anotace rRNA do seznamu
                    rrna_annotations.append({
                        'seq_id': seq_id,  # identifikator sekvence
                        'start': int(start),  # pocatecni pozice
                        'end': int(end),  # koncova pozice
                        'strand': strand,  # smer retezce
                        'attributes': attr_dict  # dalsi atributy
                    })
    return rrna_annotations


def calculate_rrna_percentage(fasta_path, gff_path, read_fasta_func):
    """
    Vypočítá procento rRNA v sadě sekvencí.

    Args:
        fasta_path (str): Cesta k souboru FASTA.
        gff_path (str): Cesta k souboru GTF/GFF.
        read_fasta_func (funkce): Funkce pro čtení souboru FASTA (dostupná v tomto programu - read_fasta).

    Returns:
        float: Procento rRNA v celkovém množství sekvencí.
    """
    # nacteni sekvenci ze souboru FASTA
    sequences = read_fasta_func(fasta_path)
    total_sequence_length = sum(len(seq) for seq in sequences.values())

    # extrahovani anotaci rRNA ze souboru GTF/GFF
    rrna_annotations = extract_rrna_from_gff(gff_path)

    # vypocet celkove delky sekvenci rRNA
    rrna_total_length = 0
    for annotation in rrna_annotations:
        start, end = annotation['start'], annotation['end']
        rrna_total_length += end - start + 1

    # vypocet a vraceni procentualniho podilu rRNA
    return (rrna_total_length / total_sequence_length) * 100


def extract_gene_name(attributes, file_type):
    """
    Extrahuje název genu z pole atributů.
    Pracuje s formáty GTF/GFF.

    Args:
        attributes (str): Řetězec obsahující atributy.
        file_type (str): Typ souboru ('GTF' nebo 'GTF'/'GFF3').

    Returns:
        str nebo None: Název genu nebo None, pokud není nalezen.
    """
    if file_type == 'GTF':
        # GTF format: gene_id "gene_name"; ...
        for attribute in attributes.split(';'):
            if attribute.strip().startswith('gene_id'):
                return attribute.split('"')[1]
    else:
        # GFF3 format: ID=gene0;Name=gene_name;...
        for attribute in attributes.split(';'):
            key, _, value = attribute.partition('=')
            if key == 'ID' or key == 'Name':
                return value
    return None


def parse_gtf_gff_for_strand(file_path):
    """
    Zpracovává soubor GTF nebo GFF3 a vrací slovník, kde jsou klíče názvy/ID genů
    a hodnoty informace o vláknech ('+' nebo '-').

    Args:
        file_path (str): Cesta k souboru GTF nebo GFF3.

    Returns:
        strand_info (dict): Slovník s názvy/ID genů jako klíči a informacemi o vláknech jako hodnotami.
    """
    strand_info = {}
    file_type = 'GTF' if file_path.endswith('.gtf') else 'GFF3'
    with open(file_path, 'r') as file:
        for line in file:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue
            gene_info = fields[8]
            strand = fields[6]
            gene_name = extract_gene_name(gene_info, file_type)
            if gene_name:
                strand_info[gene_name] = strand
    return strand_info


def reverse_complement(seq):
    """
    Vrací reverzní komplementární sekvenci pro zadanou DNA sekvenci.

    Args:
        seq (str): DNA sekvence, která má být převedena na její reverzní komplement.

    Returns:
        str: Reverzní komplementární sekvence k zadané DNA sekvenci.
    """
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return ''.join(complement[base] for base in reversed(seq))


def simulate_reads(seqs, read_length, num_reads, gc_content, operon_locations, strand_info):
    """
    Simuluje čtení na základě poskytnutých sekvencí.

    Args:
        seqs (dict): Slovník sekvencí, kde klíče jsou názvy sekvencí a hodnoty jsou sekvence.
        read_length (int): Délka čtení.
        num_reads (int): Počet čtení, které mají být simulovány pro každou sekvenci.
        gc_content (float nebo None): Obsah GC. Pokud není poskytnut, bude vypočítán.
        operon_locations (dict nebo None): Slovník s lokacemi operonů pro každou sekvenci.
        strand_info (dict nebo None): Informace o vláknech pro každou sekvenci.

    Returns:
        reads (list): Seznam čtení jako dvojic (název_sekvence, čtení).
        """
    reads = []

    # prochazeni kazde sekvence a generovani cteni
    for seq_name, seq in seqs.items():
        strand = strand_info.get(seq_name, "both") if strand_info else "both"

        # pokud neni obsah GC poskytnut, vypocita se pomoci funkce calculate_gc_content
        if gc_content is None:
            gc_content = calculate_gc_content(seq)

        if operon_locations:
            operons = operon_locations.get(seq_name, [])

            # generovani zadaneho poctu cteni pro aktualni sekvenci
            for _ in range(num_reads):
                # nahodny vyber pocatecni pozice pro cteni
                start_pos = random.randint(0, len(seq) - read_length)
                read_gc_content = gc_content

                # kontrola, zda je cteni v oblasti nejakeho operonu, a nastaveni prislusneho GC obsahu
                for start, end in operons:
                    if start <= start_pos < end:
                        read_gc_content = gc_content[(start, end)]
                        break

                read = ""
                for _ in range(read_length):
                    if random.random() < gc_content:
                        read += random.choice(["G", "C"])
                    else:
                        read += random.choice(["A", "T"])
                if strand == "reversed" or (strand == "both" and random.choice([True, False])):
                    read = reverse_complement(read)
                reads.append((seq_name, read))

        else:  # jinak simulovano bez zohledneni operonu
            for _ in range(num_reads):
                read = ""
                for _ in range(read_length):
                    if random.random() < gc_content:
                        read += random.choice(["G", "C"])
                    else:
                        read += random.choice(["A", "T"])
                if strand == "reversed" or (strand == "both" and random.choice([True, False])):
                    read = reverse_complement(read)
                reads.append((seq_name, read))
    return reads


def phred_to_ascii(phred_score):
    """
    Převod Phred score do hodnot ASCII.

    Args:
        phred_score (int): Phred quality score (běžně v rozmezí 0 až 40).

    Returns:
        str: Odpovídající hodnota ASCII.
    """
    return chr(phred_score + 33)


def generate_quality_scores(length, max_quality, min_quality):
    """
    Generuje klesající kvalitní skóre pro danou délku sekvence.
    Pokud nejsou zadány hraniční hodnoty, jsou odhadnuty na základě platformy Illumina.

    Args:
        length (int): Délka sekvence.
        max_quality (int, optional): Maximální ASCII hodnota pro kvalitní skóre na začátku sekvence.
        min_quality (int, optional): Minimální ASCII hodnota pro kvalitní skóre na konci sekvence.

    Returns:
        str: Řetězec klesajících kvalitních skóre.
    """
    # vychozi hodnoty ASCII na zaklade typickych skore kvality Illumina
    if max_quality is None:
        # pokud neni max_quality poskytnuto, nastavi se na vychozi skore vysoke kvality v ASCII
        # 40 (Phred skore) + 33 (posun pro prevod do ASCII)
        max_quality = 40 + 33
    if min_quality is None:
        # pokud není min_quality poskytnuto, nastavi se na vychozi skore nizsi kvality v ASCII
        # 30 (Phred skore) + 33 (posun pro prevod do ASCII)
        min_quality = 30 + 33

    quality_scores = []
    for i in range(length):
        quality = int(max_quality - ((max_quality - min_quality) * (i / length)))
        # nahodna variace
        quality = max(min_quality, min(max_quality, quality + random.randint(-5, 5)))
        quality_scores.append(chr(quality))

    return ''.join(quality_scores)


def write_output(reads, output_format, output_file, max_quality, min_quality):
    """
    Zapíše sekvence do souboru ve formátu FASTA nebo FASTQ.

    Args:
        reads (list of tuples): Seznam sekvencí a jejich názvů.
                                Každý prvek seznamu je tuple (seq_name, read),
                                kde 'seq_name' je název sekvence a 'read' je samotná sekvence.
        output_format (str): Formát výstupního souboru. Může být 'fasta' nebo 'fastq'.
        output_file (str): Název výstupního souboru.
        max_quality (int, optional): Maximální ASCII hodnota pro kvalitní skóre na začátku sekvence (pouze pro FASTQ).
        min_quality (int, optional): Minimální ASCII hodnota pro kvalitní skóre na konci sekvence (pouze pro FASTQ).

    Returns:
        Funkce vytvoří soubor s daným názvem a zapíše do něj sekvence.

    Pro FASTQ formát je možné nastavit maximální a minimální kvalitu čtení,
    přičemž pokud tyto hodnoty nejsou zadány, jsou odhadnuty na základě platformy Illumina.
    """
    with open(output_file, 'w') as file:
        for i, (seq_name, read) in enumerate(reads):
            if output_format.lower() == 'fasta':
                file.write(f'>{seq_name}_read{i}\n{read}\n')
            elif output_format.lower() == 'fastq':
                quality = generate_quality_scores(len(read), max_quality, min_quality)
                file.write(f'@{seq_name}_read{i}\n{read}\n+\n{quality}\n')


# Priklad simulace:

seqs = read_fasta('your_fasta_file')

operon_locations = extract_operons_from_gff('path_to_your_gff_file')

strand_info = parse_gtf_gff_for_strand('path_to_your_gff_file')

reads = simulate_reads(seqs, 75, 100, gc_content=None, operon_locations=operon_locations, strand_info=strand_info)

rRNA_percentage = calculate_rrna_percentage('your_fasta_file', 'path_to_your_gff_file', read_fasta)
print(f"rRNA Percentage: {rRNA_percentage:.2f}%")

write_output(reads, output_format='fastq', output_file='output.fastq', max_quality=None, min_quality=None)
