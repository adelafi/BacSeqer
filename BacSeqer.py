# BacSeqer - Simulátor čtení pro bakteriální RNA-Seq

# Autor: Adela Fialova (adelafialova21@gmail.com)
# Python 3.10


import random
import re
from Bio import SeqIO


def read_fasta(file_path):
    """
    Načte FASTA soubor a vrátí sekvence.

    Args:
        file_path (str): Cesta k FASTA souboru.

    Returns:
        seqs (dict): Slovník, kde klíč je název sekvence a hodnota je samotná sekvence.
    """
    with open(file_path, 'r') as file:
        seqs = {}
        current_seq = ""
        for line in file:
            if line.startswith('>'):
                if current_seq:
                    seqs[seq_name] = current_seq
                    current_seq = ""
                seq_name = line.strip()
            else:
                current_seq += line.strip()
        if current_seq:
            seqs[seq_name] = current_seq
    return seqs


def extract_sequences(fasta_file, gtf_gff_file):
    """
    Extrakce sekvencí ze souboru FASTA na základě anotací v souboru GTF/GFF.

    Args:
        fasta_file (str): Cesta k souboru ve formátu FASTA.
        gtf_gff_file (str): Cesta k souboru GTF nebo GFF s anotacemi.

    Returns:
        extracted_seqs (dict): Slovník, kde klíče jsou ID anotací a hodnoty jsou extrahované sekvence.
    """

    def read_fasta_II(fasta_file):
        """
        Načte soubor ve formátu FASTA a vrátí slovník sekvencí.

        Args:
            fasta_file (str): Cesta k souboru FASTA.

        Returns:
            seqs (dict): Slovník, kde klíče jsou identifikátory záznamů a hodnoty jsou odpovídající sekvence.
        """

        seqs = {}
        with open(fasta_file, 'r') as file:
            for record in SeqIO.parse(file, "fasta"):
                seqs[record.id] = str(record.seq)
        return seqs

    def parse_gtf_gff(gtf_gff_file):
        """
        Zpracování souboru GTF/GFF a vrácení slovníku všech umístění prvků.

        Args:
            gtf_gff_file (str): Cesta k souboru GTF/GFF.

        Returns:
            feature_locations (dict): Slovník obsahující umístění všech prvků, kde klíčem je ID.
                                      Každá hodnota je trojice obsahující chromozóm, začátek a konec prvku.
        """

        def extract_id_from_attributes(attributes):
            """Extrahuje ID."""
            for attr in attributes.split(';'):
                if attr.strip().startswith('ID='):
                    return attr.split('=')[1].strip()
            return None

        feature_locations = {}
        with open(gtf_gff_file, 'r') as file:
            for line in file:
                if line.startswith('#') or line.strip() == '':
                    continue
                parts = line.strip().split('\t')
                if len(parts) < 9:
                    continue
                chrom, source, feature_type, start, end, score, strand, phase, attributes = parts
                feature_id = extract_id_from_attributes(attributes)
                if feature_id:
                    feature_locations[feature_id] = (chrom, int(start), int(end))
        return feature_locations

    seqs = read_fasta_II(fasta_file)
    annotations = parse_gtf_gff(gtf_gff_file)
    extracted_seqs = {}
    for annotation_id, (chrom, start, end) in annotations.items():
        if chrom in seqs:
            extracted_seq = seqs[chrom][start - 1:end]
            extracted_seqs[annotation_id] = extracted_seq
    return extracted_seqs


def calculate_gc_content(seq):
    """
    Spočítá GC content pro načtené sekvence.

    Args:
        seq (str): Sekvence.

    Returns:
        float: %GC content načtené sekvence.
    """
    return (seq.count('G') + seq.count('C')) / len(seq)


def adjust_gc_content(read, gc_content):
    """
    Upraví obsah GC v čtení na požadovaný obsah GC při zachování jeho délky.

    Args:
        read (str): Sekvence, kterou chceme upravit.
        gc_content (float): Požadovaný obsah GC (v rozsahu 0 až 1).

    Returns:
        read (str): Upravená sekvence.
    """

    current_gc_content = (read.count('G') + read.count('C')) / len(read)

    # Spočítat rozdíl mezi aktuálním a požadovaným GC obsahem
    diff = gc_content - current_gc_content

    # Počet A nebo T bazí, které je třeba nahradit
    num_to_replace = int(abs(diff) * len(read))

    # Náhodně vybereme indexy A nebo T bazí, které nahradíme
    if diff > 0:
        at_indices = [i for i, base in enumerate(read) if base in ('A', 'T')]
        replace_indices = random.sample(at_indices, num_to_replace)

        # Nahradíme A nebo T na vybraných pozicích
        for idx in replace_indices:
            if read[idx] == 'A':
                read = read[:idx] + random.choice('GC') + read[idx + 1:]
            else:  # read[idx] == 'T'
                read = read[:idx] + random.choice('GC') + read[idx + 1:]
    elif diff < 0:
        gc_indices = [i for i, base in enumerate(read) if base in ('G', 'C')]
        replace_indices = random.sample(gc_indices, num_to_replace)

        # Nahradíme G nebo C na vybraných pozicích
        for idx in replace_indices:
            if read[idx] == 'G':
                read = read[:idx] + random.choice('AT') + read[idx + 1:]
            else:  # read[idx] == 'C'
                read = read[:idx] + random.choice('AT') + read[idx + 1:]
    elif diff == 0:
        pass

    return read


def extract_operons(file_path):
    """
    Extrahuje operony a jejich umístění ze souboru GTF nebo GFF.

    Args:
        file_path (str): Cesta k souboru GTF nebo GFF.

    Returns:
        operons (dict): Slovník s názvy operonů jako klíči a jejich umístěními jako hodnotami.
    """
    operons = {}
    with open(file_path, 'r') as file:
        for line in file:
            if line.startswith('#') or line.strip() == "":
                continue  # Přeskočení hlaviček a prázdných řádků

            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue  # Přeskočení neúplných řádků

            attributes = parts[8]

            # Zpracování GTF/GFF atributů
            note_match = None
            if 'Note=' in attributes:  # GFF formát
                note_match = re.search(r'Note=Operon\s+(\d+)', attributes)
            elif 'note "' in attributes:  # GTF formát
                note_match = re.search(r'note\s+"Operon\s+(\d+)"', attributes)

            if note_match:
                operon_number = note_match.group(1)
                operon_name = f"Operon {operon_number}"
                location = (parts[0], int(parts[3]), int(parts[4]))  # (seqname, začátek, konec)

                if operon_name in operons:
                    operons[operon_name].append(location)
                else:
                    operons[operon_name] = [location]

    return operons


def extract_rrna(file_path):
    """
    Zpracovává soubor GFF nebo GTF za účelem extrakce anotací rRNA. Vrátí slovník, kde klíčem je seq_id s indexem a hodnotami
    jsou start a end anotace rRNA.

    Args:
        file_path (str): Cesta k souboru GFF nebo GTF.

    Returns:
        rrna_dict (dict): Slovník, kde klíče jsou seq_id s indexem a hodnoty jsou slovníky s klíči 'start' a 'end'.
    """
    rrna_dict = {}
    rrna_index = {}

    with open(file_path, 'r') as file:
        for line in file:
            if line.startswith('#') or line.strip() == '':
                continue  # Přeskočit komentáře a prázdné řádky

            fields = line.strip().split('\t')
            if len(fields) == 9:
                seq_id, source, feature_type, start, end, score, strand, phase, attributes = fields

                if feature_type == 'rRNA':
                    # Rozlišení formátu na základě obsahu atributů
                    if 'gene_id "' in attributes:  # GTF formát
                        gene_id = next((item.split('"')[1] for item in attributes.split(';') if 'gene_id "' in item),
                                       None)
                    elif 'ID=' in attributes:  # GFF formát
                        gene_id = next((item.split('=')[1] for item in attributes.split(';') if item.startswith('ID=')),
                                       None)

                    if gene_id:
                        index = rrna_index.get(gene_id, 0) + 1
                        rrna_dict[f"{gene_id}_{index}"] = {'start': int(start), 'end': int(end)}
                        rrna_index[gene_id] = index

    return rrna_dict


def extract_cds(file_path, start_extension, end_extension):
    """
    Extrakce CDS (Coding DNA Sequences) z GTF/GFF souboru.

    Args:
        file_path (str): Cesta k GTF/GFF souboru.
        start_extension (int, optional): Počet nukleotidů přidaných na začátek CDS.
        end_extension (int, optional): Počet nukleotidů přidaných na konec CDS.

    Returns:
        cds_dict (dict): Slovník, kde klíče jsou identifikátory genů nebo transkriptů a hodnoty jsou seznamy tuplů
        (start, end) CDS.
    """
    cds_dict = {}

    with open(file_path, 'r') as file:
        for line in file:
            if line.startswith('#') or line.strip() == '':
                continue

            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue

            feature_type = fields[2]
            if feature_type == 'CDS':
                start, end = int(fields[3]), int(fields[4])
                start = max(1, start - start_extension)
                end += end_extension

                attributes = fields[8]
                gene_id = None

                # Pro GTF formát
                if "gene_id" in attributes:
                    gene_id = [attr for attr in attributes.split(';') if 'gene_id' in attr][0].split('"')[1]
                # Pro GFF3 formát
                elif "ID=" in attributes:
                    gene_id = [attr for attr in attributes.split(';') if 'ID=' in attr][0].split('=')[1]

                if gene_id:
                    if gene_id in cds_dict:
                        cds_dict[gene_id].append((start, end))
                    else:
                        cds_dict[gene_id] = [(start, end)]

    return cds_dict


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
    # Načtení sekvencí ze souboru FASTA
    sequences = read_fasta_func(fasta_path)
    total_sequence_length = sum(len(seq) for seq in sequences.values())

    # Extrahování anotací rRNA ze souboru GTF/GFF
    rrna_dict = extract_rrna(gff_path)

    # Výpočet celkové délky sekvencí rRNA
    rrna_total_length = 0
    for seq_id, rrna_info in rrna_dict.items():
        start, end = rrna_info['start'], rrna_info['end']
        rrna_total_length += end - start + 1

    # Výpočet a vrácení procentuálního podílu rRNA
    return (rrna_total_length / total_sequence_length) * 100


def parse_for_strand(file_path):
    """
    Zpracovává soubor GTF nebo GFF3 a vrací slovník, kde jsou klíče názvy/ID genů
    a hodnoty informace o vláknech ('+' nebo '-').

    Args:
        file_path (str): Cesta k souboru GTF nebo GFF3.

    Returns:
        strand_info (dict): Slovník s názvy/ID genů jako klíči a informacemi o vláknech jako hodnotami.
    """

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


def create_long_mrna(seqs, operon_locations):
    """
    Vytvoří dlouhou mRNA sekvenci na základě lokací operonů.

    Args:
        seqs (dict): Slovník sekvencí, kde klíče jsou identifikátory genů a hodnoty jsou DNA sekvence.
        operon_locations (dict): Slovník s lokacemi operonů, kde klíče jsou identifikátory operonů a hodnoty jsou
        seznamy identifikátorů genů v operonu.

    Returns:
        long_mrna_sequences (dict): Slovník s dlouhými mRNA sekvencemi pro každý operon.
    """
    long_mrna_sequences = {}

    for operon_id, gene_ids in operon_locations.items():
        long_mrna_sequence = ''
        for gene_id in gene_ids:
            gene_sequence = seqs.get(gene_id, '')
            if gene_sequence:
                long_mrna_sequence += gene_sequence

        if long_mrna_sequence:  # Zajištění, že výsledná mRNA sekvence není prázdná
            long_mrna_sequences[operon_id] = long_mrna_sequence

    return long_mrna_sequences


def reverse_complement(seq):
    """
    Vrací reverzní komplementární sekvenci pro zadanou DNA sekvenci.

    Args:
        seq (str): DNA sekvence, která má být převedena na její reverzní komplement.

    Returns:
        str: Reverzní komplementární sekvence k zadané DNA sekvenci.
    """
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', "Y": random.choice("ACGT")}
    return ''.join(complement[base] for base in reversed(seq))


def simulate_reads(seqs, read_length, num_reads, gc_content, cds_locations, operon_locations, rrna_locations,
                   rrna_percentage, strand_ori, strand_info):
    """
    Simuluje čtení na základě poskytnutých sekvencí.

    Args:
        seqs (dict): Slovník sekvencí, kde klíče jsou názvy sekvencí a hodnoty jsou sekvence.
        read_length (int): Délka čtení.
        num_reads (int): Počet čtení, která mají být simulována pro každou sekvenci.
        gc_content (float nebo None): Požadovaný obsah GC.
        cds_locations (dict nebo None): Slovník s lokacemi kódujících oblastí.
        operon_locations (dict nebo None): Slovník s lokacemi operonů.
        rrna_locations (dict nebo None): Slovník s lokacemi rRNA.
        rrna_percentage (float nebo None): Požadované procentuální zastoupení rRNA ve vstupní knihovně (0 až 100).
        strand_ori (str): Požadovaná strand orientation. (+ nebo -)
        strand_info (dict nebo None): Informace o orientaci vláken pro každou sekvenci.

    Returns:
        reads (list): Seznam čtení jako dvojic (název sekvence, čtení).
        """
    reads = []
    if rrna_percentage:
        rrna_seqs = seqs

    if cds_locations:
        for gene_id, locations in cds_locations.items():
            for start, end in locations:
                cds_seq = seqs.get(gene_id, '')[start:end]
                if cds_seq:  # Zajištění, že sekvence není prázdná
                    seqs[gene_id] = cds_seq

    if operon_locations:
        operons = create_long_mrna(seqs, operon_locations)
        to_add = len(operon_locations.keys())
        for i in range(0, to_add):
            updated_operons = {f"{key}_{i}": value for key, value in operons.items()}
            seqs.update(updated_operons)

    if rrna_percentage:
        num_dna = len(seqs)
        rrna_needed = int((num_dna * rrna_percentage) / (100 - rrna_percentage))
        rrnas = {}
        for i, id in enumerate(rrna_locations.keys()):
            id = id[:-2]
            rrna = rrna_seqs[id]
            new_key = f"rrna_{i}_{id}"
            rrnas[new_key] = rrna
        for i in range(0, rrna_needed):
            selected_key = random.choice(list(rrnas.keys()))
            seqs[selected_key] = rrnas[selected_key]

    while len(reads) != num_reads:
        seq_name = random.choice(list(seqs.keys()))
        seq = seqs[seq_name]

        # pokud neni obsah GC poskytnut, vypocita se pomoci funkce calculate_gc_content
        if gc_content is None:
            gc_content = calculate_gc_content(seq)

        # generovani zadaneho poctu cteni pro aktualni sekvenci
        if len(seq) >= read_length:
            read_start = random.randint(0, len(seq) - read_length)
            read = seq[read_start:read_start + read_length]
            read = adjust_gc_content(read, gc_content)
            if strand_ori:
                if strand_info[seq_name] != strand_ori:
                    read = reverse_complement(read)
            reads.append((seq_name, read))
        else:
            continue

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

    with open(output_file, 'w') as file:
        for i, (seq_name, read) in enumerate(reads):
            if output_format.lower() == 'fasta':
                file.write(f'>{seq_name}_read{i}\n{read}\n')
            elif output_format.lower() == 'fastq':
                quality = generate_quality_scores(len(read), max_quality, min_quality)
                file.write(f'@{seq_name}_read{i}\n{read}\n+\n{quality}\n')


