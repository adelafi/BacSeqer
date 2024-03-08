## BacSeqer

# Author: Adela Fialova
# Python 3.10.


import random
import pandas as pd
import numpy as np
from Bio import SeqIO


def read_fasta(file_path):
    """
    Načte FASTA soubor a vrátí sekvence.

    Args:
    file_path (str): Cesta k FASTA souboru.

    Returns:
    seqs: Slovník, kde klíč je název sekvence a hodnota je samotná sekvence.
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


def read_gtf_for_fasta(gtf_path):
    """
    Načte GTF soubor a vrátí lokace pro extrakci sekvencí.

    Args:
    gtf_path (str): Cesta k GTF souboru.

    Returns:
    dict: Slovník, kde klíč je název sekvence a hodnota je seznam tuple (start, end).
    """
    seq_locations = {}                                  # vytvoreni slovniku pro ukladani seznamu lokaci
    with open(gtf_path, 'r') as file:                   # otevreni GTF souboru
        for line in file:                               # prochazeni souboru po radcich
            if line.startswith('#'):                    # hledani zacatku komentare
                continue                                # preskoceni komentare
            parts = line.strip().split('\t')            # odstraneni znaku na zacatku a na konci radku, rozdeleni radku na casti pomoci tabulatoru
            seq_name = parts[0]                         # ulozeni nazvu sekvence
            start = int(parts[3])                       # prevod zacatku lokace na cele cislo
            end = int(parts[4])                         # prevod konce lokace na cele cislo

            if seq_name not in seq_locations:           # inicializace prazdneho seznamu pro nazev sekvence, pokud jeste neni ve slovniku
                seq_locations[seq_name] = []
            seq_locations[seq_name].append((start, end))# pridani zacatku a konce do seznamu lokaci pro danou sekvenci
    return seq_locations                                # vraceni slovniku lokaci


def extract_sequence(sequence, locs):
    """
    Extrahuje a spojí podsekvence na základě lokací.

    Args:
    sequence (str): Celá sekvence.
    locs (list of tuple): Seznam lokací (start, end).

    Returns:
    str: Spojená podsekvence.
    """
    extracted = ''                                      # vytvoreni prazdneho retezce pro pridavani podsekvenci
    for start, end in locs:                             # prochazeni seznamu lokaci
        extracted += sequence[start-1:end]              # extrakce podsekvence (-1 protoze GTF je 1-based)
    return extracted                                    # vraceni spojene extrahovane sekvence


def extract_sequences_from_fasta(fasta_path, locations):
    """
    Extrahuje sekvence z FASTA souboru na základě lokací z GTF.

    Args:
    fasta_path (str): Cesta k FASTA souboru.
    locations (dict): Slovník lokací získaných z GTF.

    Returns:
    dict: Slovník, kde klíč je název sekvence a hodnota je sekvence.
    """
    seqs = {}                                           # vytvoreni prazdneho slovniku pro ukladani sekvenci
    with open(fasta_path, 'r') as file:                 # otevreni FASTA souboru
        current_seq = ""                                # prazdny retezec pro aktualne ctenou sekvenci
        current_seq_name = ""                           # prazdny retezec pro nazev aktualne ctene sekvence
        for line in file:                               # prochazeni souboru radek po radku
            if line.startswith('>'):                    # hledani hlavicky
                if current_seq_name and current_seq_name in locations: # kontrola spravne nactene sekvence
                    seqs[current_seq_name] = extract_sequence(current_seq, locations[current_seq_name]) # ulozeni extrahovane sekvence do slovniku
                current_seq_name = line.strip().split()[0][1:] # nastaveni nazvu nove sekvence
                current_seq = ""                        # resetovani aktualni sekvence
            else:
                current_seq += line.strip()             # pridavani retezcu sekvence do current_seq
        if current_seq_name and current_seq_name in locations: # pro posledni sekvenci
            seqs[current_seq_name] = extract_sequence(current_seq, locations[current_seq_name])
    return seqs                                         # vraceni slovniku s extrahovanymi sekvencemi

# Příklad použití
#locations = read_gtf_for_fasta('GTFtry.gtf')
#sequences = extract_sequences_from_fasta('FASTAtry.fna', locations)


def assign_expression_to_transcripts(seqs, mu):
    """
    Přiřadí každé sekvenci v 'seqs' náhodnou expresní úroveň z geometrické distribuce.

    Args:
    seqs (dict): Slovník sekvencí s klíčem jako názvem sekvence a hodnotou jako sekvencí.
    mu (float): Průměrný počet neúspěchů před prvním úspěchem pro geometrickou distribuci.

    Returns:
    dict: Slovník s klíčem jako názvem sekvence a hodnotou jako expresní úroveň.
    """
    expression_levels = {}
    for seq_name in seqs.keys():
        # Geometrická distribuce očekává pravděpodobnost 'p' jako parametr
        p = 1 / (mu + 1)
        expression_level = np.random.geometric(p)
        expression_levels[seq_name] = expression_level
    return expression_levels


def simulate_reads(seqs, read_length, num_reads, gc_content):
    reads = []
    for seq_name, seq in seqs.items():
        for _ in range(num_reads):
            read = ""
            for _ in range(read_length):
                # Zde je zavedena pravděpodobnost pro výběr GC nebo AT s ohledem na celkový GC obsah
                if random.random() < gc_content:
                    read += random.choice(["G", "C"])
                else:
                    read += random.choice(["A", "T"])
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

# Example usage:
# score = 32
# print(f'Phred score {score} -> ASCII character: {phred_to_ascii(score)}')


def generate_quality_scores(length, max_quality=65, min_quality=61):
    """
    Generuje klesající kvalitní skóre pro danou délku sekvence.

    Args:
    length (int): Délka sekvence.
    max_quality (int): Maximální ASCII hodnota pro kvalitní skóre na začátku sekvence.
    min_quality (int): Minimální ASCII hodnota pro kvalitní skóre na konci sekvence.

    Returns:
    str: Řetězec klesajících kvalitních skóre.
    """
    quality_scores = []
    for i in range(length):
        # Lineární interpolace mezi max_quality a min_quality
        quality = int(max_quality - ((max_quality - min_quality) * (i / length)))
        # Přidání náhodné variace
        quality = max(min_quality, min(max_quality, quality + random.randint(-5, 5)))
        quality_scores.append(chr(quality))

    return ''.join(quality_scores)


def write_output(reads, output_format, output_file):
    """
    Zapíše sekvence do souboru ve formátu FASTA nebo FASTQ.

    Args:
    reads (list of tuples): Seznam sekvencí a jejich názvů.
                            Každý prvek seznamu je tuple (seq_name, read),
                            kde 'seq_name' je název sekvence a 'read' je samotná sekvence.
    output_format (str): Formát výstupního souboru. Může být 'fasta' nebo 'fastq'.
    output_file (str): Cesta k výstupnímu souboru.

    Returns:
    Funkce vytvoří soubor s daným názvem a zapíše do něj sekvence.

    Tato funkce nyní generuje náhodné kvalitní skóre pro simulovaná čtení v FASTQ formátu.
    """
    with open(output_file, 'w') as file:
        for i, (seq_name, read) in enumerate(reads):
            if output_format.lower() == 'fasta':
                file.write(f'>{seq_name}_read{i}\n{read}\n')
            elif output_format.lower() == 'fastq':
                quality = generate_quality_scores(len(read))
                file.write(f'@{seq_name}_read{i}\n{read}\n+\n{quality}\n')

# Příklad použití
seqs = read_fasta('moje.fna')
reads = simulate_reads(seqs, read_length=75, num_reads=1000, gc_content=0.44)
write_output(reads, output_format='fastq', output_file='output.fastq')