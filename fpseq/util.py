import unicodedata
import re


AA_WEIGHTS = {
    'A': 89.0932,
    'C': 121.1582,
    'D': 133.1027,
    'E': 147.1293,
    'F': 165.1891,
    'G': 75.0666,
    'H': 155.1546,
    'I': 131.1729,
    'K': 146.1876,
    'L': 131.1729,
    'M': 149.2113,
    'N': 132.1179,
    'O': 255.3134,
    'P': 115.1305,
    'Q': 146.1445,
    'R': 174.201,
    'S': 105.0926,
    'T': 119.1192,
    'U': 168.0532,
    'V': 117.1463,
    'W': 204.2252,
    'Y': 181.1885
}
WATER = 18.0153

letters_3to1 = {
    'Ala': 'A',
    'Arg': 'R',
    'Asn': 'N',
    'Asp': 'D',
    'Asx': 'B',
    'Cys': 'C',
    'Gln': 'Q',
    'Glu': 'E',
    'Glx': 'Z',
    'Gly': 'G',
    'His': 'H',
    'Ile': 'I',
    'Leu': 'L',
    'Lys': 'K',
    'Met': 'M',
    'Phe': 'F',
    'Pro': 'P',
    'Pyl': 'O',
    'Sec': 'U',
    'Ser': 'S',
    'Thr': 'T',
    'Trp': 'W',
    'Tyr': 'Y',
    'Val': 'V',
    'Xaa': 'X',
    'Xle': 'J'
}

letters_1to3 = {
    'A': 'Ala',
    'B': 'Asx',
    'C': 'Cys',
    'D': 'Asp',
    'E': 'Glu',
    'F': 'Phe',
    'G': 'Gly',
    'H': 'His',
    'I': 'Ile',
    'J': 'Xle',
    'K': 'Lys',
    'L': 'Leu',
    'M': 'Met',
    'N': 'Asn',
    'O': 'Pyl',
    'P': 'Pro',
    'Q': 'Gln',
    'R': 'Arg',
    'S': 'Ser',
    'T': 'Thr',
    'U': 'Sec',
    'V': 'Val',
    'W': 'Trp',
    'X': 'Xaa',
    'Y': 'Tyr',
    'Z': 'Glx'
}


def seq1(seq):
    """Convert protein sequence from three-letter to one-letter code. """
    return ''.join(letters_3to1.get(aa) for aa in seq)


def seq3(seq):
    """Convert protein sequence from one-letter to three-letter code. """
    return ''.join(letters_1to3.get(aa) for aa in seq)


def protein_weight(seq):
    try:
        return sum(AA_WEIGHTS[x] for x in seq) - (len(seq) - 1) * WATER
    except KeyError as e:
        raise ValueError('%s is not a valid unambiguous amino acid letter' % e)


def slugify(value, allow_unicode=False):
    """
    Convert to ASCII if 'allow_unicode' is False. Convert spaces to hyphens.
    Remove characters that aren't alphanumerics, underscores, or hyphens.
    Convert to lowercase. Also strip leading and trailing whitespace.
    """
    value = str(value)
    if allow_unicode:
        value = unicodedata.normalize('NFKC', value)
    else:
        value = unicodedata.normalize('NFKD', value).encode('ascii', 'ignore').decode('ascii')
    value = re.sub(r'[^\w\s-]', '', value).strip().lower()
    return re.sub(r'[-\s]+', '-', value)
