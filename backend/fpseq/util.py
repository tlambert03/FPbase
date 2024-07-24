import re
import unicodedata
from textwrap import wrap

AA_WEIGHTS = {
    "A": 89.0932,
    "C": 121.1582,
    "D": 133.1027,
    "E": 147.1293,
    "F": 165.1891,
    "G": 75.0666,
    "H": 155.1546,
    "I": 131.1729,
    "K": 146.1876,
    "L": 131.1729,
    "M": 149.2113,
    "N": 132.1179,
    "O": 255.3134,
    "P": 115.1305,
    "Q": 146.1445,
    "R": 174.201,
    "S": 105.0926,
    "T": 119.1192,
    "U": 168.0532,
    "V": 117.1463,
    "W": 204.2252,
    "Y": 181.1885,
}
WATER = 18.0153

letters_3to1 = {
    "Ala": "A",  # alanine
    "Arg": "R",  # arginine
    "Asn": "N",  # asparagine
    "Asp": "D",  # aspartic acid
    "Asx": "B",  # asparagine or aspartic acid
    "Cys": "C",  # cysteine
    "Gln": "Q",  # glutamine
    "Glu": "E",  # glutamic acid
    "Glx": "Z",  # glutamine or glutamic acid
    "Gly": "G",  # glycine
    "His": "H",  # histidine
    "Ile": "I",  # isoleucine
    "Leu": "L",  # leucine
    "Lys": "K",  # lysine
    "Met": "M",  # methionine
    "Phe": "F",  # phenylalanine
    "Pro": "P",  # proline
    "Pyl": "O",  # pyrrolysine (unusual)
    "Sec": "U",  # selenocysteine (unusual)
    "Ser": "S",  # serine
    "Thr": "T",  # threonine
    "Trp": "W",  # tryptophan
    "Tyr": "Y",  # tyrosine
    "Val": "V",  # valine
    "Xaa": "X",  # any amino acid
    "Xle": "J",  # leucine or isoleucine
    "TERM": "*",  # terminatino codon
}

letters_1to3 = {v: k for k, v in letters_3to1.items()}


def seq1(seq):
    """Convert protein sequence from three-letter to one-letter code."""
    return "".join(letters_3to1.get(aa) for aa in seq)


def seq3(seq):
    """Convert protein sequence from one-letter to three-letter code."""
    return "".join(letters_1to3.get(aa) for aa in seq)


def protein_weight(seq):
    try:
        return sum(AA_WEIGHTS[x] for x in seq) - (len(seq) - 1) * WATER
    except KeyError as e:
        raise ValueError(f"{e} is not a valid unambiguous amino acid letter") from e


def slugify(value, allow_unicode=False):
    """
    Convert to ASCII if 'allow_unicode' is False. Convert spaces to hyphens.
    Remove characters that aren't alphanumerics, underscores, or hyphens.
    Convert to lowercase. Also strip leading and trailing whitespace.
    """
    value = str(value)
    if allow_unicode:
        value = unicodedata.normalize("NFKC", value)
    else:
        value = unicodedata.normalize("NFKD", value).encode("ascii", "ignore").decode("ascii")
    value = re.sub(r"[^\w\s-]", "", value).strip().lower()
    return re.sub(r"[-\s]+", "-", value)


def evenchunks(string, chunksize=10):
    out = []
    for i in range(0, len(string), chunksize):
        end = i + chunksize
        out.append(string[i:end])
    return out


def chunked_lines(string, chunksize=10, chunks_per_line=5, spacer=" "):
    chunks = evenchunks(string, chunksize)
    lines = []
    while chunks:
        lines.append(spacer.join(chunks[:chunks_per_line]))
        del chunks[:chunks_per_line]
    return lines


def chunk_string(string, chunksize=10, max_width=55, space=" "):
    chunks = wrap(string, chunksize)
    line_chunks = []
    line = ""
    while chunks:
        cur_chunk = chunks.pop(0)
        if len(line) + len(cur_chunk) + 1 <= max_width:
            line += space if len(line) else ""
            line += cur_chunk
        else:
            line_chunks.append(line)
            line = cur_chunk
    line_chunks.append(line)
    return line_chunks
