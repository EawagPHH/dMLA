from itertools import groupby
from collections import defaultdict
from dataclasses import dataclass
import gzip, glob


def read_fasta(fasta):
    open_func = gzip.open if fasta.endswith('.gz') else open
    with open_func(fasta, 'rt') as file_handle:
        grouped = groupby(file_handle, lambda x: x[0] == ">")
        for cond, entry in grouped:
            if cond:
                fasta_id = next(entry).strip()
                _, seq_iter = next(grouped)
                seq = ''.join([line.strip() for line in seq_iter]).upper()
                yield ([fasta_id, seq])


def reverse_complement(sequence):
    intab = "ATGCYRSWKMBDHVN-"
    outtab = "TACGRYSWMKVHDBN-"
    trans_table = str.maketrans(intab, outtab)
    complement = sequence.translate(trans_table)
    rev_complement = complement[::-1]
    return rev_complement


annot_dict = {"GTGAATTTGGCAAACTTTTCATCATATCGCTTAAAATACG": "mcr28",
              "TTGTACTTATCTTCGTAGGTTTTCAACTTGGCAATCATCT": "mcr64",
              "CGCTTAGCACCTCTGATAGTTGGTTCGAAATTTCGATGGT": "aadA1",
              "TTCCATAGCGTTAAGGTTTCATTTAGCGCCTCAAATAGAT": "aadA21",
              "CATCGATCGCAGATCCGAACAGGTGGATTGTGTCCAGTGT": "aadA57",
              "GCGGAGCGTAACGCGGCAACGACCGTCTCCGCTCCTCCTT": "aac3II63",
              "CACCAGCCACGCCGCATCGGGCCCGATCGCCGCCATCGAT": "aac3VI23",
              "AATGCCTGGCGTGTTTGAACCATGTACACGGCTGGACCAT": "aac6I6",
              "TTGTTTTAGGGCGACTGCCCTGCTGCGTAACATCGTTGCT": "aac6I123",
              "GGAGCCTCCGCGATTTCATACGCTTCGTCTGCCCACCAAG": "ant2Ia1",
              "TCCGCTGCAATTTCGGTCGCCTGGTAGTCGCCGTGCTCGG": "aph6I29",
              "AATATTTTCACCTGAATCAGGATATTCTTCTAATACCTGG": "aph3386",
              "GACGTTCCAGGAGCAACGATGCCTGGTAGCTGTCCAGTTT": "qnrB7",
              "CCGATAAATTCAGTGCCGCTCAGGTCGGCACCTGAAAAAT": "qnrB35",
              "AGGTGAGATCACTTAAGTCTTTATGTGAAAAGTTGTGGTG": "qnrS1",
              "TTTGTTAGGGAATTGAAACTGTAGAATATCTTGGTGAATT": "ermB1",
              "CGGAGAAATCTGGCCACGACGAATCGTCGTCGAGCCAGCG": "mphA18",
              "GGTAAATTTCTCCGCCACCAGACACTATAACGTGACCGGT": "dfrA243",
              "TTGATAGCTCTTTCAAAGCATTTTCTATTGAAGGAAAAAC": "dfrA369",
              "CAAAAGTCTTGCGTCCAACCAACAGCCATTGGTTATAGGT": "dfrA15",
              "CCCGCCACCAGACACTATAACGTGATCGGTGAGTTCAGCC": "dfrA58",
              "GAAAGACTAATACATTTTCATTTGAGCTTGAAATTCCTTT": "dfrA1",
              "GCGGCGCAATACGTCTGATCTCATCGGCCGGCGATACAGG": "sul15",
              "GAATGCATAACGACGAGTTTGGCAGATGATTTCGCCAATT": "sul22",
              "GCCTAAAAAGAAGCCCATACCCGGATCAAGAATAATTCGT": "sul336",
              "AGCGCCGCCAGTGAGCCTTGCAGCTGCCCCTGACGTTCCT": "tetA13",
              "AGCCAAGGAGCAAAGATAACCTGCATTAACGCATAAAGTG": "tetB43",
              "AAAACTAAGATATGGCTCTAACAATTCTGTTCCAGCTTTT": "tetM34",
              "GGTTACGACATTCGAACCGTGCAGGATCTGCTCGGCCATT": "intI125",
              "CGAAGCGGGTAGATATCACACTCTGTCTGGCTTTTGGCTG": "uidA160",
              "ACAACGAACTGAACTGGCAGACTATCCCGCCGGGAATGGT": "uidA616",
              "ACGTTGTTGCTGAAGCAACCGGTATCTTCCTGACCGACGA": "gadph195",
              "GTTGCGCACCATATGCACGATGCATAACTGGATGCGGGCC": "astA13",
              "ATTTTTGCTTCATAAGCCGATAGAAGATTATAGGATAATG": "aatA22",
              "AAGTCCATTTCTGCCAGTGTAATTCCCTCCGACAGTCAGT": "aatA12314",
              "GCAGTGCGGAGGTCATTTGCTGTCACTCCCGACACGCCAT": "ipaH321",
              "GCAATACCTCCGGATTCCGTGAACAGGTCGCTGCATGGCT": "ipaH921",
              "CAGTACCATTAATATCATGCGGAATATTCAGAGAAAGAAT": "eae111",
              "AAATCCTGATCAATGAAGACGTTATAGCCCAGCATATTTT": "eae6271",
              "ACTCCATTAACGCCAGATATGATGAAACCAGTGAGTGACG": "stx2181",
              "GAGCACTTTGCAGTAACGGTTGCAGATTCCAGCGACTGGT": "stx23291",
              "TATTGCTTCCTTCTTTTGTGTATACATCGATTGTACAATT": "aggR11",
              "GATCAACATCTTCAGCAGTCATTACATAAGAACGCCCACT": "stx1301",
              "TAAAGGTATCGTCATCATTATATTTTGTATACTCCACCTT": "stx11581",
              "GCCCAATATACAGACCATTAATTGCAGACGTTGCGCTCAT": "bfpA82",
              "TGAAGAATATTGGCAGGAAGAATTGTTAATGGCGTTACGT": "invE132",
              "TAACGCAGAAACCTCCTGTTCATATGGGTGAGGGCTGTAT": "eltA73",
              "GTCAGATATGATGACGGATATGTTTCCACTTCTCTTAGTT": "eltA4510",
              "GTTTTATTATTCCATACACATAATTTATCAATTTTGGTCT": "eltB12",
              "CAGAACTATGTTCGGAATATCGCAACACACAAATATATAC": "eltB2318",
              "GTGTTGTTCATATTTTCTGATTTTTTTTCACTGTTGTTTT": "estIa94",
              "GTAAAATGTGTTATTCATATTTTCTGATTTTTTTTCACTG": "estIa1515",
              "GATACCTTTCGCGTCAGCGCCGTAACTTCCTTGTGCAGGT": "bfpF56",
              "ATATCTTGGATTATAACGCATCGCTGAAATCATTGCATTT": "bfpF105",
              "ACTTGATAATAAAGTAGCAATAACAAACATTCTGATTTTT": "aafA83",
              "GACAGACACCCCCCTTACTTTCCCTTCATTAAGGCCTTGC": "aap223",
              "AAAGTAGTCACCACTGTTTTCAACTTTTGAATTTCCCTTT": "aaiC17",
              "TTTTTAAGAGAGGTGAAAAAGAAGTTAAAATAGAAATTTT": "aaiC47",
              "CAAAAGATAAAAATGGTTTTTCTGGATTATTTTTAGCGAT": "arpA19",
              "CTCACCGCCATATGTCAGTAAGTGAGAAGCGAAACTGTCG": "chuA24",
              "CAAAGAAGGCGCATTCGTTCCTTTCGTCACCCTCGGTGAT": "trpA354",
              "GCGTTTCTGTCTCACCCGCAAGGACAGCGCTGGCGATATA": "Tsp19",
              "AACCGCCCTCATTCATTTCCAGCGCGCTCACAACAATATC": "yjaA12"}


bc_dict = {"CCACTAG": "450",
           "GTCTTCT": "451",
           "ACAAAGC": "452",
           "AAAAGGC": "453",
           "ACTGTGT": "454",
           "CTGGTAC": "455",
           "GTTGCTA": "456",
           "CCATTCA": "457",
           "GATATCG": "458",
           "ACAGGAT": "459",
           "CGCATAC": "460",
           "GGACCTA": "461",
           "GAACTGA": "462",
           "AGTCGTG": "463",
           "TTCCAGG": "464",
           "TAGAGCG": "465",
           "AAACCCT": "466",
           "ATCCAGT": "467",
           "TCTTCGT": "468",
           "GTGTCCT": "469",
           "CACAGAT": "470",
           "ATAGAGC": "471",
           "GGTCATG": "472",
           "CCATCTC": "473",
           "TATGCAG": "474",
           "AACCGCT": "475",
           "TGACTCT": "476",
           "AAGCACA": "477",
           "TGTGAAC": "478",
           "AACAGTG": "479",
           "GCCTATA": "480",
           "GTCCTCT": "481",
           "TAGAACG": "482",
           "TGGGTGA": "483",
           "ACACGTG": "484",
           "TGCCCAA": "485",
           "TCGTACA": "486",
           "AACCAAG": "487",
           "GTGCTAT": "488",
           "TCCTCAT": "489",
           "GTCGTAT": "490",
           "TCTGGAA": "491",
           "CAGTAGG": "492",
           "ACGTCAT": "493",
           "ACGGAAC": "494",
           "CGGATGA": "495",
           "CCGCATA": "496",
           "GTGAAGA": "497",
           "AACGTGT": "498",
           "AGAAGAC": "499",
           "TGTAGGG": "500",
           "AATGCCT": "501",
           "GCTTTCT": "502",
           "CATGAAG": "503",
           "CGATCAC": "504",
           "TCGCGTT": "505",
           "GATTGGC": "506",
           "CATCGGT": "507",
           "CTTTCCA": "508",
           "ACCAGAT": "509",
           "GTGATAC": "510",
           "TGCATCC": "511",
           "TGCAAAC": "512",
           "GCAACGA": "513",
           "GCCATAC": "514",
           "TACCTTC": "515",
           "GAATCGA": "516",
           "CTACGTT": "517",
           "GTTTCGG": "518",
           "GCAAATG": "519",
           "CTGACAC": "520",
           "GCACTTT": "521"}


@dataclass
class AnnotatedFasta:
    fa_id: str
    seq: str
    annot: str = "NA"
    rc: bool = False
    target: str = "NA"
    umi: str = "NA"
    bc: str = "NA"
    sample: str = "NA"
    def __post_init__(self):
        self.process_sequence(self.seq)
    def process_sequence(self, seq: str):
        try:
            self.annotate_sequence(seq)
        except ValueError:
            pass
        try:
            rev_seq = reverse_complement(seq)
            self.annotate_sequence(rev_seq, reverse=True)
        except ValueError:
            pass
    def annotate_sequence(self, seq: str, reverse: bool = False):
        _, trunc_seq = seq.split("TCTTTTCGCAGGCTGGAGCCCAGGTCTTCCTAT")
        trunc_seq, _ = trunc_seq.split("GAATGAGTGTGCGTGCACTC")
        self.bc, self.target = trunc_seq.split("TGGGCCCAATTTTCCGTGACAATTAATT")
        self.rc = reverse
        self.seq = seq
        self.umi = self.target[-10:]
        self.target = self.target[:-10]
        self.annot = annot_dict.get(reverse_complement(self.target), "NA")
        self.sample = bc_dict.get(self.bc, "NA")


def process_file(input_file, output_file):
    f_file = read_fasta(input_file)
    acc = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))
    for f_id, f_line in f_file:
        _, count = f_id.split("=")
        fa = AnnotatedFasta(f_id, f_line)
        acc[fa.sample][fa.annot][fa.umi] += int(count)
    with open(output_file, "w") as output:
        for sample, rest1 in acc.items():
            for target, rest2 in rest1.items():
                if sample != "NA" and target != "NA":
                    output.write(f"{sample},{target},{len(rest2)}\n")


for file in glob.glob("fasta/*.fa.gz"):
    trunc_file = file.split(".")[0].split("/")[1]
    output_file = f"{trunc_file}.csv"
    print(file)
    process_file(file, output_file)